module simulation_data_mod
  use working_prec_mod, only: ik, rk
  use constants_mod   , only: unity
  use util_mod        , only: trace
  use transition_matrix_mod, only: tm_t
  use wang_landau_mod , only: wl_t
  use wltm_mod        , only: wltm_t
  use sampling_mod    , only: sample_array_t, sample_scalar_t
  use cell_mod        , only: cell_t
  use lj_potential_mod, only: lj_potential_t
  implicit none
  private

  public :: simulation_data_t
  public :: simulation_sampling_t

  type :: simulation_data_t
    ! derived types
    type(cell_t)         :: cell  !< cell: volume, density, cut-off data, cellmatrix, coordinates
    type(lj_potential_t) :: pot2b !< lennard-jones potential
    type(tm_t)           :: tm    !< transition matrix
    type(wl_t)           :: wl    !< wang landau
    type(wltm_t)         :: wltm    !< wang landau Transition matrix

    ! thermodynamics
    character(len=3) :: ensemble !< ensemble
    real(rk) :: energy           !< system energy
    real(rk) :: running_energy   !< running energy
    real(rk) :: virial(3,3)      !< system virial
    real(rk) :: pressure         !< system pressure
    real(rk) :: pressure_in      !< constant pressure for NPT simulation
    real(rk) :: pressure_tensor(3,3)  !< system pressure tensor
    real(rk) :: temperature      !< system temperature
    real(rk) :: beta             !< inverse temperature

    ! numerics
    logical     :: is_production  !< track production
    integer(ik) :: nequilibration !< #equilibration cycles
    integer(ik) :: nproduction    !< #production cycles
    integer(ik) :: current_cycle  !< counter
    integer(ik) :: nmoves         !< #moves per cycle

    ! sampling
    character(len=50)  :: deffnm       !< defaultname for output data
    integer(ik)        :: stdout_freq  !< print every x'th cycle to stdout
  contains
    procedure, public :: init
    procedure, public :: calc_pressure_from_virial
  end type simulation_data_t

  type :: simulation_sampling_t
    type(sample_array_t)  :: h
    type(sample_array_t)  :: coord
    type(sample_scalar_t) :: pressure
    type(sample_scalar_t) :: energy
    type(sample_scalar_t) :: volume
  end type simulation_sampling_t
contains

  subroutine init (this, filename, deffnm, sample)
    class(simulation_data_t), intent(inout) :: this
    character(*)            , intent(in)    :: filename
    character(*)            , intent(in)    :: deffnm    
    type(simulation_sampling_t), intent(inout) :: sample
    real(rk) :: p_tail, u_tail

    this%deffnm = deffnm

    ! Read input from file
    call read_inputfile (this, sample, filename)


    ! Set the number of moves to number of particles
    this%nmoves = this%cell%npart

    ! Initialize pair potential
    call this%pot2b%init ([1.0_rk],[1.0_rk])

    ! Initialize cell
    ! 1_ik = lattice, 2_ik = fcc
    call this%cell%init (2_ik)
    call this%pot2b%eval_tail (this%cell%rc, this%cell%density, u_tail, p_tail)

    ! Transition matrix
    !call this%tm%init (161_ik, 5.8_rk, 9.0_rk, log(this%cell%volume))
    !call this%tm%init_from_file ('lambda.inp', 161_ik)
    !call this%wl%init_from_file ('lambda.inp', 161_ik, 0.001_rk, 0.1_rk)
    !call this%wltm%init_from_file ('lambda.inp', 161_ik, 0.00005_rk, 0.0005_rk, 0.15_rk)
    !call this%wltm%init_from_file (filename = 'lambda.inp', initial_id = 161_ik, &
    !  tol_stop_wl = 0.001_rk, tol_start_tm = 0.001_rk, p = 0.1_rk)

    print*, ''
    print*, 'Density: ', this%cell%density
    print*, 'temperature:', this%temperature
    print*, ''
  end subroutine init

  subroutine calc_pressure_from_virial (this, p_tail)
    class(simulation_data_t), intent(inout) :: this
    real(rk)                , intent(in)    :: p_tail 
    
    this%pressure_tensor = (unity * (real(this%cell%npart, rk) / this%beta) - this%virial) &
      / this%cell%volume + p_tail * unity

    this%pressure = trace (this%pressure_tensor) / 3._rk
  end subroutine calc_pressure_from_virial

  !> \brief reads input data from file(s)
  subroutine read_inputfile (sim, sample, filename)
    type(simulation_data_t) , intent(inout)    :: sim !< container type
    type(simulation_sampling_t), intent(inout) :: sample
    character(len=*)        , intent(in)    :: filename !< input file name
    character(len=100)                      :: buffer, label
    integer(ik) :: pos, ios = 0, line = 0 !< position in line
    integer(ik) :: ik_dummy
    integer(ik) :: fh !< filehandle = unit

    open (newunit = fh, file = filename, iostat = ios)
    if (ios /= 0) stop 'Error while opening inputfile!'
    ios = 0_ik

    read_file: do while (ios == 0_ik)
      read(fh, '(A)', iostat = ios) buffer
      if (ios == 0_ik) then
        line = line + 1_ik 

        ! find the first instance of whitespace (=position)
        pos = scan (buffer, ' ')
        ! check for comment (character at pos = 1: #)
        if (scan (buffer, '#') == 1) cycle read_file
        if (scan (buffer, ' ') == 1) cycle read_file
        ! split label and data
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)
        ! TODO: check for essential informations!
        ! Abort if essential parameters not set.


      case ('npart')
        read(buffer, *, iostat = ios) sim%cell%npart
      case ('nproduction')
        read(buffer, *, iostat = ios) sim%nproduction
      case ('nequilibration')
        read(buffer, *, iostat = ios) sim%nequilibration
     ! case ('nspecies')
     !   read(buffer, *, iostat = ios) sim%nspecies
      case ('ensemble')
        read(buffer, *, iostat = ios) sim%ensemble
      case ('stdout_freq')
        read(buffer, *, iostat = ios) sim%stdout_freq
      case ('press_freq')
        read(buffer, *, iostat = ios) ik_dummy
        call sample%pressure%init (trim(sim%deffnm)//'.pres.out', .true., ik_dummy)
      case ('h_freq')
        read(buffer, *, iostat = ios) ik_dummy
        call sample%h%init (trim(sim%deffnm)//'.h.out', .true., ik_dummy, [3_ik,3_ik])
      case ('ener_freq')
        read(buffer, *, iostat = ios) ik_dummy
        call sample%energy%init (trim(sim%deffnm)//'.ener.out', .true., ik_dummy)
      case ('coord_freq')
        read(buffer, *, iostat = ios) ik_dummy
        call sample%coord%init (trim(sim%deffnm)//'.xyz.out', .true., ik_dummy, [3_ik,sim%cell%npart])
      case ('vol_freq')
        read(buffer, *, iostat = ios) ik_dummy
        call sample%volume%init (trim(sim%deffnm)//'.vol.out', .true., ik_dummy)
   !   case ('init_cfg')
   !     read(buffer, *, iostat = ios) sim%initial_configuration

   !   case ('units')
   !     read(buffer, *, iostat = ios) sim%unitid

        ! Thermodynamic variables
      case ('temperature')
        read(buffer, *, iostat = ios) sim%temperature
        sim%beta = 1._rk/sim%temperature
      case ('density')
        read(buffer, *, iostat = ios) sim%cell%density
      case ('pressure')
        read(buffer, *, iostat = ios) sim%pressure_in
      case ('rc')
        read(buffer, *, iostat = ios) sim%cell%rc
      case ('rc_max')
        read(buffer, *, iostat = ios) sim%cell%rc_max
      case ('rc_const')
        call parse_logical (buffer, label, sim%cell%rc_is_constant)
      case ('isotropic')
        call parse_logical (buffer, label, sim%cell%is_isotropic)
        !read(buffer, *, iostat = ios) ik_dummy
        !if (ik_dummy == 0_ik) then
        !  sim%cell%is_isotropic = .false.
        !else if (ik_dummy == 1_ik) then
        !  sim%cell%is_isotropic = .true.
        !else
        !  stop 'Wrong input for ', label,'! Enter 0 or 1!'
        !end if
!      case ('acc_translation')
!        read(buffer, *, iostat = ios) sim%displacement_acc
!      case ('updaterate_translation')
!        read(buffer, *, iostat = ios) sim%displacement_freq

       ! Transition matrix method
      !case ('transition_matrix')
      !  read(buffer, *, iostat = ios) ik_dummy
      !  if (ik_dummy == 0_ik) then
      !    tm_active = .false.
      !  else 
      !    tm_active = .true.
      !  end if
      case default
        !                 print*, 'Skipping invalid entry in line: ', line
      end select
    end if

  end do read_file
end subroutine read_inputfile

subroutine parse_logical (buffer, label, logical_in)
  character(*), intent(in) :: buffer
  character(*), intent(in) :: label
  logical     , intent(inout) :: logical_in
  integer(ik) :: i
  
  ! read buffer
  read(buffer, *) i
  if (i == 0_ik) then
    logical_in = .false.
  else if (i == 1_ik) then
    logical_in = .true.
  else
    write(*,*) 'Wrong input for ', label,'! For logicals, enter either 0 or 1!'
    stop 
  end if
end subroutine parse_logical

end module simulation_data_mod
