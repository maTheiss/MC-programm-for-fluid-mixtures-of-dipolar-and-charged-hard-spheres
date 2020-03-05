!> \brief global variables and important subroutines  
!> \author M.Theiss 
module globals_mod
  use working_prec_mod , only: rk, ik 
  use basic_parameters_mod 
  use sim_parameters_mod  
  use read_input_mod 
  implicit none

  real(rk)   , allocatable :: particle_positions (:,:) !< list containing hs positions (dim = 3,npar)
  real(rk)   , allocatable :: charge_positions (:,:) !< containing charge positions (dim = 3,ncharges)
  real(rk)   , allocatable :: sigma_list (:) !< containing hs diameters
  real(rk)   , allocatable :: charge_list (:) !< containing charge strength
  integer(ik), allocatable :: species_list (:) !< containing species index (dim = npart)
  integer(ik), allocatable :: LENGTH (:) !< containing species index (dim = npart)

  integer(ik) :: npart !< total number of particles
  integer(ik) :: nhs !< number of hard sphere particles
  integer(ik) :: nions !< number of ions
  integer(ik) :: ndipoles !< number of dipoles
  integer(ik) :: ncharges !< number of charges (= nions + ndipoles)
    
  real(rk), parameter       :: dipole_sepa = 0.1  ! separation of charges within a dipole = 0.1*sigma 
   
  real(rk) :: boxlength !< length of the cubic simulation box
  real(rk) :: rc !< cut off radius
  real(rk) :: rc2 !< squared cut off radius
  real(rk) :: density !< read in density  [-]
  real(rk) :: volume !< simulation volume; V = boxlength**3
  real(rk) :: Temp  !< Temperature [K]
  real(rk) :: beta  !< contains inverse Temp. : 1//
  real(rk) :: frac   !< ion mole fraction 
  integer(ik):: nsamp !< sample every x cycles
  real(rk)   :: q_ions, mu_dipoles !< dimensionless charge of ions and dimensionless dipole moment
  real(rk)  :: lambda  !< perturbation parameter 
  
contains

!---------------------------------------------------------------------------------------------------------------------
!> initialize system 
!---------------------------------------------------------------------------------------------------------------------
  subroutine init_globals (read_conf, first_run, nhs_in, nions_in, ndipoles_in &
                            , density_in, Temp_in, sigma_in)
    logical,     intent(in) :: read_conf
    logical,     intent(in) :: first_run 
    integer(ik), intent(in) :: nhs_in
    integer(ik), intent(in) :: nions_in
    integer(ik), intent(in) :: ndipoles_in
    real(rk)   , intent(in) :: density_in
    real(rk)   , intent(in) :: Temp_in 
    real(rk)   , intent(in) :: sigma_in 
    integer(ik) :: i, nions_test 
    real(rk)    :: nions_test1
     
    CoulCombo = 167100.002729011_rk  !ec * ec * 1.0e10_rk / ( 4.0_rk * PI * eps0 * kB )
       
    nhs = nhs_in
    nions = nions_in
    ndipoles = ndipoles_in
    ncharges = nions + 2 * ndipoles
    npart = nions + ndipoles + nhs 
     
    density = density_in
    Temp = Temp_in 
    beta = 1.0_rk / Temp 
    
    ! set volume and box
    volume = real(npart, rk) / density   !Annahme: sig=1
    boxlength = volume ** (1._rk/3._rk)  !..dim.los (dann boxlength* = boxlength/sig) oder dim.behaftet mit sig=1 angst.  
    ! set cut off to half of the boxlength
    rc = 0.5_rk * boxlength
    rc2 = rc * rc
    
    write(*,"(A,F10.5)") 'boxlength:', boxlength 
    write(*,"(A,F10.5)") 'density  :', density  
    write(*,"(A,F10.5)") 'fraction :',real(nions,rk)/real(npart,rk) 
    
  !-----------------------------------------------------------------------------
  ! allocate basic quantities (for vectors of length Nsp)
  !-----------------------------------------------------------------------------
    allocate (particle_positions(3,npart))
    allocate (charge_positions(3,ncharges))
    allocate (sigma_list(npart))
    allocate (charge_list(ncharges))
    allocate (species_list(npart)) 
    allocate (LENGTH(npart))
    allocate (EXPX(0:Kmax, ncharges))
    allocate (EXPY(-Kmax:Kmax, ncharges))
    allocate (EXPZ(-Kmax:Kmax, ncharges))
    
    
    ! set forcefield parameters
    do i = 1, npart 
      sigma_list(i) = sigma_in
    end do 
    
    !def. species_list & LENGTH:
    !species_list: hs <-> 1, ion <-> 2, dipole <-> 3
    !LENGTH:       hs <-> 0, ion <-> 1, dipole <-> 2
    if ( nhs /= 0_ik )  then
      do i = 1, nhs
        species_list(i) = 1_ik 
        LENGTH(i) = 0_ik 
      end do 
    end if
     if ( nions /= 0_ik )  then
      do i = nhs+1, nhs + nions
        species_list(i) = 2_ik 
        LENGTH(i) = 1_ik 
      end do 
    end if
    if ( ndipoles /= 0_ik )  then 
      do i = nhs+nions+1, npart 
        species_list(i) = 3_ik 
        LENGTH(i) = 2_ik 
      end do 
    end if  
    
    
    !distributing dimensionless charges 
    !..to ions
    if ( nions /= 0 ) then 
     do i = 1, nions/2
       charge_list(2*i) = sqrt(q_ions)
     end do 
     do i = 0, nions/2-1
       charge_list(2*i+1) = -sqrt(q_ions)  
     end do 
   end if 
   !..to dipoles
   if ( ndipoles /= 0 ) then 
    do i = 1, (ncharges-nions)/2 
       charge_list(nions+2*i) = sqrt(mu_dipoles*sigma_list(npart)**2.0_rk/dipole_sepa**2.0_rk)
    end do 
    do i = 0, (ncharges-nions)/2-1   
       charge_list(nions+2*i+1) = -sqrt(mu_dipoles*sigma_list(npart)**2.0_rk/dipole_sepa**2.0_rk) 
    end do 
   end if    
     
    !initialize counter
    total_natt = 0_ik
    total_nacc = 0_ik
    natt = 0_ik
    nacc = 0_ik
    
    call identifier_files
    
     if ( .not. read_conf ) then
        !set particles on a simple cubic lattice
        !call set_lattice ()   ! wohl nicht für höhere Dichten geeignet
        call lattice_fcc   
        !init. max displacement
        max_disp(:)=0.1
        max_rot(:) =0.2
     else if ( .not. first_run ) then  
        call restarted_files 
        call read_positions
        call read_res_file(nhs, nions, ndipoles, npart, lambda, boxlength, density)  !reading max_disp/rot and check data with input file
        call check_overlap
     else   
        call read_positions
        call read_res_file(nhs, nions, ndipoles, npart, lambda, boxlength, density)  !reading max_disp/rot and check data with input file
        !additional "_res" in files 
        call restarted_files 
        call check_overlap
     end if      
    print*, ''
    print*, 'System initialized!.'
  end subroutine init_globals
  
!---------------------------------------------------------------------------------------------------------------------
!> returns charge index(number) of given particle index(number)
!---------------------------------------------------------------------------------------------------------------------
  function return_charge_index (pid_in) result (res)
    integer(ik), intent(in) :: pid_in
    integer(ik) :: res(2)

    res = 0_ik

    ! is it hs, ion or dipole?
    if ( pid_in .gt. nhs .and. pid_in .le. nhs+nions) then                  ! if (pid_in <= nhs + nions) then...   <-> res(1) sonst auch neg. für hs-teilchen
      res(1) = pid_in - nhs
    else 
      res(1) = (pid_in- (nhs + nions) -1_ik)*2_ik + nions + 1_ik   !pid_in - (nhs + nions)
      res(2) = res(1) + 1_ik
    end if
  end function return_charge_index

!---------------------------------------------------------------------------------------------------------------------
! returns squared distance of two particles
!---------------------------------------------------------------------------------------------------------------------
  !> return squared distance including nearest image convention
  real(rk) function squared_distance (pos1, pos2)
    real(rk), intent(in) :: pos1(3) !< position of particle 1
    real(rk), intent(in) :: pos2(3) !< position of particle 2
    real(rk) :: dr(3) !< dummy variable

    dr = pos2 - pos1
    ! nearest image convention
!     dr = dr - nint(dr / boxlength) * boxlength
    squared_distance = sum(dr * dr)  !vektorlänge² = x²+y²+z²

  end function squared_distance
  
     

!---------------------------------------------------------------------------------------------------------------------
!> check if particle with pid_in overlaps with other particles
!  if yes -> return true
!---------------------------------------------------------------------------------------------------------------------
  logical function is_overlapping (pid_in)
    integer(ik), intent(in) :: pid_in
    integer(ik) :: i
    real(rk)    :: d2 
    real(rk)    :: sig1, sig2, sig3
    real(rk)    :: xi, yi, zi 
    real(rk)    :: xij, yij, zij 

    is_overlapping = .false. 
    
    xi = particle_positions(1,pid_in)
    yi = particle_positions(2,pid_in)
    zi = particle_positions(3,pid_in)
    
    sig1 = sigma_list(pid_in)
    
    do i = 1, pid_in-1
         
        xij = abs( particle_positions(1,i) - xi )
        yij = abs( particle_positions(2,i) - yi )
        zij = abs( particle_positions(3,i) - zi )

        if( xij > boxlength - xij ) xij = xij - boxlength
        if( yij > boxlength - yij ) yij = yij - boxlength
        if( zij > boxlength - zij ) zij = zij - boxlength

        d2 = xij*xij + yij*yij + zij*zij
             
        sig2 = sigma_list(i)
        sig3 = (0.5_rk * (sig1 + sig2)) ** 2.0_rk
        if (d2 <= sig3) then
            is_overlapping = .true.
            return
        end if
    end do

    do i = pid_in+1, npart
    
        xij = abs( particle_positions(1,i) - xi )
        yij = abs( particle_positions(2,i) - yi )
        zij = abs( particle_positions(3,i) - zi )

        if( xij > boxlength - xij ) xij = xij - boxlength
        if( yij > boxlength - yij ) yij = yij - boxlength
        if( zij > boxlength - zij ) zij = zij - boxlength
        
        d2 = xij*xij + yij*yij + zij*zij
       
        sig2 = sigma_list(i)
        sig3 = (0.5_rk * (sig1 + sig2)) ** 2.0_rk
        if (d2 <= sig3) then
            is_overlapping = .true.
            return
        end if
    end do
    
     
    
  end function is_overlapping

!---------------------------------------------------------------------------------------------------------------------
!> check if particles from input files are overlapping 
!  if yes -> return true
!---------------------------------------------------------------------------------------------------------------------
subroutine check_overlap  
integer(ik) :: i 
  
do i = 1, npart

    if (is_overlapping(i)) then 
    
      print*, 'particles in input file are overlapping'  
      stop
    
    end if    
    
end do 

end subroutine check_overlap


!---------------------------------------------------------------------------------------------------------------------
!> set configuration to simple cubic lattice
!---------------------------------------------------------------------------------------------------------------------
  subroutine set_lattice ()
    integer(ik) :: n, i, j, k, itel, ii
    real(rk)    :: del, dx, dy, dz
    real(rk):: ran 
    integer(ik) :: chargeid(2) !< id(s) of ion (dipole) charge(s)
    
    ! output
    print*, 'Set startconfig -> cubic lattice!'

    ! set initial configuration to simple lattice    
    n = int(npart**(1._rk/3._rk)-1e-6_rk, ik) + 1_ik

    if (n == 0) n = 1
    del = boxlength /real(n, rk) ! spacing in Angst.
    print*, boxlength, del,n
    itel = 0_ik
    dx = -del
    do i = 1, n
      dx = dx + del
      dy = -del
      do j = 1, n
        dy = dy + del
        dz = -del
        do k = 1, n
          dz = dz + del
          if (itel < npart) then
            itel = itel + 1
            particle_positions(:,itel) = [dx, dy, dz]   !Angst.
          end if
        end do
      end do
    end do
    
   if (nions /= 0_ik) then 
    do i = 1, nions
        charge_positions(:,i) = particle_positions(:,i+nhs)
    end do 
   end if 
     

   if (ndipoles /= 0_ik) then 
    do ii = 1, npart 
      if ( ii .gt. nions + nhs ) then 
        chargeid = return_charge_index (ii)
        call random_number(ran)

        charge_positions(1,chargeid(1)) =  particle_positions(1,ii) + &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * sin(ran*PI) * cos(ran*2.0_rk*PI-PI)
        charge_positions(2,chargeid(1)) =  particle_positions(2,ii) + &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * sin(ran*PI) * sin(ran*2.0_rk*PI-PI)
        charge_positions(3,chargeid(1)) =  particle_positions(3,ii) + &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * cos(ran*PI)  
        charge_positions(1,chargeid(2)) =  particle_positions(1,ii) - &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * sin(ran*PI) * cos(ran*2.0_rk*PI-PI)
        charge_positions(2,chargeid(2)) =  particle_positions(2,ii) - &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * sin(ran*PI) * sin(ran*2.0_rk*PI-PI)
        charge_positions(3,chargeid(2)) =  particle_positions(3,ii) - &
                                                0.5_rk*dipole_sepa*sigma_list(ii) * cos(ran*PI)  
        
     end if 
    end do 
   end if  
   
  end subroutine set_lattice

end module globals_mod
 








