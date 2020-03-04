!> \brief Implementation of the derived type of the simulation cell.
!> \author M.Theiss 

module cell_mod
  use working_prec_mod, only: ik, rk
  use util_mod        , only: det, inverse, &
    math_euclidean_norm, math_cross_product
  implicit none
  private

  public :: cell_t

  !> Simulation cell type.
  type :: cell_t
    logical                  :: is_isotropic = .false.
    logical                  :: rc_is_constant = .true.
    integer(ik)              :: npart       !< number of particles in this cell
    real(rk), allocatable    :: coord(:,:)  !< coordinates, dim(3,npart)
    real(rk) :: h(3,3)           !< box matrix
    real(rk) :: h_old(3,3)       !< old box matrix (to reset after declined volume move)
    real(rk) :: angles(3)        !< angles(alpha, beta, gamma)
    real(rk) :: da, db, dc = 0.0 !< lengths to calculate the inscribed sphere of the box (for rc)
    real(rk) :: rc, rc2, rc_min, rc_max, rc_old !< cut-off
    real(rk) :: volume !< volume of the cell
    real(rk) :: density !< volume of the cell
  contains
    private
    procedure, public :: init_from_file !< Initialization routine.
    procedure, public :: init_lattice  !< Initialization routine.
    procedure, public :: set_position     !< set a particle into the cell at position (incl PBC)
    procedure, public :: squared_distance !< Returns squared distance between two particles
    procedure, public :: connecting_vector !< Returns vector between i and j
    procedure, public :: change_random !< Change the cell randomly (is_isotropic -> scale, else deform)
    procedure, public :: restore !< Restore the box after a declined move
    procedure         :: set_rc !< Set the cut off according to the current box size/shape.
    generic, public   :: init => init_from_file, init_lattice !< Generic Initialization call
  end type cell_t
contains

  !> Initialization routine.
  subroutine init_from_file (this, config_file)
    class(cell_t), intent(inout) :: this !< dummy argument
    character(*) , intent(in)    :: config_file
  end subroutine init_from_file

  !> Initialization routine.
  subroutine init_lattice (this, lattice_index)
    class(cell_t), intent(inout) :: this !< dummy argument
    integer(ik), optional, intent(in)    :: lattice_index
    integer(ik) :: i
    real(rk)    :: r(3)

    ! allocate array
    allocate (this%coord(3,this%npart))

    !print*, 'DENSITY MANUALLY SET! INPUT IGNORED!!'
    !this%density = real(this%npart, rk) / exp(5.8_rk)
    !this%density = real(this%npart, rk) / exp(9._rk)

    ! set h matrix
    this%h = 0.0_rk
    do i = 1,3
      this%h(i,i) = ( real(this%npart, rk)/this%density )**(1._rk/3._rk)
    end do

    ! set cut-off and volume
    ! note: if the rc is constant, this wont change rc but throw an error, if rc is too large for box
    call this%set_rc ()

    ! set lattice or generate random
    if (present(lattice_index)) then
      if (lattice_index == 1_ik) call set_lattice (this)
      if (lattice_index == 2_ik) call set_fcc (this)
      if (lattice_index == 3_ik) call set_bcc (this)
    else
      do i = 1, this%npart
        call random_number(r)
        this%coord(:,i) = r
      end do
    end if
    

    write(*,'(X,A,X,F8.4)'), 'Density:', this%density
    write(*,'(X,A,X,F4.2)'), 'Cut-off radius:', this%rc
    write(*,'(X,A,X,9(F6.3,X))'), 'H:', this%h
  end subroutine init_lattice

  !> Returns the squared distance (incl nearest image) between particle i and j
  elemental real(rk) function squared_distance (this, i, j)
    class(cell_t), intent(in) :: this !< dummy argument
    integer(ik), intent(in) :: i, j             !< particle indices
    real(rk) :: dr(3)                           !< local, store the connection vector

    ! nearest image convention
    ! positions are INTERNAL coordinates [0,1], which makes scaling with
    ! boxlength obsolete. (see ortho implementation)
    dr = this%coord(:,i) - this%coord(:,j) - NINT( this%coord(:,i) - this%coord(:,j) ) 

    ! transform INTERNAL -> CARTESIAN distance
    dr = MATMUL (this%h, dr)

    squared_distance = SUM (dr*dr)
  end function squared_distance

  !> Set the particle i to position r applying PBC for triclinic cells.
  !! The position has to be in internal coordiantes.
  subroutine set_position (this, i, r)
    class(cell_t), intent(inout) :: this !< dummy
    integer(ik), intent(in) :: i                   !< index of particle 
    real(rk)   , intent(in) :: r(3)                !< position of particle
    integer(ik) :: j                               !< loop variable

    do j = 1,3
      if ( r(j) > 1._rk ) then
        this%coord(j,i) = r(j) - 1._rk
      else if ( r(j) < 0._rk ) then
        this%coord(j,i) = r(j) + 1._rk
      else
        this%coord(j,i) = r(j)
      end if
    end do

  end subroutine set_position

  !> Return distance vector between two particles (incl. nearest image)
  pure function connecting_vector (this, i, j) result (dr)
    class(cell_t), intent(in) :: this !< dummy argument
    integer(ik), intent(in)   :: i, j !< particle indices
    real(rk)                  :: dr(3)

    ! nearest image convention
    ! positions are INTERNAL coordinates [0,1], which makes scaling with
    ! boxlength obsolete. 
    dr = this%coord(:,i) - this%coord(:,j) - nint( this%coord(:,i) - this%coord(:,j) ) 

    ! transform INTERNAL -> CARTESIAN distance
    dr = matmul (this%h, dr)
  end function connecting_vector

  subroutine set_rc (this)
    class(cell_t), intent(inout) :: this !< dummy argument

    ! store old values
    this%rc_old = this%rc

    ! calculate the radius of the smallest inscribed sphere.
    this%volume = det (this%h)
    this%da = this%volume / math_euclidean_norm( math_cross_product (this%h(:,2), this%h(:,3)))
    this%db = this%volume / math_euclidean_norm( math_cross_product (this%h(:,1), this%h(:,3)))
    this%dc = this%volume / math_euclidean_norm( math_cross_product (this%h(:,1), this%h(:,2)))
    this%rc_min = 0.5_rk * min( this%da, this%db, this%dc )

    ! if rc is constant, then just test if it is smaller than rc_min
    ! if rc is variable, set it to the inscribed sphere radius as long as does not exceed rc_max
    if (this%rc_is_constant) then
      if (this%rc > this%rc_min) stop 'RC too large for cell! Decrease RC!'
    else
      this%rc = min (this%rc_min, this%rc_max)
    end if

    this%rc2 = this%rc**2

  end subroutine set_rc

  !> Change simulation cell randomly, according to given maximum displacement.
  subroutine change_random (this, dv)
    class(cell_t), intent(inout) :: this !< dummy argument
    real(rk)     , intent(in)    :: dv   !< maximum displacement

    integer(ik) :: i
    real(rk) :: displace_ran(3,3) ! 9 random numbers for cell entries
    real(rk) :: displace_ran_s ! 1 random number for isotropic cell changes
 
    ! store old box
    this%h_old = this%h

    ! perform change depending on deformation scheme.
    if (this%is_isotropic) then
      call random_number(displace_ran_s)
      displace_ran_s = (displace_ran_s - 0.5_rk) * 2._rk*dv
      do i=1,3
        ! scale
        this%h(i,i) = this%h(i,i) + displace_ran_s
      end do
    else
      call random_number(displace_ran)
      displace_ran = (displace_ran - 0.5_rk) * 2._rk*dv
      this%h = this%h + displace_ran
    end if

    ! set/check new cut-off
    ! this will set a new cut-off unless rc_is_constant = .true.
    call this%set_rc ()
    this%density = real(this%npart, rk)/this%volume
  end subroutine change_random

  !> Restore the state of the old cell -> cell, density, volume, cut off ...
  subroutine restore (this)
    class(cell_t), intent(inout) :: this !< dummy argument

    this%h = this%h_old
    call this%set_rc ()
    this%density = real(this%npart, rk)/this%volume
  end subroutine restore


  ! set configuration to simple cubic lattice - internal coordinates
  subroutine set_lattice (cell)
    type(cell_t), intent(inout) :: cell
    integer(ik) :: n, i, j, k, itel
    real(rk)    :: del, dx, dy, dz

    ! output
    print*, 'Set startconfig -> cubic lattice!'
    print*, 'Isotropic?', cell%is_isotropic

    ! set initial configuration to simple lattice    
    n = int(cell%npart**(1._rk/3._rk)-1e-6_rk, ik) + 1_ik

    if (n == 0) n = 1
    del = 1._rk /real(n, rk) ! spacing
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
          if (itel < cell%npart) then
            itel = itel + 1
            cell%coord(:,itel) = [dx, dy, dz]
          end if
        end do
      end do
    end do
  end subroutine set_lattice

  subroutine set_fcc (cell)
    type(cell_t), intent(inout) :: cell
    integer(ik) :: n, i, j, k, itel, countatom
    real(rk)    :: del, dx, dy, dz

    ! output
    print*, 'Set startconfig -> cubic fcc!'
    print*, 'Isotropic?', cell%is_isotropic

    n = int((cell%npart/4.)**(1._rk/3._rk)-1e-6, ik)+1

    if (n == 0) n = 1
    del = 1._rk/REAL(n, rk)
    itel = 0
    dx = -del
    do i = 1, n
      dx = dx + del
      dy = -del
      do j = 1, n
        dy = dy + del
        dz = -del
        do k = 1, n
          dz = dz + del
          if (itel < cell%npart) then
            itel = itel + 1
            cell%coord(:,itel) = [dx, dy, dz]
            if (itel < cell%npart .and. k<=n .and. j<=n) then
              itel = itel + 1
              cell%coord(:,itel) = [dx, dy+0.5*del, dz+0.5*del]
            end if
            if (itel < cell%npart .and. k<=n .and. i<=n) then
              itel = itel + 1
              cell%coord(:,itel) = [dx+0.5*del, dy, dz+0.5*del]
            end if
            if (itel < cell%npart .and. j<=n .and. i<=n) then
              itel = itel + 1
              cell%coord(:,itel) = [dx+0.5*del, dy+0.5*del, dz]
            end if
          end if
        end do
      end do
    end do

    ! print lattice to file 
   ! write (23,'(A,I9)') 'MODEL',1
   ! write (23,'(A,3f9.3,3f7.2,A)')'CRYST1', cell%box(1,1), cell%box(2,2), cell%box(3,3), &
   !   90._rk, 90._rk, 90._rk, ' P 1           1'
   ! countatom = 0_ik
   ! do i = 1,cell%npart
   !   countatom = countatom + 1_ik
   !   write (23,'(A,I7,A,I12,4x,3f8.3)') 'ATOM',countatom,'  O',  &
   !     countatom, matmul(cell%box,cell%coord(:,i))
   ! end do
  end subroutine set_fcc

  subroutine set_bcc (cell)
    type(cell_t), intent(inout) :: cell
    integer(ik) :: n, i, j, k, itel, countatom
    real(rk)    :: del, dx, dy, dz

    ! output
    print*, 'Set startconfig -> cubic bcc!'
    print*, 'Isotropic?', cell%is_isotropic
    n = int((cell%npart/2._rk)**(1./3.)-1e-6, ik) + 1_ik
    !n = INT(0.5 *((sqrt(4. * npart**2.+1.)+2. * npart)**(1./3.)-1./(sqrt(4. * npart**2.+1.)+2. * npart)**(1./3.)+1.)-1e-6) + 1

    if (n == 0) n = 1
    del = 1._rk/real(n, rk)
    itel = 0
    dx = -del
    do i = 1, n
      dx = dx + del
      dy = -del
      do j = 1, n
        dy = dy + del
        dz = -del
        do k = 1, n
          dz = dz + del
          if (itel < cell%npart) then
            itel = itel + 1
            cell%coord(:,itel) = [dx, dy, dz]
          end if
        end do
        dz = -del
        do k = 1, n
          dz = dz + del
          if (itel < cell%npart .and. i<=n .and.j<=n .and.k<=n) then
            itel = itel + 1
            cell%coord(:,itel) = [dx + 0.5 * del, dy + 0.5 * del, dz + 0.5 * del]
          end if
        end do
      end do
    end do

    !! print lattice to file 
    !write (23,'(A,I9)') 'MODEL',1
    !write (23,'(A,3f9.3,3f7.2,A)')'CRYST1', cell%box(1,1), cell%box(2,2), cell%box(3,3), &
    !  90._rk, 90._rk, 90._rk, ' P 1           1'
    !countatom = 0_ik
    !do i = 1,cell%npart
    !  countatom = countatom + 1_ik
    !  write (23,'(A,I7,A,I12,4x,3f8.3)') 'ATOM',countatom,'  O',  &
    !    countatom, matmul(cell%box,cell%coord(:,i))
    !end do
  end subroutine set_bcc


end module cell_mod
