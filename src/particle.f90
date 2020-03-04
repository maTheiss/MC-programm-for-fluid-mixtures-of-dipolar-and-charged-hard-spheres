!> module containing the particle types
module particle_mod
  use working_prec_mod, only: rk, ik
  implicit none
  private
  
  public:: hs_t
  public:: ion_t
  public:: dipole_t

  !> general particle type
  type hs_t
    integer(ik) :: spid !< id of species
    real(rk)    :: pos(3) !< position of particle
  end type hs_t

  !> ion extension of hs
  type, extends(hs_t) :: ion_t
    real(rk) :: q !< charge strength
  !  real(rk) :: posi(3)
  end type ion_t

  !> dipole extension of hs
  type, extends(hs_t) :: dipole_t
    real(rk) :: q1, q2 !< partial charges
    real(rk) :: pos1(3), pos2(3) !< positions of charges
  end type dipole_t

contains


!  subroutine init_hs (hs, sig)
!    type(hs_t), intent(inout) :: hs

!  end subroutine init_hs

end module particle_mod
