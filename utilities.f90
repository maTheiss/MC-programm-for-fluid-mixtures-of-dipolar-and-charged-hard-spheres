module utilities_mod
  use working_prec_mod, only: rk, ik  
  use globals_mod
  implicit none
  private
  
  public :: outfold
 
contains
   
 
!--------------------------------------------------------------------------------------------------
! outfold
! This subroutine reconstructs a molecule which might be divided into several
! parts because of the periodic boundary conditions.
!--------------------------------------------------------------------------------------------------
subroutine outfold ( pid, stion, endion ) 
integer(ik), intent(in) :: pid
integer(ik), intent(in) :: stion  
integer(ik), intent(in) :: endion   
integer(ik)             :: i
real(rk), dimension(3)  :: xref, dx

    xref(1:3) = particle_positions(1:3,pid)  !reference

     do i = stion, endion 

        dx(1:3) = charge_positions(1:3,i) - xref(1:3)
        
        charge_positions(1:3,i) = charge_positions(1:3,i) - nint(dx(1:3)/boxlength) * boxlength 
        
        xref(1:3) = charge_positions(1:3,i)
          
     end do 

end subroutine outfold 

 

end module utilities_mod