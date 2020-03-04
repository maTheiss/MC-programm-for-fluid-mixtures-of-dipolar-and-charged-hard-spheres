 module RealMolecule_mod
  use working_prec_mod, only: rk, ik
  use globals_mod
  implicit none
  private

  public :: RealMolecule  
  public :: RealMolecule2
  public :: RealMolecule_move  
 
contains
 
 subroutine RealMolecule(  ENERGY )
    real(rk) , intent(out) :: ENERGY 
    real(rk):: d2
    real(rk):: CHARGEi 
    real(rk):: cp1(3), cp2(3)
    real(rk):: xij, yij, zij, rij
    real(rk):: xi, yi, zi
    integer(ik):: i, iii, j , chargeid(2) 
    integer(ik):: stion, endion 
    
    ENERGY = 0.0_rk  
  
  do iii = 1, npart 
   
    chargeid = return_charge_index(iii)
    stion = chargeid(1)
    endion = stion + LENGTH(iii)-1_ik
    
    do i = stion, endion
    
        xi = charge_positions(1,i)  
        yi = charge_positions(2,i)
        zi = charge_positions(3,i)
        
        CHARGEi = charge_list(i)
        
        do j = 1, stion-1_ik 
        
            xij = abs( charge_positions(1,j) - xi )
            yij = abs( charge_positions(2,j) - yi )
            zij = abs( charge_positions(3,j) - zi )

            if( xij > boxlength - xij ) xij = xij - boxlength
            if( yij > boxlength - yij ) yij = yij - boxlength
            if( zij > boxlength - zij ) zij = zij - boxlength

            rij = sqrt(xij*xij + yij*yij + zij*zij)
 
                ENERGY = ENERGY + CHARGEi * charge_list(j)  * &
                    erfc( Alpha	* rij / boxlength ) / rij
                     

        end do
        do j = endion+1, ncharges
 
            xij = abs( charge_positions(1,j) - xi )
            yij = abs( charge_positions(2,j) - yi )
            zij = abs( charge_positions(3,j) - zi )

            if( xij > boxlength - xij ) xij = xij - boxlength
            if( yij > boxlength - yij ) yij = yij - boxlength
            if( zij > boxlength - zij ) zij = zij - boxlength

            rij = sqrt(xij*xij + yij*yij + zij*zij)
 
                ENERGY = ENERGY + CHARGEi * charge_list(j) * &
                    erfc( Alpha	* rij / boxlength ) / rij
                     

        end do
    
    
    end do 
    
 end do    
ENERGY = 0.5_rk*ENERGY
     
    
 end subroutine RealMolecule

 

 subroutine RealMolecule_move(  stion, endion, ENERGY )        
    integer(ik), intent(in) :: stion, endion
    real(rk) , intent(out)  :: ENERGY 
    real(rk):: d2, sqrtd2, CHARGEi 
    real(rk):: xij, yij, zij, rij
    real(rk) :: xi, yi, zi
    integer(ik):: i, j, Start, Finish!, LENGTH 
    integer(ik) :: chargeid(2)
        
    ENERGY = 0.0_rk  
    
    do i = stion, endion 
        
        xi = charge_positions(1,i)  
        yi = charge_positions(2,i)
        zi = charge_positions(3,i)
        
        CHARGEi = charge_list(i)
        
        do j = 1, stion-1_ik
            
            xij = abs( charge_positions(1,j) - xi )
            yij = abs( charge_positions(2,j) - yi )
            zij = abs( charge_positions(3,j) - zi )

            if( xij > boxlength - xij ) xij = xij - boxlength
            if( yij > boxlength - yij ) yij = yij - boxlength
            if( zij > boxlength - zij ) zij = zij - boxlength 
        
            rij = sqrt(xij*xij + yij*yij + zij*zij)
 
                ENERGY = ENERGY + CHARGEi * charge_list(j)  * &
                    erfc( Alpha * rij / boxlength ) / rij
        
        end do 
        do j = endion+1, ncharges
            
            
            xij = abs( charge_positions(1,j) - xi )
            yij = abs( charge_positions(2,j) - yi )
            zij = abs( charge_positions(3,j) - zi )

            if( xij > boxlength - xij ) xij = xij - boxlength
            if( yij > boxlength - yij ) yij = yij - boxlength
            if( zij > boxlength - zij ) zij = zij - boxlength
            
            rij = sqrt(xij*xij + yij*yij + zij*zij)
 
                ENERGY = ENERGY + CHARGEi * charge_list(j) * &
                    erfc( Alpha * rij / boxlength ) / rij
              
        end do
    
    end do 
      
     
 
 end subroutine RealMolecule_move

 
 
 
subroutine RealMolecule2( Nc, Xc, Yc, Zc, CHARGEc, Np, Xp, Yp, Zp, CHARGEp, &
                                                  ENERGY )
 

! This routine calculates the real Ewald sum energy of Nc ionic beads 
! interacting with Np ionic beads.

! Nc is the number of ionic beads in one group.
! TYPEc contains the group identity of the Nc beads.
! Xc, Yc, Zc are the coordinates of the Nc beads.

integer, intent(in)                 :: Nc
real(rk), dimension(Nc), intent(in)  :: CHARGEc
real(rk), dimension(Nc), intent(in)     :: Xc, Yc, Zc

! Np is the number of ionic beads in another group.
! TYPEp contains the group identity of the Np beads.
! Xp, Yp, Zp are the coordinates of the Np beads.
                                                                    
integer, intent(in)     :: Np
real(rk), dimension(Np), intent(in):: Xp, Yp, Zp
real(rk), dimension(Np), intent(in):: CHARGEp
 
! ! real, dimension(Niongrs), intent(in)	:: CHARGE

! BoxSize is the length of the simulation box.

! real, intent(in):: BoxSize

! Alpha is an Ewald sum parameter, Alpha = kappa * L, for kappa in A + T.

! real, intent(in):: Alpha

! erfc is the complementary error function.

! ! real, external:: erfc

! ENERGY contains the energy of interaction between group c
! and group p for each hamiltonian.

real(rk), intent(out):: ENERGY

! Local variables

integer:: i, j 
real(rk):: xij, yij, zij
real(rk):: xi, yi, zi
real(rk):: rij


ENERGY = 0.0

do i = 1, Nc

    xi = Xc(i)
    yi = Yc(i)
    zi = Zc(i)
    

    do j = 1, Np
 
            
            xij = abs( Xp(j) - xi )
            yij = abs( Yp(j) - yi )
            zij = abs( Zp(j) - zi )

            if( xij > boxlength - xij ) xij = xij - boxlength
            if( yij > boxlength - yij ) yij = yij - boxlength
            if( zij > boxlength - zij ) zij = zij - boxlength

            rij = sqrt(xij*xij + yij*yij + zij*zij)

            if (rij .NE. 0.0 .AND. CHARGEc(i) * CHARGEp(j)  .NE. 0.0) then
            
                ENERGY = ENERGY + CHARGEc(i) * CHARGEp(j)  * &
                    erfc( Alpha	* rij / boxlength ) / rij
                    
            end if


    end do

end do

return

end subroutine RealMolecule2
 

 
end module RealMolecule_mod




