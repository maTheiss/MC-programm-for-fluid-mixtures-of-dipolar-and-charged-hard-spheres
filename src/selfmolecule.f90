 module SelfMolecule_mod
  use working_prec_mod, only: rk, ik
  use globals_mod
  implicit none
  private

  public :: SelfMolecule   
  public :: SelfMolecule2
 
contains
 
 !> Self dipole-dipole term 
 ! first term of A+T in USCH
  subroutine SelfMolecule(  ENERGY )
    real(rk) , intent(out) :: ENERGY 
    real(rk):: d2
    integer(ik):: i, j 
    integer(ik) :: chargeid(2)
    
    ENERGY = 0.0_rk  
    
   do i = 1, npart
    if ( i .gt. nhs + nions) then 
    chargeid = return_charge_index (i)
!             ENERGY = ENERGY + ( Alpha/boxlength*charge_list(nions+i*j)**2._rk/sqrt(PI)  ) 
!         d2 = squared_distance (charge_positions(:,chargeid(1)), charge_positions(:,chargeid(2))) 
        ENERGY = ENERGY + charge_list(chargeid(1)) * charge_list(chargeid(2)) * erf ( Alpha/boxlength * dipole_sepa ) / dipole_sepa
    end if 
   end do 
   
    
  end subroutine SelfMolecule
 
 
subroutine SelfMolecule2( Nc, Xc, Yc, Zc, CHARGEc, ENERGY)
integer(ik), intent(in)              :: Nc
real(rk), dimension(Nc), intent(in)  :: Xc, Yc, Zc
real(rk), dimension(Nc), intent(in)  :: CHARGEc
real(rk), intent(out)                :: ENERGY 
! Local variables
integer:: i, j 
real(rk):: xij, yij, zij
real(rk):: xi, yi, zi
real(rk):: rij

 ENERGY = 0.0

do i = 1, Nc - 1

    do j = i + 1, Nc
            
            xij = abs( Xc(j) - Xc(i) )
            yij = abs( Yc(j) - Yc(i) )
            zij = abs( Zc(j) - Zc(i) )

            if( xij > boxlength - xij ) xij = xij - boxlength
            if( yij > boxlength - yij ) yij = yij - boxlength
            if( zij > boxlength - zij ) zij = zij - boxlength

            rij = sqrt(xij*xij + yij*yij + zij*zij)
            
            ENERGY = ENERGY + CHARGEc(i)  * CHARGEc(j) * &
                        erf( Alpha * rij / boxlength ) / rij

    end do
end do

return

end subroutine SelfMolecule2


end module SelfMolecule_mod




