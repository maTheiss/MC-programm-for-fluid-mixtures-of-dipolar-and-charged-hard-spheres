
subroutine RealInteract( ENERGY )
use working_prec_mod , only: rk, ik 
use globals_mod
use realmolecule_mod
implicit none

! This routine calculates the total coulombic energy of a group of molecules.

! Nb is the total number of ionic beads in a system.
! TYPEb contains the group identity of each bead.
! Xb, Yb, Zb are the coordinates of each bead.

! integer, intent(in)									:: Nb
! integer, dimension(Nb), intent(in)					:: TYPEb
! real, dimension(Nb), intent(in)						:: Xb, Yb, Zb
! 
! ! Nmol is the number of molecules in the system.
! ! LENGTH contains the number of ionic beads in each molecule.
! 
! integer, intent(in)									:: Nmol
! integer, dimension(Nmol), intent(in)				:: LENGTH

real(rk), intent(out)   :: ENERGY

! Local variables

integer(ik) :: i, Nc, Nb
integer(ik) :: Start, Finish
real(rk)    :: ENERGY_part
integer(ik) :: chargeid(2)
! integer(ik) :: LENGTH 


ENERGY = 0.0_rk

Finish = 0_ik

if (ncharges == 0) return 

Nb = ncharges  

do i = 1, npart - 1
 
    Nc = LENGTH(i)

    if( Nc == 0 ) cycle

    Start = Finish + 1_ik        
    Finish = Start + LENGTH(i) - 1_ik
    
    if( Finish == Nb ) cycle
    
    call Realmolecule2( Nc, charge_positions(1,Start:Finish), charge_positions(2,Start:Finish), &
                                        charge_positions(3,Start:Finish),charge_list(Start:Finish), &
                                        Nb - Finish, charge_positions(1,Finish+1:Nb), & 
                                        charge_positions(2,Finish+1:Nb), &
                                        charge_positions(3,Finish+1:Nb), charge_list(Finish+1:Nb), &
                                        ENERGY_part)

    ENERGY = ENERGY + ENERGY_part 
    
end do 

return

end subroutine RealInteract


subroutine RealInteract_move( pid, ENERGY )
use working_prec_mod , only: rk, ik 
use globals_mod
use realmolecule_mod
implicit none

! This routine calculates the total coulombic energy of a group of molecules.

! Nb is the total number of ionic beads in a system.
! TYPEb contains the group identity of each bead.
! Xb, Yb, Zb are the coordinates of each bead.

! integer, intent(in)									:: Nb
! integer, dimension(Nb), intent(in)					:: TYPEb
! real, dimension(Nb), intent(in)						:: Xb, Yb, Zb
! 
! ! Nmol is the number of molecules in the system.
! ! LENGTH contains the number of ionic beads in each molecule.
! 
! integer, intent(in)									:: Nmol
! integer, dimension(Nmol), intent(in)				:: LENGTH
integer(ik), intent(in):: pid
real(rk), intent(out)   :: ENERGY

! Local variables

integer :: i, Nc, Nb
integer :: Start, Finish
real(rk)    :: ENERGY_part
integer(ik) :: chargeid(2)
! integer(ik) :: LENGTH 
real(rk),dimension(ncharges):: temp5, temp6, temp7, temp8
   
 
ENERGY = 0.0

Finish = 0

if (ncharges == 0) return 

Nb = ncharges  

! do i = 1, npart - 1

    chargeid = return_charge_index (pid)
!     LENGTH = 1_ik
!     if ( chargeid(2) /= 0 ) LENGTH = 2_ik
    
    Nc = LENGTH(pid)
 
    Start = chargeid(1)  
    Finish = Start + Nc - 1_ik
     
    if( Start == 1 ) then
    
                call Realmolecule2( Nc, charge_positions(1,Start:Finish), charge_positions(2,Start:Finish), &
                                        charge_positions(3,Start:Finish),charge_list(Start:Finish), &  !Nb - Nc
                                        Nb-Nc, charge_positions(1,Finish+1:Nb), & 
                                        charge_positions(2,Finish+1:Nb), &
                                        charge_positions(3,Finish+1:Nb), charge_list(Finish+1:Nb), &
                                        ENERGY_part)
       
    else if( Start + Nc - 1 == ncharges ) then
    
                call Realmolecule2( Nc, charge_positions(1,Start:Finish), charge_positions(2,Start:Finish), &
                                        charge_positions(3,Start:Finish),charge_list(Start:Finish), &
                                        Nb-Nc, charge_positions(1,1:Start-1), &  
                                        charge_positions(2,1:Start-1), &
                                        charge_positions(3,1:Start-1), charge_list(1:Start-1), &
                                        ENERGY_part)
    else 
    
    
                temp5( 1:Start-1 ) = charge_positions(1, 1:Start-1 )
                temp5( Start:ncharges-Nc ) = charge_positions(1, Finish+1:ncharges )
                
                temp6( 1:Start-1 ) = charge_positions(2, 1:Start-1 )
                temp6( Start:ncharges-Nc ) = charge_positions(2, Finish+1:ncharges )
                
                temp7( 1:Start-1 ) = charge_positions(3, 1:Start-1 )
                temp7( Start:ncharges-Nc ) = charge_positions(3, Finish+1:ncharges )
                
                temp8( 1:Start-1 ) = charge_list(1:Start-1) 
                temp8( Start:ncharges-Nc ) = charge_list(Finish+1:ncharges)
                
                call Realmolecule2( Nc, charge_positions(1,Start:Finish), charge_positions(2,Start:Finish), &
                                        charge_positions(3,Start:Finish),charge_list(Start:Finish), &
                                        Nb-Nc, temp5, temp6, temp7, temp8, &
                                        ENERGY_part)
    
    end if 
    ENERGY = ENERGY + ENERGY_part

! end do
  
return

end subroutine RealInteract_move









