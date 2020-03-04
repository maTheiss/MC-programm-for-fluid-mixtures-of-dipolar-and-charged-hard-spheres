module particle_swap_mult_mod 
  use working_prec_mod, only: rk, ik
  use basic_parameters_mod
  use sim_parameters_mod
  use globals_mod
  use RealMolecule_mod
  use ES_Fourier_mod
  use utilities_mod
  use Energy_mod
  implicit none
  private 

  public:: particle_swap_mult

contains
subroutine particle_swap_mult( dUREAL_sum, dUFOURIER_sum, ener_delt)

real(rk), intent(out) :: dUREAL_sum  
real(rk), intent(out) :: dUFOURIER_sum
real(rk), intent(out) :: ener_delt

integer(ik) :: MaxSwap, NSwap
integer(ik) :: i, ii, k
integer(ik) :: Mol1, Mol2
integer(ik) :: lenion1,  lenion2
integer(ik) :: chargeid(2) 
integer(ik) :: stion1, endion1, stion2, endion2

real(rk) :: x12, y12, z12  
real(rk) :: swap, zzMol1, zzMol2, acceptance 
real(rk) :: UREAL1_old, UREAL2_old, UREAL1_new, UREAL2_new
real(rk) :: U_corr1  
real(rk) :: dUFOURIER1, dUFOURIER2
  
real(rk), allocatable, dimension(:)         :: dUREAL, dUFOURIER, U_corr
real(rk), allocatable, dimension( : , : )   :: position_old  
real(rk), allocatable, dimension( : , : )   :: charge_position_old  

complex(rk), allocatable, dimension( : , : ):: EXPX1_NEW
complex(rk), allocatable, dimension( : , : ):: EXPY1_NEW, EXPZ1_NEW
complex(rk), allocatable, dimension( : , : ):: EXPX2_NEW
complex(rk), allocatable, dimension( : , : ):: EXPY2_NEW, EXPZ2_NEW

complex(rk), dimension(Nkvec)               :: SUMQEXPV_NEW, SUMQEXPV_TMP 
complex(rk), dimension(0:Kmax, ncharges)    :: EXPX_OLD
complex(rk), dimension(-Kmax:Kmax, ncharges):: EXPY_OLD, EXPZ_OLD
complex(rk), dimension(Nkvec)               :: SUMQEXPV_OLD
 
 
 dUREAL_sum  = 0.0_rk
 dUFOURIER_sum = 0.0_rk
 
 total_attempted_swap = total_attempted_swap + 1_ik
 
 MaxSwap = 5

 call random_number(swap)
 NSwap = int(swap*MaxSwap)+1

 allocate(dUREAL(NSwap))
 allocate(U_corr(NSwap))
 allocate(dUFOURIER(NSwap))
  
 dUREAL = 0.
 dUFOURIER = 0.
 U_corr = 0.0_rk
 U_corr1 = 0. 

 EXPX_OLD = EXPX
 EXPY_OLD = EXPY
 EXPZ_OLD = EXPZ

 SUMQEXPV_NEW = SUMQEXPV
 SUMQEXPV_OLD = SUMQEXPV

     
 allocate(charge_position_old(3,ncharges)) 
 allocate (position_old(3,npart)) 

 charge_position_old = charge_positions
 position_old = particle_positions 
 
kloop: do k = 1, NSwap

    call random_number (zzMol1 )
    call random_number (zzMol2 ) 

    Mol1 = int(zzMol1*npart)+1
    Mol2 = int(zzMol2*npart)+1

    if (Mol1 == Mol2) then
!     if (SPECIES(Mol1) .EQ. SPECIES(Mol2) .AND. SpRotate(SPECIES(Mol2)) == .false.) then

        dUREAL(k) = 0.0_rk
        dUFOURIER(k) = 0.0_rk 
        U_corr(k) = 0.0_rk
        cycle kloop

    end if
    
    x12 = particle_positions(1,Mol2) - particle_positions(1,Mol1)     !Xlj_new(Mol2) - Xlj_new(Mol1)
    y12 = particle_positions(2,Mol2) - particle_positions(2,Mol1)     !Ylj_new(Mol2) - Ylj_new(Mol1)
    z12 = particle_positions(3,Mol2) - particle_positions(3,Mol1)     !Zlj_new(Mol2) - Zlj_new(Mol1)
 
    lenion1 = LENGTH(Mol1)  
    lenion2 = LENGTH(Mol2)  
    
    chargeid = return_charge_index (Mol1) 
    stion1 = chargeid(1) 
    endion1 = stion1 + lenion1 - 1_ik
    
    chargeid = return_charge_index (Mol2) 
    stion2 = chargeid(1) 
    endion2 = stion2 + lenion2 - 1_ik
    
    allocate( EXPX1_NEW(0:Kmax, lenion1 ) )
    allocate( EXPY1_NEW(-Kmax:Kmax, lenion1 ) )
    allocate( EXPZ1_NEW(-Kmax:Kmax, lenion1 ) )

    allocate( EXPX2_NEW(0:Kmax, lenion2 ) )
    allocate( EXPY2_NEW(-Kmax:Kmax, lenion2 ) )
    allocate( EXPZ2_NEW(-Kmax:Kmax, lenion2 ) ) 
    
    if (lenion1 /= lenion2) then 
    
      call RealMolecule2( lenion1, charge_positions(1,stion1:endion1), &
                 charge_positions(2,stion1:endion1),charge_positions(3,stion1:endion1), &
                 charge_list(stion1:endion1), lenion2, charge_positions(1,stion2:endion2), &
                 charge_positions(2,stion2:endion2), charge_positions(3,stion2:endion2), &
                 charge_list(stion2:endion2), U_corr1 )
    else 
      U_corr1 = 0.
    end if 
    
    if (lenion1 .GT. 0) then
        call RealMolecule_move ( stion1, endion1, UREAL1_old)                     
    else
        UREAL1_old = 0.
    end if


    if (lenion2 .GT. 0) then
        call RealMolecule_move ( stion2, endion2, UREAL2_old)  
    else
        UREAL2_old = 0.
    end if
   
    ! reconstructing position of dipole charges 
    if (species_list(Mol1)==3) call outfold ( Mol1, stion1, endion1 ) 
    if (species_list(Mol2)==3) call outfold ( Mol2, stion2, endion2 )

    particle_positions(1,Mol1) = particle_positions(1,Mol1) + x12
    particle_positions(2,Mol1) = particle_positions(2,Mol1) + y12
    particle_positions(3,Mol1) = particle_positions(3,Mol1) + z12
    
    particle_positions(1,Mol2) = particle_positions(1,Mol2) - x12
    particle_positions(2,Mol2) = particle_positions(2,Mol2) - y12
    particle_positions(3,Mol2) = particle_positions(3,Mol2) - z12
    
    do i = stion1, endion1 
        
        charge_positions(1,i) = charge_positions(1,i) + x12
        charge_positions(2,i) = charge_positions(2,i) + y12
        charge_positions(3,i) = charge_positions(3,i) + z12
        
    end do 
    
    do i = stion2, endion2  
    
        charge_positions(1,i) = charge_positions(1,i) - x12
        charge_positions(2,i) = charge_positions(2,i) - y12
        charge_positions(3,i) = charge_positions(3,i) - z12

    end do
     
    ! periodic boundaries:
    do ii = stion1, endion1
      do i = 1,3  
        if (charge_positions(i,ii) > boxlength) charge_positions(i,ii)= charge_positions(i,ii) - boxlength
        if (charge_positions(i,ii) < 0.0_rk)    charge_positions(i,ii)= charge_positions(i,ii) + boxlength  
      end do  
    end do 
    do ii = stion2, endion2
      do i = 1,3  
        if (charge_positions(i,ii) > boxlength) charge_positions(i,ii)= charge_positions(i,ii) - boxlength
        if (charge_positions(i,ii) < 0.0_rk)    charge_positions(i,ii)= charge_positions(i,ii) + boxlength  
      end do  
    end do 
     
    if (lenion1 .GT. 0) then
        call RealMolecule_move ( stion1, endion1, UREAL1_new)                     
    else
        UREAL1_new = 0.
    end if


    if (lenion2 .GT. 0) then
        call RealMolecule_move ( stion2, endion2, UREAL2_new)  
    else
        UREAL2_new = 0.
    end if
                                            
                                            
    if (lenion1 .GT. 0) then 

        call Fourier_Move2( lenion1, charge_positions(1,stion1:endion1), charge_positions(2,stion1:endion1), &
                                                charge_positions(3,stion1:endion1), &
                                                charge_list(stion1:endion1),  &
                                                EXPX(:,stion1:endion1), &
                                                EXPY(:,stion1:endion1), EXPZ(:,stion1:endion1), EXPX1_NEW(:,1:lenion1), &
                                                EXPY1_NEW(:,1:lenion1), EXPZ1_NEW(:,1:lenion1), &
                                                SUMQEXPV, SUMQEXPV_TMP, dUFOURIER1 )
    
                EXPX(:,stion1:endion1) = EXPX1_NEW(:,1:lenion1)
                EXPY(:,stion1:endion1) = EXPY1_NEW(:,1:lenion1)
                EXPZ(:,stion1:endion1) = EXPZ1_NEW(:,1:lenion1)        
!                 SUMQEXPV = SUMQEXPV_NEW      

    else

        dUFOURIER1 = 0.
        SUMQEXPV_TMP = SUMQEXPV !_NEW
        
    end if


    if (lenion2 .GT. 0) then
    
        call Fourier_Move2( lenion2, charge_positions(1,stion2:endion2), charge_positions(2,stion2:endion2), &
                                                charge_positions(3,stion2:endion2), &
                                                charge_list(stion2:endion2),  &
                                                EXPX(:,stion2:endion2), &
                                                EXPY(:,stion2:endion2), EXPZ(:,stion2:endion2), EXPX2_NEW(:,1:lenion2), &
                                                EXPY2_NEW(:,1:lenion2), EXPZ2_NEW(:,1:lenion2), &
                                                SUMQEXPV_TMP, SUMQEXPV, dUFOURIER2 )
    
                EXPX(:,stion2:endion2) = EXPX2_NEW(:,1:lenion2)
                EXPY(:,stion2:endion2) = EXPY2_NEW(:,1:lenion2)
                EXPZ(:,stion2:endion2) = EXPZ2_NEW(:,1:lenion2)
!                 SUMQEXPV = SUMQEXPV_NEW      
     
    else

        dUFOURIER2 = 0.

    end if   

    dUFOURIER(k) = (dUFOURIER1 + dUFOURIER2) !* CoulCombo    
    dUREAL(k) = (UREAL1_new + UREAL2_new - UREAL1_old - UREAL2_old) !* CoulCombo
    U_corr(k) = 2.0_rk*U_corr1 !* CoulCombo
   
    deallocate( EXPX1_NEW)
    deallocate( EXPY1_NEW)
    deallocate( EXPZ1_NEW)

    deallocate( EXPX2_NEW)
    deallocate( EXPY2_NEW)
    deallocate( EXPZ2_NEW) 
    
end do kloop        
 
dUREAL_sum = sum(dUREAL(1:NSwap))
dUFOURIER_sum = sum(dUFOURIER(1:NSwap))

ener_delt = dUREAL_sum + dUFOURIER_sum + sum(U_corr)
   
! accept or reject move
    call random_number(acceptance)
    if (log(acceptance) <  - lambda * ener_delt) then
        total_accepted_swap = total_accepted_swap + 1_ik 
!         print*, 'Accepted', ener_delt  !; pause
    else
        ! set back    
        charge_positions = charge_position_old 
        particle_positions = position_old 
        EXPX = EXPX_OLD
        EXPY = EXPY_OLD
        EXPZ = EXPZ_OLD 
        SUMQEXPV  = SUMQEXPV_OLD
        ener_delt = 0.0_rk
        dUREAL_sum = 0.0_rk
        dUFOURIER_sum = 0.0_rk
!        print*, 'Rejected'!; pause
    end if 
    
    deallocate(position_old) 
    deallocate(charge_position_old)
    deallocate(dUREAL)
    deallocate(dUFOURIER)
    deallocate(U_corr)

return

end subroutine particle_swap_mult

end module particle_swap_mult_mod 


