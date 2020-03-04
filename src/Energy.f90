 module Energy_mod
  use working_prec_mod, only: rk, ik
  use basic_parameters_mod
  use sim_parameters_mod 
  use SelfMolecule_mod 
  use globals_mod
  use ES_Fourier_mod
  use RealMolecule_mod
  implicit none
  private

  public :: TotalEnergy  
  public :: DispRot_Energy
 
contains
 
  subroutine TotalEnergy( UREAL, USCH,  USDIP, UFOURIER, USURF, Utot)  
  
real(rk),intent(out):: UREAL
real(rk),intent(out):: USCH
real(rk),intent(out):: USDIP
real(rk),intent(out):: UFOURIER
real(rk),intent(out):: USURF, Utot
  
complex(rk), dimension(Nkvec)           :: SUMQEXPV_NEW 
complex(rk), allocatable, dimension(:,:):: EXPX_NEW
complex(rk), allocatable, dimension(:,:):: EXPY_NEW, EXPZ_NEW
! complex(rk), dimension(0:Kmax, ncharges):: EXPX_NEW
! complex(rk), dimension(-Kmax:Kmax, ncharges):: EXPY_NEW, EXPZ_NEW

real(rk), dimension(npart):: USELF_MOL
 
real(rk):: dUFOURIER
integer(ik) ::Nc, i, stion, endion, chargeid(2)  !LENGTH 
     
   
   ! real-term 
!     call RealMolecule(  UREAL  )   
!     UREAL = UREAL * CoulCombo 
!     print*, '------'
!     print*, UREAL/real(npart,rk)/Temp 
!     UREAL = 0.0_rk
    call RealInteract( UREAL )
!     UREAL = UREAL * CoulCombo
!     print*, UREAL/real(npart,rk)/Temp
!     print*, '------'
!     stop
     
        
    !----------
    ! Fourier-term  
    UFOURIER = 0.0_rk
    
    EXPX = (0.0_rk,0.0_rk)
    EXPY = (0.0_rk,0.0_rk)
    EXPZ = (0.0_rk,0.0_rk)
    SUMQEXPV = ( 0.0_rk )
    
    do i=1, npart

        chargeid = return_charge_index (i) 
!         
        Nc = LENGTH(i) 
        
        stion = chargeid(1)    
        endion =  stion +LENGTH(i)-1_ik 
        
        allocate(EXPX_NEW(0:Kmax,1:Nc))
        allocate(EXPY_NEW(-Kmax:Kmax,1:Nc))
        allocate(EXPZ_NEW(-Kmax:Kmax,1:Nc))
        
        call Fourier_Move2( Nc, charge_positions(1,stion:endion), charge_positions(2,stion:endion), &
                                                charge_positions(3,stion:endion), &
                                                charge_list(stion:endion),  &
                                                EXPX(:,stion:endion), &
                                                EXPY(:,stion:endion), EXPZ(:,stion:endion), EXPX_NEW(:,1:Nc), &
                                                EXPY_NEW(:,1:Nc), EXPZ_NEW(:,1:Nc), &
                                                SUMQEXPV, SUMQEXPV_NEW, dUFOURIER )

                UFOURIER = UFOURIER + dUFOURIER
                EXPX(:,stion:endion) = EXPX_NEW(:,1:Nc)
                EXPY(:,stion:endion) = EXPY_NEW(:,1:Nc)
                EXPZ(:,stion:endion) = EXPZ_NEW(:,1:Nc)
                SUMQEXPV = SUMQEXPV_NEW
        
        deallocate(EXPX_NEW)
        deallocate(EXPY_NEW)
        deallocate(EXPZ_NEW)
        
    end do
!     UFOURIER = UFOURIER !* CoulCombo 
!         print*,'Fourier1',UFOURIER/real(npart,rk)/Temp
        
    !----------
    ! Fourier-term  
! !     UFOURIER = 0.0_rk
! !     EXPX = (0.0_rk,0.0_rk)
! !     EXPY = (0.0_rk,0.0_rk)
! !     EXPZ = (0.0_rk,0.0_rk)
! ! 
! !     SUMQEXPV = ( 0.0_rk )
! !     
! !     do i=1, ncharges       
! !      
! !             
! !             call Fourier_Move( i,    EXPX(:,i), EXPY(:,i), EXPZ(:,i), EXPX_NEW(:,i) , &
! !                                                     EXPY_NEW(:,i) , EXPZ_NEW(:,i) , &
! !                                                     SUMQEXPV, SUMQEXPV_NEW, dUFOURIER )
! !                  
! !                     UFOURIER = UFOURIER + dUFOURIER
! !                     EXPX(:,i) = EXPX_NEW(:,i) 
! !                     EXPY(:,i) = EXPY_NEW(:,i) 
! !                     EXPZ(:,i) = EXPZ_NEW(:,i) 
! !                     SUMQEXPV = SUMQEXPV_NEW
! ! 
! !     end do 
! !     UFOURIER = UFOURIER * CoulCombo
! !     print*,'Fourier2',UFOURIER/real(npart,rk)/Temp


    ! Surface-term
    USURF = 0.0_rk !hat keinen Einfluss (bei unendlich gr. Boxlänge)

    ! Self-term: dipole
    USDIP = 0.0_rk
    if ( ndipoles /= 0_ik)  call SelfMolecule( USDIP ) 
!     USDIP = USDIP * CoulCombo 
! !     !--------
! !     print*,'1',USDIP/real(npart,rk) 
! !     
! !     USELF_MOL(:) = 0.0_rk
! !     do i=1, npart
! !     
! !         chargeid = return_charge_index (i)
! !          
! !         Nc = LENGTH(i)
! !             
! !             stion = chargeid(1)   !STARTion(i)
! !             endion =  stion +LENGTH(i)-1_ik  !STARTion(i)+LENGTH-1
! ! 
! !             if ( LENGTH(i)  .GT. 1 ) then 
! !             
! !                     call SelfMolecule2(Nc, charge_positions(1,stion:endion), charge_positions(2,stion:endion), &
! !                                                     charge_positions(3,stion:endion), &
! !                                                     charge_list(stion:endion), &
! !                                                     USELF_MOL(i))
! ! 
! ! !                     USELF_MOL(i) = USELF_MOL(i) * CoulCombo
! ! 
! !         else
! ! 
! !             USELF_MOL(i) = 0.0
! ! 
! !         end if
! ! 
! !     end do   
! !     USDIP = sum(USELF_MOL)  
! !     print*,'2',USDIP/real(npart,rk)
! !     !-------- 
    
    
    ! Self-term: charges (ions and dipoles)
    ! summation of the square of all the charges
    USCH = 0.0_rk 
    do i=1, ncharges  
        USCH = USCH + charge_list(i)** 2.0_rk 
    end do 
    USCH = Alpha / sqrt(PI) / boxlength * USCH  
!     USCH = USCH * CoulCombo 
 
    Utot = UREAL - USCH -  USDIP + UFOURIER + USURF
!     
!     print*,'Utot'!, Utot/real(npart,rk)/Temp
!     print*,UREAL/real(npart,rk)/Temp, UFOURIER/real(npart,rk)/Temp
!     print*,USCH/real(npart,rk)/Temp, USDIP/real(npart,rk)/Temp
  end subroutine TotalEnergy
  
  
 !---------------------------------------------------------------------------------------------------------------------------------------
 !  for delta energy
!---------------------------------------------------------------------------------------------------------------------------------------
  subroutine DispRot_Energy( pid, Nc, stion, endion, UREAL, USCH,  USDIP, UFOURIER, USURF) 
                             
integer(ik), intent(inout):: pid  ! particle-ID
integer(ik), intent(inout):: Nc   != LENGTH(pid) 
integer(ik), intent(inout):: stion
integer(ik), intent(inout):: endion   
  
complex(rk), dimension(0:Kmax,Nc):: EXPX_NEW
complex(rk), dimension(-Kmax:Kmax,Nc):: EXPY_NEW, EXPZ_NEW
complex(rk), dimension(0:Kmax):: EXPX_NEW_CH
complex(rk), dimension(-Kmax:Kmax):: EXPY_NEW_CH, EXPZ_NEW_CH
complex(rk), dimension(Nkvec):: SUMQEXPV_NEW

! complex(rk), dimension(0:Kmax, stion:endion):: EXPX_old
! complex(rk), dimension(-Kmax:Kmax, stion:endion):: EXPY_old, EXPZ_old
! complex(rk),dimension(Nkvec)                 :: SUMQEXPV_old

real(rk),intent(out):: UREAL
real(rk),intent(out):: USCH
real(rk),intent(out):: USDIP
real(rk),intent(out):: UFOURIER
real(rk),intent(out):: USURF
real(rk):: dUFOURIER
! logical :: isoverlapping

! integer(ik) :: chargeid(2)
integer(ik) :: i 

! EXPX_old(:,stion:endion) = EXPX(:,stion:endion)
! EXPY_old(:,stion:endion) = EXPY(:,stion:endion)
! EXPZ_old(:,stion:endion) = EXPZ(:,stion:endion)
! SUMQEXPV_old = SUMQEXPV 

    UREAL = 0.0_rk
   ! real-term 
    call RealMolecule_move(  stion, endion,  UREAL)
!     UREAL = UREAL !* CoulCombo
!     print*, '------'
!     print*, UREAL/real(npart,rk)/Temp
!     UREAL = 0.0_rk 
!     call RealInteract_move(  pid,  UREAL) 
!     UREAL = UREAL * CoulCombo
!     print*, UREAL/real(npart,rk)/Temp
!     print*, '------'
!     stop
     
    ! Fourier-term  
!     UFOURIER = 0.0_rk       
!     i = stion 
!             
!             call Fourier_Move( i, EXPX(:,i), EXPY(:,i), EXPZ(:,i), EXPX_NEW_CH(:), &
!                                                     EXPY_NEW_CH(:), EXPZ_NEW_CH(:), &
!                                                     SUMQEXPV, SUMQEXPV_NEW, dUFOURIER )
!                  
!                     UFOURIER = UFOURIER + dUFOURIER 
!                     EXPX(:,i) = EXPX_NEW_CH(:)
!                     EXPY(:,i) = EXPY_NEW_CH(:)
!                     EXPZ(:,i) = EXPZ_NEW_CH(:)
!                     SUMQEXPV = SUMQEXPV_NEW
!     if ( endion /= stion ) then 
!     i = endion       
!             call Fourier_Move( i, EXPX(:,i), EXPY(:,i), EXPZ(:,i), EXPX_NEW_CH(:), &
!                                                     EXPY_NEW_CH(:), EXPZ_NEW_CH(:), &
!                                                     SUMQEXPV, SUMQEXPV_NEW, dUFOURIER )
!                  
!                     UFOURIER = UFOURIER + dUFOURIER 
!                     EXPX(:,i) = EXPX_NEW_CH(:)
!                     EXPY(:,i) = EXPY_NEW_CH(:)
!                     EXPZ(:,i) = EXPZ_NEW_CH(:)
!                     SUMQEXPV = SUMQEXPV_NEW
!     end if  
!     UFOURIER = UFOURIER * CoulCombo
!     print*,'-------------------'
!     print*, UFOURIER/real(npart,rk)/Temp
    
    UFOURIER = 0.0_rk    
     
    
        !zum Vgl. mit Fourier_Move 
!     EXPX(:,stion:endion) = EXPX_old(:,stion:endion)
!     EXPY(:,stion:endion) = EXPY_old(:,stion:endion)
!     EXPZ(:,stion:endion) = EXPZ_old(:,stion:endion)
!     SUMQEXPV = SUMQEXPV_old  
    
    call Fourier_Move2( Nc, charge_positions(1,stion:endion), charge_positions(2,stion:endion), &
                                            charge_positions(3,stion:endion), &
                                            charge_list(stion:endion),  &
                                            EXPX(:,stion:endion), &
                                            EXPY(:,stion:endion), EXPZ(:,stion:endion), EXPX_NEW(:,1:Nc), &
                                            EXPY_NEW(:,1:Nc), EXPZ_NEW(:,1:Nc), &
                                            SUMQEXPV, SUMQEXPV_NEW, dUFOURIER )

            UFOURIER = dUFOURIER  !UFOURIER + 
            EXPX(:,stion:endion) = EXPX_NEW(:,1:Nc)
            EXPY(:,stion:endion) = EXPY_NEW(:,1:Nc)
            EXPZ(:,stion:endion) = EXPZ_NEW(:,1:Nc)
            SUMQEXPV = SUMQEXPV_NEW
             
     
!     UFOURIER = UFOURIER !* CoulCombo 
    
!     print*,'DR',UREAL/real(npart,rk)/Temp, UFOURIER/real(npart,rk)/Temp
     
    
    USURF = 0.0_rk
!     call SurfaceMolecule ( USURF_in ) 
     
    ! Self-term: dipole
    ! not necessary for ener_delt (not dependend  of position) 
    USDIP = 0.0_rk 
    
    ! Self-term: charges
    ! not necessary for ener_delt (not dependend of position)
    USCH = 0.0_rk 
    
  end subroutine DispRot_Energy

end module Energy_mod
