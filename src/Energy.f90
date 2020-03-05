!> \brief total energy  
!> \author M.Theiss 
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
 
  
!---------------------------------------------------------------------------------------------------------------------------------------
!> total energy 
!---------------------------------------------------------------------------------------------------------------------------------------
subroutine TotalEnergy( UREAL, USCH,  USDIP, UFOURIER, USURF, Utot)  
  
real(rk),intent(out):: UREAL
real(rk),intent(out):: USCH
real(rk),intent(out):: USDIP
real(rk),intent(out):: UFOURIER
real(rk),intent(out):: USURF, Utot
  
complex(rk), dimension(Nkvec)           :: SUMQEXPV_NEW 
complex(rk), allocatable, dimension(:,:):: EXPX_NEW
complex(rk), allocatable, dimension(:,:):: EXPY_NEW, EXPZ_NEW
real(rk), dimension(npart):: USELF_MOL
 
real(rk):: dUFOURIER
integer(ik) ::Nc, i, stion, endion, chargeid(2)  !LENGTH 
     
   
    ! real-term 
    call RealInteract( UREAL )
        
    ! Fourier-term  
    UFOURIER = 0.0_rk
    
    EXPX = (0.0_rk,0.0_rk)
    EXPY = (0.0_rk,0.0_rk)
    EXPZ = (0.0_rk,0.0_rk)
    SUMQEXPV = ( 0.0_rk )
    
    do i=1, npart

        chargeid = return_charge_index (i) 
         
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

    ! Surface-term
    USURF = 0.0_rk !hat keinen Einfluss (bei unendlich gr. Boxlänge)

    ! Self-term: dipole
    USDIP = 0.0_rk
    if ( ndipoles /= 0_ik)  call SelfMolecule( USDIP ) 
   
    
    ! Self-term: charges (ions and dipoles)
    ! summation of the square of all the charges
    USCH = 0.0_rk 
    do i=1, ncharges  
        USCH = USCH + charge_list(i)** 2.0_rk 
    end do 
    USCH = Alpha / sqrt(PI) / boxlength * USCH  
!     USCH = USCH * CoulCombo 
 
    Utot = UREAL - USCH -  USDIP + UFOURIER + USURF

end subroutine TotalEnergy
  
  
!---------------------------------------------------------------------------------------------------------------------------------------
!>  for delta energy
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

real(rk),intent(out):: UREAL
real(rk),intent(out):: USCH
real(rk),intent(out):: USDIP
real(rk),intent(out):: UFOURIER
real(rk),intent(out):: USURF
real(rk):: dUFOURIER
integer(ik) :: i 


   ! real-term 
    UREAL = 0.0_rk
    call RealMolecule_move(  stion, endion,  UREAL)
    
    ! Fourier-term  
    UFOURIER = 0.0_rk    
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
             
    USURF = 0.0_rk
     
    ! Self-term: dipole
    ! not necessary for ener_delt 
    USDIP = 0.0_rk 
    
    ! Self-term: charges
    ! not necessary for ener_delt 
    USCH = 0.0_rk 
    
end subroutine DispRot_Energy

end module Energy_mod
