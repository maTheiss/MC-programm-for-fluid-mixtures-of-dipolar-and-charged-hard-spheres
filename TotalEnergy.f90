 module TotalEnergy_mod
  use working_prec_mod, only: rk, ik
  use globals_mod
  implicit none
  private

  public :: TotalEnergy  
 
contains
 
  subroutine TotalEnergy( Kmax, Nkvec, KX, KY, KZ, CONST, &
				   EXPX, EXPY, EXPZ, SUMQEXPV, SUMQX, SUMQY, SUMQZ, &
                                   UREAL, USCH,  USDIP, UFOURIER, USURF) 
                                   
! Kmax is an Ewald sum parameter.
! Nkvec is the number of k-vectors used in the Fourier sum.
! KX, KY, KZ contain the vector identity of the Nkvec vectors.
! CONST contains the constant part of the Fourier summation for a given Nkvec.

integer, intent(in)										:: Kmax
integer, intent(in)										:: Nkvec
integer, dimension(Nkvec), intent(in)					:: KX, KY, KZ
real, dimension(Nkvec), intent(in)						:: CONST

! EXPX contains the value of exp( i*kx*x ) for a given kx and ion.

complex, dimension(0:Kmax, ncharges), intent(out) 			:: EXPX
complex, dimension(-Kmax:Kmax, ncharges), intent(out) 		:: EXPY, EXPZ

! SUMQEXPV contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian.

complex, dimension(Nkvec), intent(inout)				 :: SUMQEXPV

! SUMQX is the summation of qi * xi for all ions in the box.

real, intent(inout)									:: SUMQX, SUMQY, SUMQZ


    real(rk),intent(out):: UREAL
    real(rk),intent(out):: USCH
    real(rk),intent(out):: USDIP
    real(rk),intent(out):: UFOURIER
    real(rk),intent(out):: USURF
    
   CoulCombo = ec * ec * 1.0_rke10_rk / ( 4.0_rk * PI * eps0 * kB )
   
   USURF = 0._rk
   ! real-term 
    call RealMolecule(  UREAL  )
    UREAL = UREAL * CoulCombo
    
    ! Self-term: dipole
    call SelfMolecule( USDIP ) 
    USDIP = USDIP * CoulCombo 
    
    ! Fourier-term 
!     call FourierMolecule ( UFOURIER ) 
    UFOURIER = 0.0_rk
    EXPX = ( 0.0_rk, 0.0_rk )
    EXPY = ( 0.0_rk, 0.0_rk )
    EXPZ = ( 0.0_rk, 0.0_rk )

    SUMQEXPV = ( 0.0_rk )

! !     ii = 0_ik 
    maxcharges = ncharges  
! !     if ( ndpioles /= 0 ) maxcharges = ncharges - 1_ik    ! -1 nur wenn auch Dipole enthalten sind
    do i=1, maxcharges       
    
! !             if ( i .le. nions ) then 
! !                 LENGTHion = 1_ik 
! !             else 
! !                 LENGTHion = 2_ik
! !             end if
            
            
!             stion = STARTion(i)
!             endion = STARTion(i) + LENGTHion(i)-1
                                                           !starting ionic bead of each molecule   =! 1 ??
                                                ! LENGTHion = 1(ion) or 2(dipole) or altogether 
                                                                               ! LENGTHion= contains the number of ionic beads in each molecule.
                                                                                                            !=Ns is the number of ions being displaced, rotated, or transferred.
                                                                                                            ! For a volume change Ns is the number of ions in the simulation box.
                                                                                                            !Xion(stion:endion) entspr. 1 oder 2 Positionen !?
! !             stion = 1_ik
! !              if ( i .le. nions ) then 
! !                 endion = 1_ik 
! !                 Xion(stion:endion) = charge_positions(1,i)
! !                 Yion(stion:endion) = charge_positions(2,i)
! !                 Zion(stion:endion) = charge_positions(3,i)
! !             else 
! !                 endion = 2_ik 
! !                 Xion(stion:endion) = charge_positions(1,i+ii:i+ii+1)
! !                 Yion(stion:endion) = charge_positions(2,i+ii:i+ii+1)
! !                 Zion(stion:endion) = charge_positions(3,i+ii:i+ii+1)
! !                 ii = ii + 1_ik
! !             end if
! !             LENGTHion = endion
! !             CHARGE = charge_list(i)
            
            EXPX_CH(:) = EXPX(:,i) 
            EXPY_CH(:) = EXPY(:,i)
            EXPZ_CH(:) = EXPZ(:,i)
            EXPX_NEW_CH(:) = EXPX_NEW(:,i)
            EXPY_NEW_CH(:) = EXPY_NEW(:,i)
            EXPZ_NEW_CH(:) = EXPZ_NEW(:,i)
            
            call Fourier_Move( i,    Kmax, Nkvec, KX, KY, KZ, CONST, EXPX_CH(:), &
                                                    EXPY_CH(:), EXPZ_CH(:), EXPX_NEW_CH(:), &
                                                    EXPY_NEW_CH(:), EXPZ_NEW_CH(:), &
                                                    SUMQEXPV, SUMQEXPV_NEW, dUFOURIER )
                
!             call Fourier_Move( LENGTHion, Xion(stion:endion), Yion(stion:endion), &
!                                             Zion(stion:endion), &
!                                             TYPEion(stion:endion), Niongrs, CHARGE,  &
!                                                     boxlength, Kmax, Nkvec, KX, KY, KZ, CONST, EXPX(:,stion:endion), &
!                                                     EXPY(:,stion:endion), EXPZ(:,stion:endion), EXPX_NEW(:,stion:endion), &
!                                                     EXPY_NEW(:,stion:endion), EXPZ_NEW(:,stion:endion), &
!                                                     SUMQEXPV, SUMQEXPV_NEW, dUFOURIER )

                    UFOURIER = UFOURIER + dUFOURIER
                    EXPX(:,i) = EXPX_NEW_CH(:)
                    EXPY(:,i) = EXPY_NEW_CH(:)
                    EXPZ(:,i) = EXPZ_NEW_CH(:)
                    SUMQEXPV = SUMQEXPV_NEW

    end do
    UFOURIER = UFOURIER * CoulCombo
    
    
!     call SurfaceMolecule ( USURF_in ) 
     
    
    ! Self-term: charges
    USCH = 0._rk
    
    do i=1, nions
        USCH = USCH + charge_list(i)** 2._rk
    end do
    
    USCH = Alpha / sqrt(PI) / boxlength * USCH
    USCH = USCH * CoulCombo
    
  end subroutine TotalEnergy

end subroutine TotalEnergy_mod
