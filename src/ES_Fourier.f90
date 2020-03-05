!> \brief fourier contribution to total energy  
!> \author M.Theiss 
 module ES_Fourier_mod
  use working_prec_mod, only: rk, ik
  use basic_parameters_mod
  use sim_parameters_mod 
  use globals_mod 
  implicit none
  private 
  
  public :: Fourier_Setup  
  public :: Fourier_Move   
  public :: Fourier_Move2
 
contains

!---------------------------------------------------------------------------------------------------------------------------------------
!> fourier setup  
!---------------------------------------------------------------------------------------------------------------------------------------
subroutine Fourier_Setup(  ) 
integer(ik):: i, j, k, Ksq 
real(rk):: b

b = PI * PI / ( Alpha * Alpha )

Nkvec = 0_ik

do i = 0, Kmax
    do j = -Kmax, Kmax
        do k = -Kmax, Kmax
        
            Ksq = i*i + j*j + k*k
            
                        if( ( Ksq /= 0_ik ) .AND. ( Ksq < Kmax * Kmax ) ) then 
                            Nkvec = Nkvec + 1_ik 
                        end if
                        
            end do
    end do
end do

allocate( CONST(Nkvec) ) 
allocate( KX(Nkvec) ) 
allocate( KY(Nkvec) ) 
allocate( KZ(Nkvec) ) 

 

Nkvec = 0_ik

do i = 0, Kmax
    do j = -Kmax, Kmax
            do k = -Kmax, Kmax
                    
                    Ksq = i*i + j*j + k*k
                    
                    if( ( Ksq /= 0 ) .AND. ( Ksq < Kmax * Kmax ) ) then
                    
                            Nkvec = Nkvec + 1_ik

                            CONST(Nkvec) = 1.0_rk / ( Pi * real( i*i + j*j + k*k, rk ) ) * &
                                                        exp( -b * real( i*i + j*j + k*k, rk ) )

                            if( i == 0 ) CONST(Nkvec) = 0.5_rk * CONST(Nkvec)

                            KX(Nkvec) = i
                            KY(Nkvec) = j
                            KZ(Nkvec) = k

                    end if
            end do
    end do
end do

end subroutine Fourier_Setup

 
!---------------------------------------------------------------------------------------------------------------------------------------
!> fourier move  
!---------------------------------------------------------------------------------------------------------------------------------------
subroutine Fourier_Move  ( i, EXPX_CH, EXPY_CH, EXPZ_CH, EXPX_NEW, &
                                                    EXPY_NEW, EXPZ_NEW, &
                                                    SUMQEXPV_in, SUMQEXPV_NEW, DU_FOURIER )

integer(ik), intent(in) :: i 

complex(rk), dimension(0:Kmax), intent(in)      :: EXPX_CH
complex(rk), dimension(-Kmax:Kmax), intent(in)  :: EXPY_CH, EXPZ_CH
complex(rk), dimension(0:Kmax), intent(out)     :: EXPX_NEW
complex(rk), dimension(-Kmax:Kmax), intent(out) :: EXPY_NEW, EXPZ_NEW

! SUMQEXPV contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian before a move.
! SUMQEXPV_NEW contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian after a move.

complex(rk), dimension(Nkvec), intent(in)   :: SUMQEXPV_in
complex(rk), dimension(Nkvec), intent(out)  :: SUMQEXPV_NEW
 
real(rk), intent(out)   :: DU_FOURIER

! Local variables

integer(ik) :: k
integer(ik) :: kkx, kky, kkz 
real(rk)    :: xi, yi, zi
real(rk)    :: bx, by, bz 
complex(rk) :: expv_diff



DU_FOURIER = 0.0_rk

EXPX_NEW = (1.0_rk,0.0_rk)
EXPY_NEW = (1.0_rk,0.0_rk)
EXPZ_NEW = (1.0_rk,0.0_rk)
 
! do i = 1, Ns

    xi = charge_positions(1,i)  !Xn(i)
    yi = charge_positions(2,i)  !Yn(i)
    zi = charge_positions(3,i)  !Zn(i)
    
    do k = 1, Kmax

            bx = real(k,rk) * 2.0_rk * PI / boxlength * xi
            by = real(k,rk) * 2.0_rk * PI / boxlength * yi
            bz = real(k,rk) * 2.0_rk * PI / boxlength * zi
            
            EXPX_NEW( k ) = cmplx( cos( bx ), sin( bx ),rk )
            EXPY_NEW( k ) = cmplx( cos( by ), sin( by ),rk )
            EXPZ_NEW( k ) = cmplx( cos( bz ), sin( bz ),rk )
            EXPY_NEW(-k ) = conjg( EXPY_NEW( k ) ) 
            EXPZ_NEW(-k ) = conjg( EXPZ_NEW( k ) )

    end do

! end do

SUMQEXPV_NEW = SUMQEXPV_in

! do i = 1, Ns
 

    do k = 1, Nkvec

            kkx = KX(k)
            kky = KY(k)
            kkz = KZ(k) 

            expv_diff = EXPX_NEW( kkx ) * EXPY_NEW( kky ) *  &
                                    EXPZ_NEW( kkz ) - &
                                    EXPX_CH( kkx ) * EXPY_CH( kky ) * &
                                    EXPZ_CH( kkz )

            SUMQEXPV_NEW( k ) = SUMQEXPV_NEW( k ) +  expv_diff * charge_list(i)   

    end do

! end do

do k = 1, Nkvec

    DU_FOURIER = DU_FOURIER + CONST(k) / boxlength *   &
                        ( conjg( SUMQEXPV_NEW( k ) ) * SUMQEXPV_NEW( k ) - &
                          conjg( SUMQEXPV_in( k ) ) * SUMQEXPV_in( k ) )

end do



end subroutine Fourier_Move



!---------------------------------------------------------------------------------------------------------------------------------------
!> fourier move #2 (faster)  
!---------------------------------------------------------------------------------------------------------------------------------------
subroutine Fourier_Move2 ( Ns, Xn, Yn, Zn, CHARGEc,  &
                                    EXPX_in, EXPY_in, EXPZ_in, EXPX_NEW, EXPY_NEW, EXPZ_NEW, &
                                    SUMQEXPV_in, SUMQEXPV_NEW, DU_FOURIER )
integer(ik), intent(in):: Ns
real(rk), dimension(Ns), intent(in):: Xn, Yn, Zn
real(rk), dimension(Ns), intent(in):: CHARGEc

complex(rk), dimension(0:Kmax, Ns), intent(in) :: EXPX_in
complex(rk), dimension(-Kmax:Kmax, Ns), intent(in) :: EXPY_in, EXPZ_in
complex(rk), dimension(0:Kmax, Ns), intent(out):: EXPX_NEW
complex(rk), dimension(-Kmax:Kmax, Ns), intent(out):: EXPY_NEW, EXPZ_NEW

! SUMQEXPV contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian before a move.
! SUMQEXPV_NEW contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian after a move.

complex(rk), dimension(Nkvec), intent(in) :: SUMQEXPV_in
complex(rk), dimension(Nkvec), intent(out):: SUMQEXPV_NEW

! DU_FOURIER is the change in fourier energy for a given move.

real(rk), intent(out) :: DU_FOURIER

! Local variables
 
integer(ik) :: i, k
integer(ik) :: kkx, kky, kkz 
real(rk)    :: xi, yi, zi
real(rk)    :: bx, by, bz 
complex(rk) :: expv_diff

DU_FOURIER = 0.0_rk

EXPX_NEW = ( 1.0_rk, 0.0_rk )
EXPY_NEW = ( 1.0_rk, 0.0_rk )
EXPZ_NEW = ( 1.0_rk, 0.0_rk )

do i = 1, Ns

    xi = Xn(i)
    yi = Yn(i)
    zi = Zn(i)
    
    do k =1, Kmax

            bx = real(k,rk) * 2.0_rk * Pi / boxlength * xi
            by = real(k,rk) * 2.0_rk * Pi / boxlength * yi
            bz = real(k,rk) * 2.0_rk * Pi / boxlength * zi
             
            EXPX_NEW( k, i ) = cmplx( cos( bx ), sin( bx ),rk )
            EXPY_NEW( k, i ) = cmplx( cos( by ), sin( by ),rk )
            EXPZ_NEW( k, i ) = cmplx( cos( bz ), sin( bz ),rk )
            EXPY_NEW(-k, i ) = conjg( EXPY_NEW( k, i ) ) 
            EXPZ_NEW(-k, i ) = conjg( EXPZ_NEW( k, i ) )

    end do

end do

SUMQEXPV_NEW = SUMQEXPV_in

do i = 1, Ns
 

    do k = 1, Nkvec

            kkx = KX(k)
            kky = KY(k)
            kkz = KZ(k)

            expv_diff = EXPX_NEW( kkx, i ) * EXPY_NEW( kky, i ) *  &
                                    EXPZ_NEW( kkz, i ) - &
                                    EXPX_in( kkx, i ) * EXPY_in( kky, i ) * &
                                    EXPZ_in( kkz, i )

            SUMQEXPV_NEW( k ) = SUMQEXPV_NEW( k ) + expv_diff * CHARGEc( i)

    end do

end do

do k = 1, Nkvec

    DU_FOURIER = DU_FOURIER + CONST(k) / boxlength *   &
                                ( conjg( SUMQEXPV_NEW( k ) ) * SUMQEXPV_NEW( k ) - &
                                  conjg( SUMQEXPV_in( k ) ) * SUMQEXPV_in( k ) )
 
end do

return

end subroutine Fourier_Move2


end module ES_Fourier_mod



