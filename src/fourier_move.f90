
subroutine Fourier_Move( Ns, Xn, Yn, Zn, TYPEs, Niongrs, CHARGE,  &
						 BoxSize, Kmax, Nkvec, KX, KY, KZ, CONST, EXPX, &
						 EXPY, EXPZ, EXPX_NEW, EXPY_NEW, EXPZ_NEW, &
						 SUMQEXPV, SUMQEXPV_NEW, DU_FOURIER )

implicit none

! This routine calculates the change in the Fourier term of the 
! Ewald sum for a displacement, rotation, tranfer, or volume change.

! ************************  PLEASE NOTE **********************************!
! For a particle transfer, EXPX = EXPY = EXPZ = ( 0.0, 0.0 )
! For the energy change of the DONOR phase after a particle 
! transfer ... CHARGE must be entered as the NEGATIVE of CHARGE, -CHARGE.
! ************************************************************************!

! Ns is the number of ions being displaced, rotated, or transferred.
! For a volume change Ns is the number of ions in the simulation box.

! Xn, Yn, Zn are the new positions of the ions after a displacement, 
! rotation, or volume change.  For a transfer Xn, Yn, Zn are the new
! positions in the receiver phase or the old positions in the donor phase.
! For a volume change you recalculate the energy from scatch.
! TYPEs contains the group identity of the Ns ions.

integer, intent(in)										:: Ns
real, dimension(Ns), intent(in)							:: Xn, Yn, Zn
integer, dimension(Ns), intent(in)						:: TYPEs

! Niongrs is the number of ion groups.

integer, intent(in)										:: Niongrs

! CHARGE contains the charge for a given group and hamiltonian.

real, dimension(Niongrs), intent(in)						:: CHARGE

! Boxsize is the length of the simulation box.

real, intent(in)											:: BoxSize

! Kmax is an Ewald sum parameter.
! Nkvec is the number of k-vectors used in the Fourier sum.
! KX, KY, KZ contain the vector identity of the Nkvec vectors.
! CONST contains the constant part of the Fourier summation for a given Nkvec.

integer, intent(in)										:: Kmax
integer, intent(in)										:: Nkvec
integer, dimension(Nkvec), intent(in)					:: KX, KY, KZ
real, dimension(Nkvec), intent(in)						:: CONST

! EXPX contains the value of exp( i*kx*x ) for a given kx and ion before a move.
! EXPX_NEW contains the value of exp( i*kx*x ) for a given kx and ion after a move.

complex, dimension(0:Kmax, Ns), intent(in)				:: EXPX
complex, dimension(-Kmax:Kmax, Ns), intent(in)			:: EXPY, EXPZ
complex, dimension(0:Kmax, Ns), intent(out)				:: EXPX_NEW
complex, dimension(-Kmax:Kmax, Ns), intent(out)			:: EXPY_NEW, EXPZ_NEW

! SUMQEXPV contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian before a move.
! SUMQEXPV_NEW contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian after a move.

complex, dimension(Nkvec), intent(in)				 	:: SUMQEXPV
complex, dimension(Nkvec), intent(out)					:: SUMQEXPV_NEW

! DU_FOURIER is the change in fourier energy for a given move.

real, intent(out)											:: DU_FOURIER

! Local variables

integer													:: i, k
integer													:: kkx, kky, kkz
integer													:: typess
real														:: xi, yi, zi
real														:: bx, by, bz
real, parameter											:: Pi = 3.14159265359
complex													:: expv_diff



DU_FOURIER = 0.0

EXPX_NEW = ( 1.0, 0.0 )
EXPY_NEW = ( 1.0, 0.0 )
EXPZ_NEW = ( 1.0, 0.0 )

do i = 1, Ns

	xi = Xn(i)
	yi = Yn(i)
	zi = Zn(i)
	
	do k =1, Kmax

		bx = real(k) * 2.0 * Pi / BoxSize * xi
		by = real(k) * 2.0 * Pi / BoxSize * yi
		bz = real(k) * 2.0 * Pi / BoxSize * zi
		
		EXPX_NEW( k, i ) = cmplx( cos( bx ), sin( bx ) )
		EXPY_NEW( k, i ) = cmplx( cos( by ), sin( by ) )
		EXPZ_NEW( k, i ) = cmplx( cos( bz ), sin( bz ) )
		EXPY_NEW(-k, i ) = conjg( EXPY_NEW( k, i ) ) 
		EXPZ_NEW(-k, i ) = conjg( EXPZ_NEW( k, i ) )

	end do

end do

SUMQEXPV_NEW = SUMQEXPV

do i = 1, Ns

	typess = TYPEs(i)

	do k = 1, Nkvec

		kkx = KX(k)
		kky = KY(k)
		kkz = KZ(k)

		expv_diff = EXPX_NEW( kkx, i ) * EXPY_NEW( kky, i ) *  &
					EXPZ_NEW( kkz, i ) - &
					EXPX( kkx, i ) * EXPY( kky, i ) * &
					EXPZ( kkz, i )

		SUMQEXPV_NEW( k ) = SUMQEXPV_NEW( k ) + &
							   expv_diff * CHARGE( typess)

	end do

end do

do k = 1, Nkvec

	DU_FOURIER = DU_FOURIER + CONST(k) / BoxSize *   &
					( conjg( SUMQEXPV_NEW( k ) ) * SUMQEXPV_NEW( k ) - &
					  conjg( SUMQEXPV( k ) ) * SUMQEXPV( k ) )

end do

return

end subroutine Fourier_Move





