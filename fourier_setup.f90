
subroutine Fourier_Setup( Kmax,  CONST, KX, KY, KZ, Nkvec )
 implicit none

! This routine calculates the constant part of the Fourier energy of the ES.

! Kmax is the number of terms to be used in the complex space summations.

integer, intent(in)							:: Kmax
									
! Alpha is an Ewald sum parameter, Alpha = kappa * L, for kappa in A + T.

! ! real, intent(in)							:: Alpha

! Nkvec is the number of k vectors generated.

integer, intent(out)						:: Nkvec

! CONST contains the constant part of the 
! Fourier summation for a given vector n, CONST(n).

real, dimension(Nkvec), intent(out)			:: CONST

! KX, KY, KZ contain the vector identity of the Nkvec vectors.

integer, dimension(Nkvec), intent(out)		:: KX, KY, KZ

! Local Stuff

integer										:: i, j, k
real											:: b, Ksq
real, parameter								:: Pi = 3.14159265359

b = Pi * Pi / ( Alpha * Alpha )

Nkvec = 0

do i = 0, Kmax
	do j = -Kmax, Kmax
		do k = -Kmax, Kmax
			
			Ksq = i*i + j*j + k*k
			
			if( ( Ksq /= 0 ) .AND. ( Ksq < Kmax * Kmax ) ) then
			
				Nkvec = Nkvec + 1

				CONST(Nkvec) = 1.0 / ( Pi * real( i*i + j*j + k*k ) ) * &
							   exp( -b * real( i*i + j*j + k*k ) )

				if( i == 0 ) CONST(Nkvec) = 0.5 * CONST(Nkvec)

				KX(Nkvec) = i
				KY(Nkvec) = j
				KZ(Nkvec) = k

			end if
		end do
	end do
end do

end subroutine Fourier_Setup





