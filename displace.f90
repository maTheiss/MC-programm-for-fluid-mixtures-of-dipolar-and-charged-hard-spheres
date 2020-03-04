
subroutine Disp_Rot( MaxSp, Nmol, Nlj, Nion, Xlj, Ylj, Zlj, TYPElj, &
					 Xion, Yion, Zion, TYPEion, BoxSize, Rcut, DXYZ, DROT, &
					 LENGTHlj, LENGTHion, SPECIES, STARTlj, STARTion, &
					 PROB_TR, PROB_SP, NTemp, BETA, Nham, LNWSTATE, &
					 Nljgrs, EPS, SIG, CP, ALP, RMAX, MASSlj, Niongrs, CHARGE, &
					 MASSion, Alpha, Kmax, Nkvec, KX, KY, KZ, CONST, &
					 EXPX, EXPY, EXPZ, SUMQEXPV, SUMQX, SUMQY, SUMQZ, ULJ, &
					 UFOURIER, UREAL, USURF, DispOrRot, SpeciesID, Success, Seed )

implicit none

! Nmol is the number of molecules in the simulation box.
! MaxSp is the maximum number of species in the system.
! Nlj is the number of LJ beads in the simulation box.
! Nion is the number of ionic beads in the simulation box.

integer, intent(in)									:: MaxSp
integer, dimension(0:MaxSp), intent(in)				:: Nmol
integer, intent(in)									:: Nlj
integer, intent(in)									:: Nion

! Xlj, Ylj, and Zlj are the coordinates of the LJ beads.
! Xion, Yion, and Zion are the coordinates of the ionic beads.

real, dimension(Nlj), intent(inout)					:: Xlj, Ylj, Zlj
real, dimension(Nion), intent(inout)				:: Xion, Yion, Zion

! TYPElj contains the group identity of each LJ bead.
! TYPEion contains the group identity of each ionic bead.

integer, dimension(Nlj), intent(in)					:: TYPElj
integer, dimension(Nion), intent(in)				:: TYPEion

! BoxSize is the length of the simulation box.

real, intent(in)									:: BoxSize
real, intent(in)									:: Rcut

! DXYZ is the maximum displacement for each species.
! DROT is the maximum rotation for each species.

real, dimension(MaxSp), intent(in)					:: DXYZ, DROT

! LENGTHlj contains the number of LJ beads in each molecule.
! LENGTHion	contains the number of ionic beads in each molecule.
! SPECIES contains the species identity of each molecule.
! STARTlj contains the starting LJ bead number for each molecule.
! STARTion contains the starting ionic bead number for each molecule.

integer, dimension(Nmol(0)), intent(in)				:: LENGTHlj, LENGTHion
integer, dimension(Nmol(0)), intent(in)				:: SPECIES
integer, dimension(Nmol(0)), intent(in)				:: STARTlj, STARTion

! PROB_TR contains the accumulative probability of translation or rotation.
! PROB_SP contains the accumulative probability of selecting a given species.

real, dimension(2), intent(in)						:: PROB_TR
real, dimension(MaxSp), intent(in)					:: PROB_SP

! NTemp is the number of temperatures being used.

integer, intent(in)									:: NTemp

! BETA contains the reciprical temperature.

real, dimension(NTemp,1), intent(in)				:: BETA

! Nham is the number of hamiltonians being used.

integer, intent(in)									:: Nham

! LNWSTATE contains the weight of each temperature and hamiltonian.

real, dimension(NTemp, Nham), intent(inout)			:: LNWSTATE

! Nljgrs is the number of LJ groups in the system.
! EPS is a rank 3 array containing the eps_ij parameters for each hamiltonian.
! SIG is a rank 3 array containing the sigma_ij parameters for each hamiltonian.
! CP contains the C_ij parameters for each hamiltonian.
! ALP contains the alpha_ij parameters for each hamiltonian.
! RMAX contains the Rmax_ij parameters for each hamiltonian.
									
integer, intent(in)									:: Nljgrs
real, dimension(Nljgrs, Nljgrs, Nham), intent(in)	:: EPS, SIG, CP, ALP, RMAX

! MASSlj contains the mass of the LJ groups.

real, dimension(Nljgrs), intent(in)					:: MASSlj

! Niongrs is the number of ionic groups in the system.
! CHARGE is a rank 2 array containing the charge of group i for each hamiltonian.
									
integer, intent(in)									:: Niongrs
real, dimension(Niongrs, Nham), intent(in)			:: CHARGE

! MASSion contains the mass of the ionic groups.

real, dimension(Niongrs), intent(in)				:: MASSion

! Alpha is an Ewald sum parameter, Alpha = kappa * L, for kappa in A + T.

real, intent(in)									:: Alpha

! Kmax is an Ewald sum parameter.
! Nkvec is the number of k-vectors used in the Fourier sum.
! KX, KY, KZ contain the vector identity of the Nkvec vectors.
! CONST contains the constant part of the Fourier summation for a given Nkvec.

integer, intent(in)									:: Kmax
integer, intent(in)									:: Nkvec
integer, dimension(Nkvec), intent(in)				:: KX, KY, KZ
real, dimension(Nkvec), intent(in)					:: CONST

! EXPX contains the value of exp( i*kx*x ) for a given kx and ion.

complex, dimension(0:Kmax, Nion), intent(inout)		:: EXPX
complex, dimension(-Kmax:Kmax, Nion), intent(inout)	:: EXPY, EXPZ

! SUMQEXPV contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian.

complex, dimension(Nkvec, Nham), intent(inout)	 	:: SUMQEXPV

! SUMQX is the summation of qi * xi for all ions in the box.

real, dimension(Nham), intent(inout)				:: SUMQX, SUMQY, SUMQZ

! ULJ is the LJ energy of the system without the long range correction.
! UFOURIER is the coulombic fourier energy of the system.
! UREAL is the coulombic real energy of the system.
! USURF is the coulombic surface energy of the system.

real, dimension(Nham), intent(inout)				:: ULJ, UFOURIER
real, dimension(Nham), intent(inout)				:: UREAL, USURF

! DispOrRot is a flag to indicate whether a displacement or rotation was attempted.
! DispOrRot = 1 for Displacement.
! DispOrRot = 2 for Rotation.

integer, intent(out)								:: DispOrRot

! SpeciesID gives the identity of the species that was displaced or rotated.

integer, intent(out)								:: SpeciesID

! Success is a logical indicating whether the move was successful or not.

logical, intent(out)								:: Success

! Seed is the current random number generator seed value.

integer, intent(inout)								:: Seed 

real, external										:: ran2

! Local Variables.

integer												:: Mol, MolSpecies
integer												:: i, Count
integer												:: axis
integer												:: lenlj, stlj, endlj
integer												:: lenion, stion, endion
integer, dimension( Nlj )							:: temp4
integer, dimension( Nion )							:: temp8
integer, dimension( Nlj + Nion )					:: temp12

real												:: r, CoulCombo
real												:: Largest, PiRatio
real												:: dx, dy, dz, dtheta
real												:: xcom, ycom, zcom
real, allocatable, dimension( : )					:: Xlj_old, Ylj_old, Zlj_old
real, allocatable, dimension( : )					:: Xlj_new, Ylj_new, Zlj_new
real, dimension( Nlj )								:: temp1, temp2, temp3
real, dimension( Nion )								:: temp5, temp6, temp7
real, dimension( Nlj + Nion )						:: temp9, temp10, temp11
real, dimension( Nljgrs + Niongrs )					:: temp13
real, allocatable, dimension( : )					:: Xion_new, Yion_new, Zion_new
real, allocatable, dimension( : )					:: Xion_old, Yion_old, Zion_old
real, allocatable, dimension( : )					:: DELTAX, DELTAY, DELTAZ
real, dimension(2,2)								:: M
real, dimension(2)									:: T
real, dimension(NTemp, Nham)						:: BETA_HAM
real, dimension(Nham)								:: ULJ_old, ULJ_new, dULJ
real, dimension(Nham)								:: UREAL_old, UREAL_new, dUREAL
real, dimension(Nham)								:: SUMQX_NEW, SUMQY_NEW, SUMQZ_NEW
real, dimension(Nham)								:: dUSURF
real, dimension(Nham)								:: dUFOURIER
real, dimension(Nham,1)								:: dU

real, parameter										:: Pi = 3.14159265359
real, parameter										:: ec = 1.60217733e-19
real, parameter										:: eps0 = 8.854187817e-12
real, parameter										:: kB = 1.380658e-23

complex, allocatable, dimension( : , : )			:: EXPX_NEW
complex, allocatable, dimension( : , : )			:: EXPY_NEW, EXPZ_NEW
complex, dimension(Nkvec, Nham)					 	:: SUMQEXPV_NEW


Success = .False.

r = ran2(Seed)
SpeciesID = 0
i = 1

do while ( SpeciesID == 0 )

	if( r < PROB_SP(i) ) SpeciesID = i

	i = i + 1

end do


MolSpecies = int( Nmol( SpeciesID ) * ran2(Seed) ) + 1

i = 0
Count = 0

do while ( Count < MolSpecies )
	
	i = i + 1
	
	if( SPECIES(i) == SpeciesID ) Count = Count + 1

end do


Mol = i

lenlj = LENGTHlj( Mol )
stlj = STARTlj( Mol )
endlj = stlj + lenlj - 1

lenion = LENGTHion( Mol )
stion = STARTion( Mol )
endion = stion + lenion - 1

allocate( Xlj_old(lenlj) )
allocate( Ylj_old(lenlj) )
allocate( Zlj_old(lenlj) )

allocate( Xlj_new(lenlj) )
allocate( Ylj_new(lenlj) )
allocate( Zlj_new(lenlj) )

Xlj_old( 1:lenlj ) = Xlj( stlj:endlj )
Ylj_old( 1:lenlj ) = Ylj( stlj:endlj )
Zlj_old( 1:lenlj ) = Zlj( stlj:endlj )

if( lenion > 0 ) then

	allocate( Xion_old(lenion) )
	allocate( Yion_old(lenion) )
	allocate( Zion_old(lenion) )

	allocate( Xion_new(lenion) )
	allocate( Yion_new(lenion) )
	allocate( Zion_new(lenion) )

	allocate( DELTAX(lenion) )
	allocate( DELTAY(lenion) )
	allocate( DELTAZ(lenion) )

	allocate( EXPX_NEW(0:Kmax, lenion ) )
	allocate( EXPY_NEW(-Kmax:Kmax, lenion ) )
	allocate( EXPZ_NEW(-Kmax:Kmax, lenion ) )
	
	Xion_old( 1:lenion ) = Xion( stion:endion )
	Yion_old( 1:lenion ) = Yion( stion:endion )
	Zion_old( 1:lenion ) = Zion( stion:endion )

end if

if( stlj == 1 )	then

	if( Nmol(0) == 1 ) then
	
		ULJ_old	= 0.0

	else

		call e6molecule_cutoff( lenlj, Xlj_old, Ylj_old, Zlj_old, TYPElj(1:lenlj), &
						 Nlj - lenlj, Xlj(endlj+1:Nlj), Ylj(endlj+1:Nlj), &
						 Zlj(endlj+1:Nlj), TYPElj(endlj+1:Nlj), Nham, &
						 Nljgrs, EPS, SIG, CP, ALP, RMAX, BoxSize, Rcut, ULJ_old )

	end if

else if( stlj + lenlj - 1 == Nlj ) then

	call e6molecule_cutoff( lenlj, Xlj_old, Ylj_old, Zlj_old, TYPElj(stlj:endlj), &
					 Nlj - lenlj, Xlj(1:stlj-1), Ylj(1:stlj-1), Zlj(1:stlj-1), &
					 TYPElj(1:stlj-1), Nham, Nljgrs, EPS, SIG, CP, ALP, &
					 RMAX, BoxSize, Rcut, ULJ_old )

else

	temp1( 1:stlj-1 ) = Xlj( 1:stlj-1 )
	temp1( stlj:Nlj-lenlj ) = Xlj( endlj+1:Nlj )
	
	temp2( 1:stlj-1 ) = Ylj( 1:stlj-1)
	temp2( stlj:Nlj-lenlj ) = Ylj( endlj+1:Nlj)
	
	temp3( 1:stlj-1 ) = Zlj( 1:stlj-1)
	temp3( stlj:Nlj-lenlj ) = Zlj( endlj+1:Nlj)
	
	temp4( 1:stlj-1 ) = TYPElj( 1:stlj-1)
	temp4( stlj:Nlj-lenlj ) = TYPElj( endlj+1:Nlj)

	call e6molecule_cutoff( lenlj, Xlj_old, Ylj_old, Zlj_old, TYPElj(stlj:endlj), &
					 Nlj - lenlj, temp1, temp2, temp3, temp4, &
					 Nham, Nljgrs, EPS, SIG, CP, ALP, RMAX, &
					 BoxSize, Rcut, ULJ_old )

end if

if( ran2(Seed) < PROB_TR(1) ) then
	
	DispOrRot = 1

	dx = ( 2.0 * ran2(Seed) - 1.0 ) * DXYZ( SpeciesID ) * BoxSize
	dy = ( 2.0 * ran2(Seed) - 1.0 ) * DXYZ( SpeciesID ) * BoxSize
	dz = ( 2.0 * ran2(Seed) - 1.0 ) * DXYZ( SpeciesID ) * BoxSize

	Xlj_new = Xlj_old + dx
	Ylj_new = Ylj_old + dy
	Zlj_new = Zlj_old + dz

	if( lenion > 0 ) then

		Xion_new = Xion_old + dx
		Yion_new = Yion_old + dy
		Zion_new = Zion_old + dz
	
	end if

else 

	DispOrRot = 2

	call Outfold( lenlj, lenion, BoxSize, Xlj_old, Ylj_old, Zlj_old, &
				  Xion_old, Yion_old, Zion_old )

	if( lenion == 0 ) then

		call CenterOfMass( lenlj, Xlj_old, Ylj_old, Zlj_old,  &
						   TYPElj( stlj:endlj ), Nljgrs, MASSlj, & 
						   xcom, ycom, zcom )

	else

		temp9( 1:lenlj ) = Xlj_old( 1:lenlj )
		temp9( lenlj+1:lenlj+lenion ) = Xion_old( 1:lenion )
		
		temp10( 1:lenlj ) = Ylj_old( 1:lenlj )
		temp10( lenlj+1:lenlj+lenion ) = Yion_old( 1:lenion )
		
		temp11( 1:lenlj ) = Zlj_old( 1:lenlj )
		temp11( lenlj+1:lenlj+lenion ) = Zion_old( 1:lenion )
		
		temp12( 1:lenlj ) = TYPElj( stlj:endlj )	 
		temp12( lenlj+1:lenlj+lenion ) = TYPEion( 1:lenion ) + Nljgrs
		
		temp13( 1:Nljgrs ) = MASSlj( 1:Nljgrs )			
		temp13( Nljgrs+1:Nljgrs+Niongrs ) = MASSion( 1:Niongrs )

		call CenterOfMass( lenlj + lenion, temp9, temp10, temp11, temp12, &
						   Nljgrs + Niongrs, temp13, xcom, ycom, zcom )
	end if

	dtheta = ( 2.0 * ran2(Seed) - 1.0 ) * Pi * DROT( SpeciesID )

	M = reshape( (/ cos( dtheta ), -sin( dtheta ), sin( dtheta ), cos( dtheta ) /), &
			     (/ 2, 2 /) )

	axis = int( 3.0 * ran2(Seed) ) + 1

	
	Select Case ( axis )
		
		case ( 1 )			! Rotation in y-z plane.

			do i = 1, lenlj
				
				T(1) = Ylj_old(i) - ycom
				T(2) = Zlj_old(i) - zcom

				T = matmul( M, T )

				Ylj_new(i) = T(1) + ycom
				Zlj_new(i) = T(2) + zcom
				Xlj_new(i) = Xlj_old(i)

			end do

			if( lenion > 0 ) then

				do i = 1, lenion
					
					T(1) = Yion_old(i) - ycom
					T(2) = Zion_old(i) - zcom

					T = matmul( M, T )

					Yion_new(i) = T(1) + ycom
					Zion_new(i) = T(2) + zcom
					Xion_new(i) = Xion_old(i)
			
				end do
		
			end if

		case ( 2 )			! Rotation in x-z plane.

			do i = 1, lenlj
				
				T(1) = Zlj_old(i) - zcom
				T(2) = Xlj_old(i) - xcom

				T = matmul( M, T )

				Zlj_new(i) = T(1) + zcom
				Xlj_new(i) = T(2) + xcom
				Ylj_new(i) = Ylj_old(i)

			end do

			if( lenion > 0 ) then

				do i = 1, lenion

					T(1) = Zion_old(i) - zcom
					T(2) = Xion_old(i) - xcom

					T = matmul( M, T )

					Zion_new(i) = T(1) + zcom
					Xion_new(i) = T(2) + xcom
					Yion_new(i) = Yion_old(i)
			
				end do
		
			end if

		case ( 3 )			! Rotation in x-y plane.

			do i = 1, lenlj

				T(1) = Xlj_old(i) - xcom
				T(2) = Ylj_old(i) - ycom

				T = matmul( M, T )

				Xlj_new(i) = T(1) + xcom
				Ylj_new(i) = T(2) + ycom
				Zlj_new(i) = Zlj_old(i)

			end do

			if( lenion > 0 ) then

				do i = 1, lenion

					T(1) = Xion_old(i) - xcom
					T(2) = Yion_old(i) - ycom

					T = matmul( M, T )

					Xion_new(i) = T(1) + xcom
					Yion_new(i) = T(2) + ycom
					Zion_new(i) = Zion_old(i)
			
				end do
		
			end if
	
	end select

end if


do i = 1, lenlj
	
	if( Xlj_new(i) > BoxSize )  Xlj_new(i) = Xlj_new(i) - &
								BoxSize * aint( Xlj_new(i) / BoxSize )
	if( Ylj_new(i) > BoxSize )  Ylj_new(i) = Ylj_new(i) - &
								BoxSize * aint( Ylj_new(i) / BoxSize )
	if( Zlj_new(i) > BoxSize )  Zlj_new(i) = Zlj_new(i) - &
								BoxSize * aint( Zlj_new(i) / BoxSize )

	if( Xlj_new(i) < 0.0 )  Xlj_new(i) = Xlj_new(i) - &
							BoxSize * aint( Xlj_new(i) / BoxSize - 1 )
	if( Ylj_new(i) < 0.0 )  Ylj_new(i) = Ylj_new(i) - &
							BoxSize * aint( Ylj_new(i) / BoxSize - 1 )
	if( Zlj_new(i) < 0.0 )  Zlj_new(i) = Zlj_new(i) - &
							BoxSize * aint( Zlj_new(i) / BoxSize - 1 )

end do


if( stlj == 1 )	then

	if( Nmol(0) == 1 ) then

		ULJ_new = 0.0

	else

		call e6molecule_cutoff( lenlj, Xlj_new, Ylj_new, Zlj_new, TYPElj(1:lenlj), &
						 Nlj - lenlj, Xlj(endlj+1 : Nlj), Ylj(endlj+1 : Nlj), &
						 Zlj(endlj+1 : Nlj), TYPElj(endlj+1	: Nlj), Nham, &
						 Nljgrs, EPS, SIG, CP, ALP, RMAX, BoxSize, Rcut, ULJ_new	)

	end if

else if( stlj + lenlj - 1 == Nlj ) then

	call e6molecule_cutoff( lenlj, Xlj_new, Ylj_new, Zlj_new, TYPElj(stlj:endlj), &
					 Nlj - lenlj, Xlj(1:stlj-1), Ylj(1:stlj-1), Zlj(1:stlj-1), &
					 TYPElj(1:stlj-1), Nham, Nljgrs, EPS, SIG, CP, ALP, RMAX, &
					 BoxSize, Rcut, ULJ_new )

else

	call e6molecule_cutoff( lenlj, Xlj_new, Ylj_new, Zlj_new, TYPElj(stlj:endlj), &
					 Nlj - lenlj, temp1, temp2, temp3, temp4, &
					 Nham, Nljgrs, EPS, SIG, CP, ALP, RMAX, BoxSize, Rcut, ULJ_new )

end if

dULJ = ULJ_new - ULJ_old

if( lenion > 0 ) then

CoulCombo = ec * ec * 1.0e10 / ( 4.0 * Pi * eps0 * kB )

	do i = 1, lenion
	
		if( Xion_new(i) > BoxSize )	Xion_new(i) = Xion_new(i) - &
									BoxSize * aint( Xion_new(i) / BoxSize )
		if( Yion_new(i) > BoxSize ) Yion_new(i) = Yion_new(i) - &
									BoxSize * aint( Yion_new(i) / BoxSize )
		if( Zion_new(i) > BoxSize ) Zion_new(i) = Zion_new(i) - &
									BoxSize * aint( Zion_new(i) / BoxSize )

		if( Xion_new(i) < 0.0 ) Xion_new(i) = Xion_new(i) - &
								BoxSize * aint( Xion_new(i) / BoxSize - 1 )
		if( Yion_new(i) < 0.0 ) Yion_new(i) = Yion_new(i) - &
								BoxSize * aint( Yion_new(i) / BoxSize - 1 )
		if( Zion_new(i) < 0.0 ) Zion_new(i) = Zion_new(i) - &
								BoxSize * aint( Zion_new(i) / BoxSize - 1 )
	
	end do

	
	do i = 1, lenion
	
		if( Xion_old(i) > BoxSize ) Xion_old(i) = Xion_old(i) - &
									BoxSize * aint( Xion_old(i) / BoxSize )
		if( Yion_old(i) > BoxSize ) Yion_old(i) = Yion_old(i) - &
									BoxSize * aint( Yion_old(i) / BoxSize )
		if( Zion_old(i) > BoxSize ) Zion_old(i) = Zion_old(i) - &
									BoxSize * aint( Zion_old(i) / BoxSize )

		if( Xion_old(i) < 0.0 ) Xion_old(i) = Xion_old(i) - &
								BoxSize * aint( Xion_old(i) / BoxSize - 1 )
		if( Yion_old(i) < 0.0 ) Yion_old(i) = Yion_old(i) - &
								BoxSize * aint( Yion_old(i) / BoxSize - 1 )
		if( Zion_old(i) < 0.0 ) Zion_old(i) = Zion_old(i) - &
								BoxSize * aint( Zion_old(i) / BoxSize - 1 )
	
	end do

	if( stion == 1 ) then

		if( Nmol(0) == 1 ) then

			UREAL_old = 0.0

		else
	
			call RealMolecule( lenion, Xion_old, Yion_old, Zion_old, TYPEion(1:lenion), &
							   Nion - lenion, Xion(endion+1 : Nion), &
							   Yion(endion+1 : Nion), Zion(endion+1 : Nion), &
							   TYPEion(endion+1 : Nion), Nham, Niongrs, CHARGE, &
							   BoxSize, Alpha, UREAL_old )

		end if

	else if( stion + lenion - 1 == Nion ) then

		call RealMolecule( lenion, Xion_old, Yion_old, Zion_old, TYPEion(stion:endion), &
						   Nion - lenion, Xion(1:stion-1), Yion(1:stion-1), Zion(1:stion-1), &
						   TYPEion(1:stion-1), Nham, Niongrs, CHARGE, BoxSize, Alpha, UREAL_old )
	
	else

		temp5( 1:stion-1 ) = Xion( 1:stion-1 )
		temp5( stion:Nion-lenion ) = Xion( endion+1:Nion )
		
		temp6( 1:stion-1 ) = Yion( 1:stion-1 )
		temp6( stion:Nion-lenion ) = Yion( endion+1:Nion )
		
		temp7( 1:stion-1 ) = Zion( 1:stion-1 )
		temp7( stion:Nion-lenion ) = Zion( endion+1:Nion )
		
		temp8( 1:stion-1 ) = TYPEion( 1:stion-1 )
		temp8( stion:Nion-lenion ) = TYPEion( endion+1:Nion )

		call RealMolecule( lenion, Xion_old, Yion_old, Zion_old, TYPEion(stion:endion), &
						   Nion - lenion, temp5, temp6, temp7, temp8, &
						   Nham, Niongrs, CHARGE, BoxSize, Alpha, UREAL_old )

	end if
	
	if( stion == 1 ) then

		if( Nmol(0) == 1 ) then

			UREAL_new = 0.0

		else
	
			call RealMolecule( lenion, Xion_new, Yion_new, Zion_new, TYPEion(1:lenion), &
							   Nion - lenion, Xion(endion+1 : Nion), Yion(endion+1 : Nion), &
							   Zion(endion+1 : Nion), TYPEion(endion+1 : Nion), Nham, &
							   Niongrs, CHARGE, BoxSize, Alpha, UREAL_new )

		end if

	else if( stion + lenion - 1 == Nion ) then

		call RealMolecule( lenion, Xion_new, Yion_new, Zion_new, TYPEion(stion:endion), &
						   Nion - lenion, Xion(1:stion-1), Yion(1:stion-1), Zion(1:stion-1), &
						   TYPEion(1:stion-1), Nham, Niongrs, CHARGE, BoxSize, Alpha, UREAL_new )

	else
	
		call RealMolecule( lenion, Xion_new, Yion_new, Zion_new, TYPEion(stion:endion), &
						   Nion - lenion, temp5, temp6, temp7, temp8, &
						   Nham, Niongrs, CHARGE, BoxSize, Alpha, UREAL_new )

	end if

	dUREAL = UREAL_new - UREAL_old

	DELTAX = Xion_new - Xion_old
	DELTAY = Yion_new - Yion_old
	DELTAZ = Zion_new - Zion_old

	call Surf_Move( lenion, DELTAX, DELTAY, DELTAZ, &
					TYPEion( stion:endion ), Nham, Niongrs, CHARGE, BoxSize, &
					SUMQX, SUMQY, SUMQZ, SUMQX_NEW, SUMQY_NEW, SUMQZ_NEW, dUSURF )

	call Fourier_Move( lenion, Xion_new, Yion_new, Zion_new, TYPEion( stion:endion ), &
					   Nham, Niongrs, CHARGE, BoxSize, Kmax, Nkvec, KX, KY, KZ, &
					   CONST, EXPX(:,stion:endion), EXPY(:,stion:endion), &
					   EXPZ(:,stion:endion), EXPX_NEW(:,1:lenion), &
					   EXPY_NEW(:,1:lenion), EXPZ_NEW(:,1:lenion), &
					   SUMQEXPV, SUMQEXPV_NEW, dUFOURIER )

	dUREAL = dUREAL	* CoulCombo
	dUSURF = dUSURF	* CoulCombo
	dUFOURIER = dUFOURIER * CoulCombo

else

	dUREAL = 0.0
	dUSURF = 0.0
	dUFOURIER = 0.0

end if


dU(:,1) = - ( dULJ(:) + dUREAL(:) + dUSURF(:) + dUFOURIER(:) )


BETA_HAM = matmul( BETA, transpose( dU ) )

Largest = maxval( LNWSTATE + BETA_HAM )

PiRatio = log( sum( exp( LNWSTATE + BETA_HAM - Largest ) ) ) + Largest

if( log( ran2(Seed) ) < PiRatio ) then

	Success = .True.

	Xlj( stlj:endlj ) = Xlj_new( 1:lenlj )
	Ylj( stlj:endlj ) = Ylj_new( 1:lenlj )
	Zlj( stlj:endlj ) = Zlj_new( 1:lenlj )

	ULJ = ULJ + dULJ

	LNWSTATE = LNWSTATE + BETA_HAM - PiRatio

	if( lenion > 0 ) then

		Xion( stion:endion ) = Xion_new( 1:lenion )
		Yion( stion:endion ) = Yion_new( 1:lenion )
		Zion( stion:endion ) = Zion_new( 1:lenion )

		EXPX(:,stion:endion) = EXPX_NEW(:,1:lenion)
		EXPY(:,stion:endion) = EXPY_NEW(:,1:lenion)
		EXPZ(:,stion:endion) = EXPZ_NEW(:,1:lenion)

		SUMQEXPV = SUMQEXPV_NEW

		SUMQX = SUMQX_NEW
		SUMQY = SUMQY_NEW
		SUMQZ = SUMQZ_NEW

		UFOURIER = UFOURIER + dUFOURIER
		UREAL = UREAL + dUREAL
		USURF = USURF + dUSURF

	end if

end if

deallocate( Xlj_old )
deallocate( Ylj_old )
deallocate( Zlj_old )

deallocate( Xlj_new )
deallocate( Ylj_new )
deallocate( Zlj_new )

if( lenion > 0 ) then

	deallocate( Xion_old )
	deallocate( Yion_old )
	deallocate( Zion_old )

	deallocate( Xion_new )
	deallocate( Yion_new )
	deallocate( Zion_new )

	deallocate( DELTAX )
	deallocate( DELTAY )
	deallocate( DELTAZ )

	deallocate( EXPX_NEW )
	deallocate( EXPY_NEW )
	deallocate( EXPZ_NEW )

end if

return

end subroutine Disp_Rot








