
subroutine CreateLatticeKRZ( Nljgrs, Niongrs, Nsp, &
				   MaxSp, MaxBeads,&
				   SIG, CHARGE,  &
				   LENLJ, LENION, BEADTYPE, GROUPTYPE, &
				   BoxSize, &
				   Nlj, Nion, MaxMol, MaxLJ, MaxIon, Nmol, &
				   Xlj, Ylj, Zlj, TYPElj, Xion, Yion, Zion, TYPEion, Xint, Yint, Zint,&
				   SPECIES, LENGTHlj, LENGTHion, STARTlj, STARTion,	&
				   Alpha, Kmax, Nkvec, KX, KY, KZ, CONST, &
				   EXPX, EXPY, EXPZ, SUMQEXPV, SUMQX, SUMQY, SUMQZ, &
				   UREAL, UFOURIER, USURF, USELF_CH, &
				   USELF_MOL, Seed )
				  
implicit none

! Nljgrs is the number of Lennard-Jones groups in the simulation.
! Niongrs is the number of ionic groups in the simulation.
! Nsp is the number of species in the simulation.

integer, intent(in)										:: Nljgrs
integer, intent(in)										:: Niongrs
integer, intent(in)										:: Nsp

! MaxSp is the maximum number of species in the simulation.
! MaxBeads is the maximum number of LJ and ionic beads in a molecule.
! MaxInt is the maximum number of integer parameters for a CB method.
! MaxReal is the maximum number of real parameters for a CB method.

integer, intent(in)										:: MaxSp
integer, intent(in)										:: MaxBeads


! SIG contains the sigma_ij parameters for each hamiltonian.

real, dimension(Nljgrs,Nljgrs), intent(in)				:: SIG


! CHARGE contains the charge of a bead for each hamiltonian.

real, dimension(Niongrs), intent(in)						:: CHARGE


! LENLJ is the LJ length of each species.
! LENION is the ionic length of each species.
! BEADTYPE indicates whether a bead is LJ or ionic.
! GROUPTYPE indicates the group identity of each bead.

integer, dimension(MaxSp), intent(in)					:: LENLJ 
integer, dimension(MaxSp), intent(in)					:: LENION
character*5, dimension(MaxBeads,MaxSp), intent(in)		:: BEADTYPE
integer, dimension(MaxBeads,MaxSp), intent(in)			:: GROUPTYPE

real, dimension(MaxBeads,MaxSp), intent(in)				:: Xint, Yint, Zint


! BoxSize contains the length of the simulation boxes.

real, intent(in)											:: BoxSize

! Nlj is the number of LJ beads in the simulation boxes.
! Nion is the number of ionic beads in the simulation boxes.

integer, intent(out)										:: Nlj
integer, intent(out)										:: Nion

! MaxMol is the maximum number of molecules per box.
! MaxLJ is the maximum number of LJ beads per box.
! MaxIon is the maximum number of ionic beads per box.

integer, intent(in)										:: MaxMol
integer, intent(in)										:: MaxLJ
integer, intent(in)										:: MaxIon

! Nmol is the number of molecules in the simulation boxes.
! Nmol(0) ==> total number of molecules.
! Nmol(1:MaxSp) ==> number of molecules of that species.

integer, dimension(0:MaxSp), intent(in)					:: Nmol

! Xlj, Ylj, and Zlj are the coordinates of the LJ beads.
! Xion, Yion, and Zion are the coordinates of the ionic beads.

real, dimension( MaxLJ ), intent(out)						:: Xlj, Ylj, Zlj
real, dimension( MaxIon ), intent(out)					:: Xion, Yion, Zion

! TYPElj contains the group identity of each LJ bead.
! TYPEion contains the group identity of each ionic bead.

integer, dimension( MaxLJ ), intent(out)					:: TYPElj
integer, dimension( MaxIon ), intent(out)				:: TYPEion

! SPECIES contains the species identity of each molecule.
! LENGTHlj contains the number of LJ beads in each molecule.
! LENGTHion	contains the number of ionic beads in each molecule.
! STARTlj contains the starting LJ bead number for each molecule.
! STARTion contains the starting ionic bead number for each molecule.

integer, dimension(MaxMol), intent(out)				:: SPECIES
integer, dimension(MaxMol), intent(out)				:: LENGTHlj, LENGTHion
integer, dimension(MaxMol), intent(out)				:: STARTlj, STARTion

! Alpha is an Ewald sum parameter, Alpha = kappa * L, for kappa in A + T.

real, intent(in)										:: Alpha

! Kmax is an Ewald sum parameter.
! Nkvec is the number of k-vectors used in the Fourier sum.
! KX, KY, KZ contain the vector identity of the Nkvec vectors.
! CONST contains the constant part of the Fourier summation for a given Nkvec.

integer, intent(in)										:: Kmax
integer, intent(in)										:: Nkvec
integer, dimension(Nkvec), intent(in)					:: KX, KY, KZ
real, dimension(Nkvec), intent(in)						:: CONST

! EXPX contains the value of exp( i*kx*x ) for a given kx and ion.

complex, dimension(0:Kmax, MaxIon), intent(out) 			:: EXPX
complex, dimension(-Kmax:Kmax, MaxIon), intent(out) 		:: EXPY, EXPZ

! SUMQEXPV contains the summation of qi*exp(i*(kx*x + ky*y + kz*z)) 
! for a given k-vector and hamiltonian.

complex, dimension(Nkvec), intent(inout)				 :: SUMQEXPV

! SUMQX is the summation of qi * xi for all ions in the box.

real, intent(inout)									:: SUMQX, SUMQY, SUMQZ

! UFOURIER is the coulombic fourier energy of the system.
! UREAL is the coulombic real energy of the system.
! USURF is the coulombic surface energy of the system.
! USELF_CH is the summation of the square of all the charges.
! USELF_MOL is the self energy of a given molecule.

real, intent(inout)									:: UREAL, UFOURIER
real, intent(inout)									:: USURF, USELF_CH
real, dimension(MaxMol), intent(inout)				:: USELF_MOL

! Seed is the current random number generator seed value.

integer, intent(inout)									:: Seed

real, external											:: ran2

! Local stuff

integer											:: i, j, k, counter, counter2
integer											:: x, y, z
integer											:: NmolPerAxis
integer											:: LatticePlace
integer											:: NFreePlaces
integer											:: MolNo
integer, dimension(0:Nsp)							:: Nmol_new
integer											:: sp
integer											:: SpeciesID
integer											:: stion, endion
integer											:: LatNumb
integer											:: nBeads
integer											:: NLayers
integer											:: MaxMolecs

logical											:: Ions = .false.

real												:: CoulCombo
real												:: MolDist
real												:: MinDist
real												:: MinLayerDist, LayerDist
real												:: tmp
real												:: Random
real, dimension(Nsp)								:: PROB_SP
real												:: dUFOURIER
real												:: theta, phi
real, allocatable, dimension(:)					:: Xtmp, Ytmp, Ztmp


logical, allocatable, dimension(:)				:: Occupied

real, allocatable, dimension(:)					:: Xlat, Ylat, Zlat
real, allocatable, dimension(:)					:: Xlayer, Ylayer

real, parameter									:: Pi = 3.14159265359
real, parameter									:: ec = 1.60217733e-19
real, parameter									:: eps0 = 8.854187817e-12
real, parameter									:: kB = 1.380658e-23

complex, dimension(0:Kmax, MaxIon)				:: EXPX_new
complex, dimension(-Kmax:Kmax, MaxIon)			:: EXPY_new, EXPZ_new
complex, dimension(Nkvec)							:: SUMQEXPV_new

CoulCombo = ec * ec * 1.0e10 / ( 4.0 * Pi * eps0 * kB )

do i = 1, Nsp

	if( LENION(i) > 0 ) Ions = .true.

end do

! Minimaler Teilchenabstand
MinDist = maxval(SIG + 0.00001)

! Anzahl der Teilchen, die in einer Achsrichtung Platz jaben
NmolPerAxis = nint(BoxSize / MinDist - 0.5)

! Tatsächlicher Teilchenabstand
MolDist = BoxSize / real(NmolPerAxis)

! Minimaler Layerabstand
MinLayerDist = sqrt(MinDist ** 2.0 - 0.5 * MolDist ** 2.0)

! Anzahl an Layern, die Platz haben
NLayers = nint(BoxSize / MinLayerDist - 0.5)

! Tatsächlicher Layerabstand
LayerDist = BoxSize / real(NLayers)

! Maximale Anzahl an Gitterplätzen
MaxMolecs = NmolPerAxis ** 2 * NLayers
if (MaxMolecs .LT. Nmol(0)) then

	write(*,*) 'Zu wenig Gitterplätze'
	write(*,*) 'Plätze:', MaxMolecs
	write(*,*) 'Moleküle:', Nmol(0)
	stop


end if


allocate(Xlat(MaxMolecs))
allocate(Ylat(MaxMolecs))
allocate(Zlat(MaxMolecs))

allocate(Xlayer(NmolPerAxis**2))
allocate(Ylayer(NmolPerAxis**2))

allocate(Occupied(MaxMolecs))


! Create lattice coordinates

! Create X and Y coordinates of one layer

counter = 0
do x=1, NmolPerAxis

	do y=1, NmolPerAxis

		counter = counter + 1

		Xlayer(counter) = (x-1) * MolDist
		Ylayer(counter) = (y-1) * MolDist

	end do

end do

! Stack layers in z coordinate
counter = 0

do z=1, NLayers

	counter2 = 0

	do i=1, NmolPerAxis

		do j=1, NmolPerAxis

			counter2 = counter2 + 1
			counter = counter + 1

			Zlat(counter) = (z-1) * LayerDist

			! Wenn z ungerade, (Layer 0, 2, 4, 6), dann nehm Koordinaten X und Y normal
			if( mod(z,2) .NE. 0) then
				Xlat(counter) = Xlayer(counter2)
				Ylat(counter) = Ylayer(counter2)
			else
		! Wenn z gerade, dann verschieb Layer um MolDist/2 in x- und in y-Richtung
				Xlat(counter) = Xlayer(counter2) + 0.5 * MolDist
				Ylat(counter) = Ylayer(counter2) + 0.5 * MolDist
			end if

		end do

	end do

end do

counter = 0

!do z=1, NmolPerAxis
!
!	do y=1, NmolPerAxis
!
!		do x=1, NmolPerAxis
!
!			counter = counter + 1
!
!			Xlat(counter) = (x-1) * MolDist
!			Ylat(counter) = (y-1) * MolDist
!			Zlat(counter) = (z-1) * MolDist
!
!		end do
!
!	end do
!
!end do

PROB_SP(1:Nsp) = real( Nmol(1:Nsp) )

do j = 2, Nsp

	PROB_SP(j) = PROB_SP(j-1) + PROB_SP(j)

end do

PROB_SP = PROB_SP / PROB_SP(Nsp)


Nmol_new = 0
Nlj = 0
Nion = 0
MolNo = 1
NFreePlaces = MaxMolecs
Occupied = .false.

do while (Nmol_new(0) .LT. Nmol(0))

	! Random Species

	Random = ran2(Seed)

	sp = 0
	j = 1

	do while ( sp == 0 )

		if( Random < PROB_SP(j) ) sp = j

		j = j + 1

	end do

	if( Nmol_new(sp) == Nmol(sp) ) cycle

	! Random position on lattice

	LatticePlace = nint(ran2(Seed) * (NFreePlaces-1)) + 1

	! Find the corresponding free position on lattice

	k = 1
	counter = 0
	do while (counter .LT. LatticePlace)

		if (.NOT. Occupied(k)) counter = counter + 1

		k = k + 1

	end do
	LatNumb = k - 1

	LENGTHlj(MolNo) = 0
	LENGTHion(MolNo) = 0

	STARTlj(MolNo) = Nlj + 1
	STARTion(MolNo) = Nion + 1

	nBeads = LENLJ(sp) + LENION(sp)

	allocate(Xtmp(nBeads))
	allocate(Ytmp(nBeads))
	allocate(Ztmp(nBeads))

	Xtmp = Xint(:,sp)
	Ytmp = Yint(:,sp)
	Ztmp = Zint(:,sp)

	do j= 1, 3

		theta = ran2(Seed) * 2 * Pi
    	phi = ran2(Seed) * 2 * Pi

    	call Rotate(nBeads, Xtmp, Ytmp, Ztmp, theta, phi)

	end do

	do j = 1,  nBeads

		if (BEADTYPE(j,sp) .EQ. 'LJ') then


			Nlj = Nlj + 1
			LENGTHlj(MolNo) = LENGTHlj(MolNo) + 1

			Xlj(Nlj) = Xlat(LatNumb) + Xtmp(j)
			Ylj(Nlj) = Ylat(LatNumb) + Ytmp(j)
			Zlj(Nlj) = Zlat(LatNumb) + Ztmp(j)

			TYPElj(Nlj) = GROUPTYPE(j,sp)

		else if (BEADTYPE(j,sp) .EQ. 'Ion') then

			Nion = Nion + 1
			LENGTHion(MolNo) = LENGTHion(MolNo) + 1

			Xion(Nion) = Xlat(LatNumb) + Xtmp(j)
			Yion(Nion) = Ylat(LatNumb) + Ytmp(j)
			Zion(Nion) = Zlat(LatNumb) + Ztmp(j)

			TYPEion(Nion) = GROUPTYPE(j,sp)

		end if

	end do

	deallocate(Xtmp)
	deallocate(Ytmp)
	deallocate(Ztmp)

	SPECIES(MolNo) = sp

	Nmol_new(0) = Nmol_new(0) + 1
	Nmol_new(sp) = Nmol_new(sp) + 1
	NFreePlaces = NFreePlaces - 1
	Occupied(LatNumb) = .true.
	MolNo = MolNo + 1

end do

do i = 1, Nlj

	if( Xlj(i) > BoxSize )  Xlj(i) = Xlj(i) - &
						BoxSize * aint( Xlj(i) / BoxSize )
	if( Ylj(i) > BoxSize )  Ylj(i) = Ylj(i) - &
						BoxSize * aint( Ylj(i) / BoxSize )
	if( Zlj(i) > BoxSize )  Zlj(i) = Zlj(i) - &
						BoxSize * aint( Zlj(i) / BoxSize )

	if( Xlj(i) < 0.0 )  Xlj(i) = Xlj(i) - &
						BoxSize * aint( Xlj(i) / BoxSize - 1 )
	if( Ylj(i) < 0.0 )  Ylj(i) = Ylj(i) - &
						BoxSize * aint( Ylj(i) / BoxSize - 1 )
	if( Zlj(i) < 0.0 )  Zlj(i) = Zlj(i) - &
						BoxSize * aint( Zlj(i) / BoxSize - 1 )

end do

do i = 1, Nion

	if( Xion(i) > BoxSize )	Xion(i) = Xion(i) - &
						BoxSize * aint( Xion(i) / BoxSize )
	if( Yion(i) > BoxSize ) Yion(i) = Yion(i) - &
						BoxSize * aint( Yion(i) / BoxSize )
	if( Zion(i) > BoxSize ) Zion(i) = Zion(i) - &
						BoxSize * aint( Zion(i) / BoxSize )

	if( Xion(i) < 0.0 ) Xion(i) = Xion(i) - &
						BoxSize * aint( Xion(i) / BoxSize - 1 )
	if( Yion(i) < 0.0 ) Yion(i) = Yion(i) - &
						BoxSize * aint( Yion(i) / BoxSize - 1 )
	if( Zion(i) < 0.0 ) Zion(i) = Zion(i) - &
						BoxSize * aint( Zion(i) / BoxSize - 1 )

end do

call HSOverlap_all( Nlj, Xlj, Ylj, Zlj, TYPElj, Nmol(0), LENGTHlj, &
					   Nljgrs, SIG, &
					   BoxSize)



do i=1, Nmol(0)

	SpeciesID = SPECIES(i)

	stion = STARTion(i)
	endion = STARTion(i)+LENGTHion(i)-1

	if ( lenion(SpeciesID) .GT. 1 ) then

		call SelfMolecule(LENGTHion(i), Xion(stion:endion), Yion(stion:endion), Zion(stion:endion), &
 						TYPEion(stion:endion), Niongrs, CHARGE, &
					 	BoxSize, Alpha, USELF_MOL(i))

		USELF_MOL(i) = USELF_MOL(i) * CoulCombo

    else

    	USELF_MOL(i) = 0.0

    end if

end do


UFOURIER = 0.0
EXPX = ( 0.0, 0.0 )
EXPY = ( 0.0, 0.0 )
EXPZ = ( 0.0, 0.0 )

SUMQEXPV = ( 0.0 )

do i=1, Nmol(0)

	SpeciesID = SPECIES(i)
	stion = STARTion(i)
	endion = STARTion(i) + LENGTHion(i)-1

	call Fourier_Move( LENGTHion(i), Xion(stion:endion), Yion(stion:endion), &
 					Zion(stion:endion), &
 					TYPEion(stion:endion), Niongrs, CHARGE,  &
						 BoxSize, Kmax, Nkvec, KX, KY, KZ, CONST, EXPX(:,stion:endion), &
						 EXPY(:,stion:endion), EXPZ(:,stion:endion), EXPX_NEW(:,stion:endion), &
						 EXPY_NEW(:,stion:endion), EXPZ_NEW(:,stion:endion), &
						 SUMQEXPV, SUMQEXPV_NEW, dUFOURIER )

		UFOURIER = UFOURIER + dUFOURIER
		EXPX(:,stion:endion) = EXPX_NEW(:,stion:endion)
		EXPY(:,stion:endion) = EXPY_NEW(:,stion:endion)
		EXPZ(:,stion:endion) = EXPZ_NEW(:,stion:endion)
		SUMQEXPV = SUMQEXPV_NEW

end do
UFOURIER = UFOURIER * CoulCombo

call RealInteract( Nion, Xion, Yion, Zion, TYPEion, Nmol(0), LENGTHion, &
					     Niongrs, CHARGE, BoxSize, Alpha, UREAL )
UREAL = UREAL * CoulCombo


USELF_CH = 0.0

do i=1, Nion

	USELF_CH = USELF_CH + CHARGE(TYPEion(i)) ** 2

end do
USELF_CH = alpha / sqrt(pi) / BoxSize * USELF_CH
USELF_CH = USELF_CH * CoulCombo


deallocate(Xlat)
deallocate(Ylat)
deallocate(Zlat)

deallocate(Xlayer)
deallocate(Ylayer)

deallocate(Occupied)


return

end subroutine CreateLatticeKRZ






