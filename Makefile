#gcc-compiler 
#-----------------------------------------------------------------------------------------------------------------
#debugging 
#FCOMPFLAGS    = -O2 -fimplicit-none -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace
#FCOMPFLAGS    = -O2 -Wconversion
#running
#FCOMPFLAGS    = -O3 -march=native -fimplicit-none -Wall -Wline-truncation -fwhole-file -std=f2008 
#-pg 
#-fopenmp 
 
#-fopenmp   ..parallel
				#gmon.out-file: Befehl: gprof main gmon.out > analysis.txt

#-O3  -traceback  -check all       
#FCOMPFLAGS    =	-O3  -traceback  -check all -profile-functions -profile-loops=all -profile-loops-report=2
#-fp-stack-check 
#FCOMPFLAGS    =	 -fast -parallel -par-report3
#FCOMPFLAGS    =	 -assume realloc_lhs -profile-functions -profile-loops=all -profile-loops-report=2 
#FCOMPFLAGS    =   -g -traceback -check all 
FCOMPFLAGS    = -Ofast

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

#LD            = ifort	
#FC            = ifort
LD            = gfortran		
FC            =	gfortran	

MAKEFILE      =	Makefile
PROGRAM       =	main

SRC	      =	working_prec.f90 \
		basic_parameters.f90 \
		sim_parameters.f90 \
		read_input.f90 \
		read_results_file.f90 \
		ran_gen.f90 \
		globals.f90 \
		utilities.f90 \
		adjust_displacement.f90 \
		read_positions.f90 \
		rw_files.f90 \
		lattice_fcc.f90 \
		ES_Fourier.f90 \
		realmolecule.f90 \
		realinteract.f90 \
		selfmolecule.f90 \
		Energy.f90 \
		displace_particle.f90 \
		rotate_particle.f90 \
		particle_swap_multiple.f90 \
		epsilon.f90 \
		nvt.f90 \
		main.f90 \
		particle.f90
	
	
#		transition_matrix.f90 \
		wang_landau.f90 \
		wltm.f90 \
		simulation_data.f90 \
		constants.f90       \
		util.f90 \
		cell.f90 \
		sampling.f90 \
		lj_potential.f90 \
		mc_moves.f90 \
		volume_move.f90 \
		npt.f90 


OBJS	      =	$(SRC:%.f90=%.o) 
# ran_uniform.o sstmm.o

%.o : %.mod

.SUFFIXES:  .f90 .o

.f90.o:

	$(FC) $(FFLAGS) -c $< -o $@

all:		$(PROGRAM)

$(PROGRAM)::	$(OBJS) $(MAKEFILE)
		@$(LD) $(LDFLAGS) $(OBJS) -o $(PROGRAM)

clean:;		@rm -f $(OBJS) $(PROGRAM) core *~ *.mod
