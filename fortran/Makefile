# Defining variables ===============================
FC = gfortran
FCMPI = mpif90
OMPFLAGS = -fopenmp
MODFLAG = -J

SRCDIR = ./src/
#OBJDIR = ./obj/
MODDIR = ./mod/
#LIBDIR = $(SRCDIR)lib/
OUTDIR = ./out/

FCFLAGS = -Wall $(MODFLAG)$(MODDIR) -I$(MODDIR) $(OMPFLAGS) -O2

PROGRAM = spatialDecompose spatialDecompose_mpi spatialDecompose_openmp bridge com corr_test

vpath %.f90 $(SRCDIR)
#vpath %.o $(SRCDIR)
#vpath % $(OUTDIR)

# Program ==================================
.PHONY : all 
all : $(PROGRAM)

# ---------------------------
spatialDecompose_openmp : spatialDecompose_openmp.o utility.o g96.o octave_save.o
	$(FC) $(OMPFLAGS) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F))

spatialDecompose_openmp.o : spatialDecompose_openmp.f90 utility.o g96.o octave_save.o
	$(FCMPI) $(OMPFLAGS) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

spatialDecompose_mpi : spatialDecompose_mpi.o utility.o g96.o octave_save.o
	$(FCMPI) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F))

spatialDecompose_mpi.o : spatialDecompose_mpi.f90 utility.o g96.o octave_save.o
	$(FCMPI) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

spatialDecompose : spatialDecompose.o utility.o g96.o octave_save.o
	$(FC) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F))

spatialDecompose.o : spatialDecompose.f90 utility.o g96.o octave_save.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

bridge: bridge.o utility.o g96.o octave_save.o
	$(FC) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F))

bridge.o : bridge.f90 utility.o g96.o octave_save.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

com: com.o utility.o trr.o octave_save.o
	$(FC) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F))

com.o : com.f90 utility.o trr.o octave_save.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

corr_test: corr_test.o correlation.o fftw.o
	$(FC) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F)) -L$(HOME)/local/fftw-3.3.2/lib -lfftw3 -lm

corr_test.o : corr_test.f90 
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

# Library ===================================
octave_save.o : octave_save.f90 utility.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

g96.o : g96.f90 utility.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $<

trr.o : trr.f90 utility.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $<

utility.o : utility.f90
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

correlation.o : correlation.f90
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

fftw.o : fftw.f90
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< -I$(HOME)/local/fftw-3.3.2/include
#==================================================
.PHONY: clean
clean :
	rm -f $(SRCDIR)*.o $(MODDIR)*.mod
	rm -f $(addprefix $(OUTDIR), $(PROGRAM) )
	rm -f *~ $(SRCDIR)*~

