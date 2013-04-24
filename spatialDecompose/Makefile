# Defining variables ===============================
FC = gfortran
FCMPI = mpif90
MODFLAG = -J
FCFLAGS = -Wall $(MODFLAG)$(MODDIR)

SRCDIR = ./src/
#OBJDIR = ./obj/
MODDIR = ./mod/
#LIBDIR = $(SRCDIR)lib/
OUTDIR = ./out/

PROGRAM = spatialDecompose spatialDecompose_mpi

vpath %.f90 $(SRCDIR)
#vpath %.o $(SRCDIR)
#vpath % $(OUTDIR)

# Program ==================================
.PHONY : all 
all : $(PROGRAM)
# ---------------------------
spatialDecompose_mpi : spatialDecompose_mpi.o utility.o g96.o octave_save.o
	$(FCMPI) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F))

spatialDecompose_mpi.o : spatialDecompose_mpi.f90 utility.o g96.o octave_save.o
	$(FCMPI) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

spatialDecompose : spatialDecompose.o utility.o g96.o octave_save.o
	$(FC) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F))

spatialDecompose.o : spatialDecompose.f90 utility.o g96.o octave_save.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

# Library ===================================
octave_save.o : octave_save.f90 utility.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

g96.o : g96.f90 utility.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $<

utility.o : utility.f90
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 
#==================================================
.PHONY: clean
clean :
	rm -f $(SRCDIR)*.o $(MODDIR)*.mod
	rm -f $(addprefix $(OUTDIR), $(PROGRAM) )
	rm -f *~ $(SRCDIR)*~

