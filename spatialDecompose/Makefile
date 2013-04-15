# Defining variables ===============================
FC = gfortran
MODFLAG = -J
FCFLAGS = $(MODFLAG)$(MODDIR)

SRCDIR = ./src/
#OBJDIR = ./obj/
MODDIR = ./mod/
#LIBDIR = $(SRCDIR)lib/
OUTDIR = ./out/

PROGRAM = spatialDecompose

vpath %.f90 $(SRCDIR)
#vpath %.o $(SRCDIR)
#vpath % $(OUTDIR)

# Program ==================================
.PHONY : all 
all : $(PROGRAM)
# ---------------------------
spatialDecompose : spatialDecompose.o utility.o g96.o octave_save.o
	$(FC) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F))

spatialDecompose.o : spatialDecompose.f90 utility.o g96.o octave_save.o
	$(FC) -o $(SRCDIR)$(@F) -c $< $(FCFLAGS)

# Library ===================================
octave_save.o : octave_save.f90 utility.o
	$(FC) -o $(SRCDIR)$(@F) -c $< $(FCFLAGS)

g96.o : g96.f90 utility.o
	$(FC) -o $(SRCDIR)$(@F) -c $< $(FCFLAGS)

utility.o : utility.f90
	$(FC) -o $(SRCDIR)$(@F) -c $< $(FCFLAGS)
#==================================================
.PHONY: clean
clean :
	rm -f $(SRCDIR)*.o $(MODDIR)*.mod
	rm -f $(addprefix $(OUTDIR), $(PROGRAM) )
	rm -f *~ $(SRCDIR)*~

