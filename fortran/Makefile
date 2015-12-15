# ===== Essential modification =====
FC = ifort
H5PFC = h5pfc
FCMPI = mpif90
MODFLAG = -module

## XTC Library
XDRFILEDIR = $(HOME)/local
XDRFILELINK = -L$(XDRFILEDIR)/lib -lxdrfile

## HDF5 Library
HDF5_INCLUDE_DIR = $(HOME)/local/include

## FFTW3 Library
FFTWDIR = $(HOME)/local
FFTWLINK = $(FFTWDIR)/lib/libfftw3.a
## For dynamic linking
# FFTWLINK = -L$(FFTWDIR)/lib -lfftw3 -lm


# ===== Advanced modification =====
SRCDIR = ./src/
MODDIR = ./mod/
OUTDIR = ./out/
BINDIR = ../bin/

FCFLAGS = $(MODFLAG) $(MODDIR) -O2 -heap-arrays 32768
## For debug
# FCFLAGS = $(MODFLAG) $(MODDIR) -g -check all -traceback


# ===== Basically NO modification is needed below this line =====
PROGRAM = decond trjconv2com
DIR = $(MODDIR) $(OUTDIR)
vpath %.f90 $(SRCDIR)

# ===== Program =====
.PHONY : all dir
all : dir $(PROGRAM)
	ln -sf ../fortran/$(OUTDIR)* $(BINDIR)

dir :
	mkdir -p $(DIR)

# ---------------------------
trjconv2com : utility.o xdr.o top.o trjconv2com.o
	$(FC) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F)) $(XDRFILELINK)

trjconv2com.o : trjconv2com.f90
	$(FC) $(FCFLAGS) -c $< -o $(SRCDIR)$(@F)

decond : mpiproc.o spatial_dec.o energy_dec.o correlation.o utility.o xdr.o top.o decond.o
	$(H5PFC) -o $(OUTDIR)$(@F) $(addprefix $(SRCDIR), $(^F)) $(FFTWLINK) $(XDRFILELINK)

decond.o : decond.f90
	$(H5PFC) $(FCFLAGS) -c $< -o $(SRCDIR)$(@F)

# ===== Library =====
spatial_dec.o : spatial_dec.f90 mpiproc.o top.o utility.o
	$(FCMPI) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $<

energy_dec.o : energy_dec.f90 mpiproc.o utility.o
	$(FCMPI) $(FCFLAGS) -I$(HDF5_INCLUDE_DIR) -o $(SRCDIR)$(@F) -c $<

mpiproc.o : mpiproc.f90
	$(FCMPI) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $<

top.o : top.f90 utility.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $<

xdr.o : xdr.f90 utility.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $<

utility.o : utility.f90
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

correlation.o : correlation.f90 fftw.o
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< 

fftw.o : fftw.f90
	$(FC) $(FCFLAGS) -o $(SRCDIR)$(@F) -c $< -I$(FFTWDIR)/include
#==================================================
.PHONY: clean
clean :
	rm -f $(SRCDIR)*.o $(MODDIR)*.mod
	rm -f $(addprefix $(OUTDIR), $(PROGRAM) )
	rm -f *~ $(SRCDIR)*~
	rm -f $(addprefix $(BINDIR), $(PROGRAM) )

