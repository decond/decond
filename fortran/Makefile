# ===== Basic modification =====
# Hopefully it will be enough by just modifying this block

# MPI Fortran compiler, use h5pfc if available
FC = h5pfc

# -module for Intel ifort, -J for GNU gfortran, -M for Fujitsu frt
MODFLAG = -module

# Compiler flags:
# -Kfast for Fujitsu compiler (K computer)
# -DFORTRAN2008 for compilers support Fortran 2008 specification
FCFLAG =

# XTC Library directory (where lib/ and include/ are located inside)
XDRFILE_DIR = $(HOME)/local

# HDF5 Library directory (where lib/ and include/ are located inside)
HDF5_DIR = $(HOME)/local

# FFTW3 Library directory (where lib/ and include/ are located inside)
FFTW_DIR = $(HOME)/local


# ===== Advanced modification =====
# Modify this block if the above modification only is not enough

# XTC Library
XDRFILE_LINK = $(XDRFILE_DIR)/lib/libxdrfile.a
# XDRFILE_LINK = -L$(XDRFILE_DIR)/lib -lxdrfile

# HDF5 Library
HDF5_INCL = -I$(HDF5_DIR)/include
# Manually setup HDF5_LINK if FC is not h5pfc
ifneq ($(FC), h5pfc)
    HDF5_LINK = $(HDF5_DIR)/lib/libhdf5hl_fortran.a \
                $(HDF5_DIR)/lib/libhdf5_fortran.a \
                $(HDF5_DIR)/lib/libhdf5_hl.a \
                $(HDF5_DIR)/lib/libhdf5.a \
                -lz -ldl -lm -Wl,-rpath -Wl,$(HDF5_DIR)/lib
endif

# FFTW3 Library
FFTW_INCL = -I$(FFTW_DIR)/include
FFTW_LINK = $(FFTW_DIR)/lib/libfftw3.a
# FFTWLINK = -L$(FFTW_DIR)/lib -lfftw3


# ===== Basically NO modification is needed below this line =====
SRCDIR = ./src
BINDIR = ../bin

FCFLAGS = $(MODFLAG) $(SRCDIR) $(HDF5_INCL) $(FFTW_INCL) $(FCFLAG)
# FCFLAGS = $(MODFLAG) $(SRCDIR) $(HDF5_INCL) $(FFTW_INCL) -g -check all -traceback

FLFLAGS = $(HDF5_LINK) $(FFTW_LINK) $(XDRFILE_LINK)


# ===== Program =====
DEC = decond
TRJCONV = trjconv2com
PROGRAM = $(TRJCONV) $(DEC)

TRJCONV_DEPS = xdr.o top.o utility.o varpars.o
TRJCONV_OBJS = trjconv2com.o $(TRJCONV_DEPS)

DECOND_DEPS = manager.o mpiproc.o spatial_dec.o energy_dec.o correlation.o utility.o xdr.o top.o fftw.o varpars.o xyz.o
DECOND_OBJS = decond.o $(DECOND_DEPS)

# ---------------------------
.PHONY : all
all : $(addprefix $(BINDIR)/, $(PROGRAM))

$(BINDIR)/$(TRJCONV) : $(addprefix $(SRCDIR)/, $(TRJCONV_OBJS))
	$(FC) $^ $(FLFLAGS) -o $@ 

$(BINDIR)/$(DEC) : $(addprefix $(SRCDIR)/, $(DECOND_OBJS))
	$(FC) $^ $(FLFLAGS) -o $@


# ---------------------------
$(SRCDIR)/trjconv2com.o : $(addprefix $(SRCDIR)/, $(TRJCONV_DEPS))

$(SRCDIR)/decond.o : $(addprefix $(SRCDIR)/, $(DECOND_DEPS))


# ===== Library =====
$(SRCDIR)/manager.o : $(addprefix $(SRCDIR)/, mpiproc.o spatial_dec.o energy_dec.o top.o utility.o varpars.o xdr.o correlation.o xyz.o)

$(SRCDIR)/spatial_dec.o : $(addprefix $(SRCDIR)/, mpiproc.o top.o utility.o)

$(SRCDIR)/energy_dec.o : $(addprefix $(SRCDIR)/, mpiproc.o utility.o)

$(SRCDIR)/mpiproc.o : $(addprefix $(SRCDIR)/, varpars.o)

$(SRCDIR)/xdr.o : $(addprefix $(SRCDIR)/, varpars.o utility.o)

$(SRCDIR)/xyz.o : $(addprefix $(SRCDIR)/, varpars.o utility.o)

$(SRCDIR)/varpars.o : $(SRCDIR)/top.o

$(SRCDIR)/top.o : $(SRCDIR)/utility.o

$(SRCDIR)/correlation.o : $(SRCDIR)/fftw.o


# ===================
$(SRCDIR)/%.o : $(SRCDIR)/%.F90
	$(FC) -c $(FCFLAGS) $< -o $@


#==================================================
.PHONY: clean
clean :
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod
	rm -f $(addprefix $(BINDIR)/, $(PROGRAM) )

