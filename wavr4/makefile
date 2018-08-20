## Sun specific modifications: 26/10/02
##
##  link make.inc to either make.win or make.unx
##
#include $(MAKEINC)
#!include make.win
include make.alpha
#include make.alpha.omp
#include make.sun32
#include make.sun64
#include make.sun64.omp

# make.inc must define the platform specific library calls,
# compile and link options. Variables which need to be defined:
# FC     f90 free source
# FC1    f77 fixed source (f90 might work as well)
# OPTS   high optimisation/less accurate compiler options
# OPTS0  low optimisation/more accurate compiler options
#        [primarily for computing a DVR grid. Since it's done
#         only once, it should not be performance sensitive]
# LIBS   libraries: specify where BLAS and LAPACK are
# OPTLIB any other optional libraries [empty]

# choose debug or release compilation
#OPTS  = $(OPT_DEB)
#OPTS0 = $(OPT_DEB)
OPTS  = $(OPT_REL)
OPTS0 = $(OPT_REL0)

##  define main program
MAIN = main

# =========================================================================
# list of object files (must be in the same order as source files)
# =========================================================================

OBJFILES1 = types.$(OBJ) \
            param.$(OBJ) \
            workarrays.$(OBJ) \
            drc3jj.$(OBJ) \
            la_auxmod.$(OBJ) \
            la_erinfo.$(OBJ) \
            la_f77_double.$(OBJ) \
            la_f95_double.$(OBJ) \
            la_dspev.$(OBJ) \
            la_dspevd.$(OBJ) \
            la_dspevx.$(OBJ) \
            la_dsyev.$(OBJ) \
            la_dsyevd.$(OBJ) \
            la_dsyevr.$(OBJ) \
            la_dsyevx.$(OBJ)

OBJFILES2 = base_lib.$(OBJ) \
            potential.$(OBJ) \
            laguer_grid.$(OBJ)

OBJFILES3 = grid3d.$(OBJ) \
            laguer_dvr.$(OBJ) \
            expansion.$(OBJ)

OBJFILES4 = angles.$(OBJ) \
            mainj.$(OBJ)

OBJFILES = $(OBJFILES1) \
           $(OBJFILES2) \
           $(OBJFILES3) \
           $(OBJFILES4) \
           $(MAIN).$(OBJ)


# =========================================================================
# compilation and linking
# =========================================================================

$(MAIN)$(EXE): $(OBJFILES)
	$(LNK) $(OUTL) $(OPTOMP) $(OBJFILES) $(OPTLIB) $(LIBS)

#	$(LNK) $(OUTL)$(MAIN)$(EXE) $(OBJFILES) $(OPTLIB) $(LIBS)
# =========================================================================
types.$(OBJ):         $(COMMON)types.f90
	$(FC) $(OPTS) $(COMMON)types.f90

param.$(OBJ):         $(COMMON)param.f90
	$(FC) $(OPTS) $(COMMON)param.f90

workarrays.$(OBJ):    $(CORE)workarrays.f90
	$(FC) $(OPTS) $(OPTOMP) $(CORE)workarrays.f90

drc3jj.$(OBJ):         $(COMMON)drc3jj.f
	$(FC1) $(OPTS) $(COMMON)drc3jj.f

la_auxmod.$(OBJ):     $(LAPACK95)la_auxmod.f90
	$(FC) $(OPTS) $(LAPACK95)la_auxmod.f90

la_erinfo.$(OBJ):     $(LAPACK95)la_erinfo.f90
	$(FC) $(OPTS) $(LAPACK95)la_erinfo.f90

la_f77_double.$(OBJ): $(LAPACK95)la_f77_double.f90
	$(FC) $(OPTS) $(LAPACK95)la_f77_double.f90

la_f95_double.$(OBJ): $(LAPACK95)la_f95_double.f90
	$(FC) $(OPTS) $(LAPACK95)la_f95_double.f90

la_dspev.$(OBJ):      $(LAPACK95)la_dspev.f90
	$(FC) $(OPTS) $(LAPACK95)la_dspev.f90

la_dspevd.$(OBJ):     $(LAPACK95)la_dspevd.f90
	$(FC) $(OPTS) $(LAPACK95)la_dspevd.f90

la_dspevx.$(OBJ):     $(LAPACK95)la_dspevx.f90
	$(FC) $(OPTS) $(LAPACK95)la_dspevx.f90

la_dsyev.$(OBJ):      $(LAPACK95)la_dsyev.f90
	$(FC) $(OPTS) $(LAPACK95)la_dsyev.f90

la_dsyevd.$(OBJ):     $(LAPACK95)la_dsyevd.f90
	$(FC) $(OPTS) $(LAPACK95)la_dsyevd.f90

la_dsyevr.$(OBJ):     $(LAPACK95)la_dsyevr.f90
	$(FC) $(OPTS) $(LAPACK95)la_dsyevr.f90

la_dsyevx.$(OBJ):     $(LAPACK95)la_dsyevx.f90
	$(FC) $(OPTS) $(LAPACK95)la_dsyevx.f90

# =========================================================================
base_lib.$(OBJ):      $(COMMON)base_lib.f90 $(OBJFILES1)
	$(FC) $(OPTS) $(COMMON)base_lib.f90

# use FC if it's free F90 and FC1 if it's fixed style
potential.$(OBJ):     potential.f90 $(OBJFILES1)
	$(FC) $(OPTS) $(OPTOMP) potential.f90

laguer_grid.$(OBJ):    $(COMMON)laguer_grid.f90 $(OBJFILES1)
	$(FC) $(OPTS0) $(COMMON)laguer_grid.f90

# =========================================================================
grid3d.$(OBJ):         $(COMMON)grid3d.f90 $(OBJFILES1) $(OBJFILES2)
	$(FC) $(OPTS0) $(COMMON)grid3d.f90

laguer_dvr.$(OBJ):     $(COMMON)laguer_dvr.f90 $(OBJFILES1) $(OBJFILES2)
	$(FC) $(OPTS0) $(COMMON)laguer_dvr.f90

expansion.$(OBJ):     $(COMMON)expansion.f90 $(OBJFILES1) $(OBJFILES2)
	$(FC) $(OPTS) $(OPTOMP) $(COMMON)expansion.f90

# =========================================================================
angles.$(OBJ):        $(COMMON)angles.f90 $(OBJFILES1) $(OBJFILES2) $(OBJFILES3)
	$(FC) $(OPTS) $(OPTOMP) $(COMMON)angles.f90

mainj.$(OBJ):         $(CORE)mainj.f90 $(OBJFILES1) $(OBJFILES2) $(OBJFILES3)
	$(FC) $(OPTS) $(CORE)mainj.f90

# =========================================================================
$(MAIN).$(OBJ):       $(CORE)$(MAIN).f90 $(OBJFILES1) $(OBJFILES2) $(OBJFILES3) $(OBJFILES4)
	$(FC) $(OPTS) $(OPTOMP) $(CORE)$(MAIN).f90

#.$(SUF).$(OBJ): modules.$(OBJ)
#	$(FC) $(OPTS) $<

all: $(MAIN)$(EXE)
	$(MAIN)$(EXE)
# <input.txt >output.txt

clean:
	$(DEL) *.$(OBJ) *.mod
#  *.dat

zip: clean
	zip -r -9 $(MAIN).zip *
