##
##  link make.inc to either make.win or make.unx
##
#!include make.win
#include make.sun64.omp
#include make.sun64
#include make.alpha
include make.intel
#include $(MAKEINC)

# choose debug or release compilation
# OPTS = $(OPT_DEB)
OPTS = $(OPT_REL)

##  define main program
MAIN = main

# =========================================================================
# do not depend on any module
#
OBJFILES1 = types.$(OBJ) param.$(OBJ)
#
# depend on OBJFILES1
#
OBJFILES2 = base_lib.$(OBJ)
#
# depend on OBJFILES2
#
OBJFILES3 = workarrays.$(OBJ)
#
# depend on OBJFILES3
#
OBJFILES4 = angles.$(OBJ)
#
# depend on all
#
OBJFILES = $(MAIN).$(OBJ) $(OBJFILES1) $(OBJFILES2) $(OBJFILES3) $(OBJFILES4)
#
# =========================================================================
#
# list of f90 files in compilation order
F90FILES1 = types.f90 param.f90
F90FILES2 = base_lib.f90
F90FILES3 = workarrays.f90
F90FILES4 = angles.f90
F90FILES = $(F90FILES1) $(F90FILES2) $(F90FILES3) $(F90FILES4)  $(MAIN).f90

$(MAIN)$(EXE): $(OBJFILES)
	$(LNK) $(OUTL) $(OPTOMP) $(OBJFILES) $(OPTLIB) $(LIBS)

# =========================================================================
types.$(OBJ): types.f90
	$(FC) $(OPTS) types.f90

# =========================================================================
param.$(OBJ):       param.f90
	$(FC) $(OPTS) param.f90

# =========================================================================
base_lib.$(OBJ):    $(OBJFILES1) base_lib.f90
	$(FC) $(OPTS) base_lib.f90

# =========================================================================
workarrays.$(OBJ):    $(OBJFILES1) workarrays.f90
	$(FC) $(OPTS) workarrays.f90

# =========================================================================
angles.$(OBJ): $(OBJFILES1) $(OBJFILES2) $(OBJFILES3) angles.f90
	$(FC) $(OPTS) angles.f90

# =========================================================================
$(MAIN).$(OBJ): $(OBJFILES1) $(OBJFILES2) $(OBJFILES3) $(OBJFILES4) $(MAIN).f90
	$(FC) $(OPTS) $(OPTOMP) $(MAIN).f90

#.$(SUF).$(OBJ): modules.$(OBJ)
#	$(FC) $(OPTS) $<

all: $(MAIN)$(EXE)
	$(MAIN)$(EXE)
# <input.txt >output.txt

clean:
	$(DEL) *.$(OBJ) *.mod $(MAIN)$(EXE) debug.pdb
#  *.dat

zip: clean
	zip -9 $(MAIN).zip *

