WAVR4  QUICK START

All platform/project dependent settings should be in make.* file

1. Edit your make.* file
   1) Define compiler, linker and options. Note that
      laguer* files should be compiled without compromising
      accuracy [controlled by OPT_DEB0 and OPT_REL0 variables]
   2) Adjust variables:
      COMMON   (default COMMON/ ),
      LAPACK95 (default LAPACK95/ ) and
      CORE (either 3-12/ or 12-3/ depending on the truncation order)
   3) make sure you have BLAS and LAPACK (current version 3.1),
      and the respective link options are correct

2. Copy your potential.f90 file to the root dir
   The potential function can be anything as long as it is
   in the form v(q1,q2,q3,theta1,theta,phi) and
   compatible with the interface provided in the sample files.

3. Run make (run "make clean" first when attemting to re-compile)

4. Prepare the input file (see 00docs.txt) and run your job

NOTE: direct calls to BLAS were usually avoided in the code
because they make the code less clear and I did not notice a
big gain in performance on Alpha. Instead use compiler option
to replace Fortran 90/95 intrinsic matmul and dot_product by
BLAS calls during compilation time (if possible).


CONTAINS:

CORE files:

12-3/     truncation order (q1,q2),q3
          main.f90 [h6 in stage 3 is prepared on disk, i.e. h6 swapped]

3-12/     not present in this version

COMMON/   files common to 12-3 and 3-12

LAPACK95/ standard f95 interfaces for LAPACK (win,unix)

LAPACK95.sun64/   f95 interfaces using 64-bit integers
          to be used when compiled in 64 bit using Sun performance lib
