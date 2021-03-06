# WAVR4

Wide Amplitude Vibration-Rotation 4-atomic code

## Quick Start

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


## Contains

CORE files:

12-3/     truncation order (q1,q2),q3
          main.f90 [h6 in stage 3 is prepared on disk, i.e. h6 swapped]

3-12/     not present in this version

COMMON/   files common to 12-3 and 3-12

LAPACK95/ standard f95 interfaces for LAPACK (win,unix)

LAPACK95.sun64/   f95 interfaces using 64-bit integers
          to be used when compiled in 64 bit using Sun performance lib


## Full Description


IN Kozin, MM Law, JM Huston, J Tennyson, JCP 118, 4896 (2003);
IN Kozin, MM Law, J Tennyson, JM Huston,  CPC XXX, XXX (200X).

The program computes ro-vibrational energy levels together with
eigenfunctions (if requested) of a four-atomic molecule.
The Hamiltonian employs polyspherical coordinates, see
M. Mladenovic JCP 112, p1070 (2000) [1].
The primary basis is a product of coupled 3D angular basis [1] and
three DVR functions for each radial coordinate. The DVR functions
are based on Morse-Laguerre or spherical oscillator functions,
see J. Tennyson & B.T. Sutcliffe JCP 77, p4065 (1982) and
JMS 101, 71 (1983) [2].

There are 4 choices of internal coordinates to choose from [1].
For atom numbering refer to the figure below (adapted from [1]):

```
1. Jacobi vectors          2. Radau vectors

 (1)    (4)                      (1)
  ^      ^                        ^ q1
  | q1   | q3                 q2  |
  |      |                 (2)<---* (4)
  *------*->(3)                   |
  |   q2                          v q3
 (2)                             (3)

3. Diatom-diatom vectors   4. Orthogonalised satellite vectors

 (1)        (4)              (1)     (2)
  ^          ^                 ^     ^
  | q1       | q2            q1 \   / q2
  |          |                   \ /
  *--------->*                    V
  |    q3    |
 (2)        (3)              (3)---->(4)
                                  q3
```

The inversion symmetry is always present and implemented.
Permutational symmetries which result in vector inversion,
(e.g. P(1-2) in case 1 above which transforms vector q1 to -q1)
can be enabled or disabled via input file (see below).
Permutation symmetries involving two or more vectors are
not implemented in this version (e.g. P(1-4) in case 3).


## ALGORITHM

```
Initiation: read input file, generate 3D DVR grid (q1,q2,q3), compute
         expansion of potential for every grid point if necessary,
         allocate work arrays, read in eigenvectors if available.
Stage 1: solve 3D FBR problem for every 3D DVR grid point.
Stage 2: solve 5D problem (3D+(q1,q2)) for every q3 grid point.
Stage 3: solve 6D problem (5D+q3).
         End of pure vibrational job otherwise repeat for every K.
Stage 4 (final): solve ro-vibration problem
```

## USER INPUT

the program can be typically run as

$ main.exe <input.txt >output.txt

OR simply

$ main.exe

if input.txt and output.txt are assigned to units 5 & 6 in main.f90

The user must modify input.txt file. The format is a free form
but line positions are fixed. Comment lines start with "!".
Almost all parameters are explained in the comment lines.

```
LINE  1: title
         -- A80 format, 80-character title
LINE  2: !
LINE  3: mass1,mass2,mass3,mass4
         -- floats, nuclear masses in atomic units
LINE  4: !
LINE  5: opt
         -- integer, defines coordinate system: 
         = 1, Jacobi; 
         = 2, Radau; 
         = 3, diatom-diatom;
         = 4, orthogonal-satellite vectors.
         See figure above or Fig. 1 in Ref [1]
LINE  6: !
LINE  7: igq(1),re(1), we(1), De(1), nn(1)
LINE  8: igq(2),re(2), we(2), De(2), nn(2)
LINE  9: igq(3),re(3), we(3), De(3), nn(3)
         -- parameters for the radial basis in q1, q2, and q3.
         igq(1:3), integer: 
         = 1, Morse oscilator-like functions;
         = 2, spherical oscillator functions.
         re(1:3), float: r_e for Morse oscilator-like functions
         and eta for spherical oscillator functions.
         we(1:3), float: w_e.
         De(1:3), float: D_e (not used for spherical harmonics).
         nn(1:3), integer: number of functions (i.e.\ grid size).
         NOTE: If the nn(i) = 1 then the i-th coordinate is frozen
         to the reference value given below (see LINE 26)
LINE 10: !
LINE 11: angular_problem_only_flag, optimize_flag, test_flag, expansion_flag
         -- integer flags.
         angular_problem_only_flag:
         = 0 (default), normal operation;
         > 0, only 3D angular problem is computed for the reference radial
         triple point specified in LINE 26.
         optimize_flag: not used.
         test_flag:
         = 0, default;
         = 1 or 2, controls output level
         giving more and more test output.
         expansion_flag:
         = 0, compute 3D FBR angular matrix elements using quadrature
         integration;
         = 1, compute angular matrix
         elements through the expansion of the angular potential at
         every radial triple grid point (stored as array)
         V= SUM a_{j,l,k}*C(j,l,k,K=0,kappa=-1)
         where C(j,l,k,K,kappa) are basis functions.
         Then the integrals of V can be expressed using 3J-symbols.
LINE 12: !
LINE 13: jmax, lmax, kmax, jrmax, krmax, j_parity_max, l_parity_max, jl_parity_max
         -- all integers, define the bending-rotation basis.
         the bending (angular) basis (j,l,k), rotational basis (J,K)
         jmax: maximum j.
         lmax: maximum l.
         kmax: maximum k.
         jrmax: maximum J.
         krmax: maximum K.
         symmetry of the system (default 0 0 0)
         j-parity = 1    if q1 -> -q1  (e.g. P(HH)   in H2--He2)
         l-parity = 1    if q2 -> -q2  (e.g. P(HeHe) in H2--He2)
         (j+l)-parity =1 if q3 -> -q3  (e.g. P(CC)   in H-CC-H)
LINE 14: !
LINE 15: ne1, ne2, ne3
         -- all integers, specify maximum quantum numbers for angular 
         expansion in j, l, and k respectively.
         (typically doubled numbers of the maximum quantum
         numbers specified for the angular basis)
LINE 16: !
LINE 17: nt1, nt, nphi, iauto
         -- all integers. The first three numbers specify the size of the
         angular quadrature in theta_1, theta and phi respectively
         (corresponding to j, l, and k respectively).
         iauto: 
         = 0, WAVR4 uses user supplied numbers;
         = 1, 2, or 3, WAVR4 automatically defines minimum, double or
         quadruple grid sizes as described in the notes below.
LINE 18: !
LINE 19: enzero
         -- zero point energy
LINE 20: !
LINE 21: icut0, icut1, icut2, icut3, icut4
         -- all integers.
         icutN defines the number of eigenvectors computed
         during every diagonalisation at stage N.
         icutN: 
         = 0, all eigenvectors are found.
         > 0, only icutN lowest eigenvectors are found.
         < 0, encutN controls the number of computed eigenvalues 
         and defines maximum energy for computed eigenvalues.
LINE 22: encut0, encut1, encut2, encut3, encut4
         -- all floats.
         encutN defines the energy cutoff used at stage N.
         The selection based on energy is always applied.
         To disable it, set the respective energy to a large number.
         Note 1: stage 0 is reserved for future use.
         Note 2: if your cutoff is based on the number of levels
         then encut1 and encut2 are still in effect.
         This is useful for removing very high eigenstates
         (particularly those on the edges of the radial grid).
         WARNING: special considerations regarding stages 3 and 4.
         If you set icut3/4 = 0 then all eigenvalues will be computed
         but no eigenvectors.  (use for test purposed only).
         Most importantly, because hamiltonian matrices in stages 3/4
         are created in a packed form, storage for eigenvectors must be
         allocated separately. In practice this means that for a start,
         it is safer to use truncation based on numbers
         (i.e.\ setting positive icut3/4) rather than based on
         energy encut3/4. If icut3/4 are negative,
         abs(icut3/4) are used to reserve the storage for eigenvectors
         and it could happen that there are more than
         abs(icut3/4) states below encut3/4.
         However it is still possible to use the energy cut off:
         icut3/4 must be sufficiently big but not too big to make
         the packed form inefficient (i.e. less than 10% of the matrix size).
LINE 23: margin0, margin1, margin2, margin3, margin4
         -- all floats. The only used now is margin1 which is used to
         skip 3D angular calculation at a radial triple point
         if minimum potential energy computed on the angular grid is higher
         than margin1.
LINE 24: imargin0, imargin1, imargin2, imargin3, imargin4
         -- all integers, not used.
LINE 25: !
LINE 26: q1,q2,q3,theta1,theta,phi
         -- all floats. These define some reference configuration
         (e.g.\ equilibrium).
LINE 27: !
LINE 28: stage_flag
         -- integer, from 0 (default) to 4. If the number is greater
         then zero it indicates up to which step the eigenvectors are
         already available and so the program will attempt to read these 
         eigenvectors. It is the responsibility of the user to make sure
         that the eigenvectors correspond to the same basis and potential.
LINE 29: !
LINE 30: oner_flag
         -- integer.
         = 0, 1/q3^2 computed in DVR approximation
         = 1, 1/q3^2 computed in analytically
         (this is only valid for igq(3) = 2, see LINES 7-9).
```

## Notes

Angular potential expansion size must be present but it is only
used when needed.

```
stage 1:  fixed q1,q2,q3,K (i.e. pure bending problem)
stage 2:  fixed q3,K  or q1,q2,K
stage 3:  fixed K
stage 4:  final
```

When "iauto" is 0 the angular quadrature size is determined by the
number of points specified. If iauto is equal to 1 the program uses
the minimum numbers of grid points. These are printed in the output.

```
iauto = 0 : the grid size to evaluate the integrals or expansion is
            determined by the input (whatever it is)

iauto > 0 : the grid size is taken to be dependent on n1, n2, n3
            n1 is jmax or ne1 (expansion size)
            n2 is lmax or ne2 (expansion size)
            n3 is kmax or ne3 (expansion size)

iauto = 1 : minimal (single) grid size
            theta1: nt1  = n1+1;
            theta : nt2  = n2+1;
            phi   : nphi = (n3+1)/2  if n3 odd;
                           n3/2+1    if n3 is even

iauto = 2 : double grid size
            theta1: nt1  = 2*n1+1;
            theta : nt2  = 2*n2+1;
            phi   : nphi = n3+1

iauto = 3 : quadruple grid size
            theta1: nt1  = 4*n1+1;
            theta : nt2  = 4*n2+1;
            phi   : nphi = 2*n3+1
```
