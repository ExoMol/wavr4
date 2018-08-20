MODULE param
  USE types
  IMPLICIT NONE
!   read in parameters
  INTEGER(I4B) :: opt
  REAL(DP) :: mass1,mass2,mass3,mass4
  REAL(DP) :: muq(3), mur,mud1,mud2, msum
  REAL(DP) :: re(3), we(3), De(3)
  INTEGER(I4B) :: lmax,kmax,jmax, jrmax, krmax, mmax
  INTEGER(I4B) :: j_parity_max, l_parity_max, jl_parity_max
  INTEGER(I4B) :: nt1,nt,nphi,nagrid
  INTEGER(I4B) :: ne1,ne2,ne3
  INTEGER(I4B) :: nn(3), nsmax
  INTEGER(I4B) :: igq(3)
!   variable parameters
  REAL(DP) :: qe(6)
  INTEGER(I4B) :: namax
!  INTEGER(I4B) :: na, jr, p, jp, namax
!  INTEGER(I4B) :: j_parity, l_parity, jl_parity

  REAL(DP) :: enzero

!    flags:
  INTEGER(I4B) :: angular_problem_only_flag
  INTEGER(I4B) :: optimize_flag
  INTEGER(I4B) :: test_flag
  INTEGER(I4B) :: expansion_flag

! HHSYM indicates extra (H-H) permutation symmetry
!  LOGICAL(LGT), PARAMETER :: HHSYM = .TRUE.
  LOGICAL(LGT), PARAMETER :: HHSYM = .FALSE.

END MODULE
