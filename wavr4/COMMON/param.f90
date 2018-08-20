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
  INTEGER(I4B) :: ne1,ne2,ne3, nexp, nbin1, nbin2
  INTEGER(I4B) :: nn(3), nsmax
  INTEGER(I4B) :: igq(3)
!  cut off related
  INTEGER(I4B) :: icut0, icut1, icut2, icut3, icut4
  INTEGER(I4B) :: imargin0,imargin1,imargin2,imargin3,imargin4
  REAL(DP) :: enzero,encut0,encut1,encut2,encut3,encut4
  REAL(DP) :: margin0,margin1,margin2,margin3,margin4
!   variable parameters
  REAL(DP) :: qe(6)
  INTEGER(I4B) :: na, jr, p, jp, namax
  INTEGER(I4B) :: j_parity, l_parity, jl_parity
!    flags:
  INTEGER(I4B) :: angular_problem_only_flag
  INTEGER(I4B) :: optimize_flag
  INTEGER(I4B) :: test_flag
  INTEGER(I4B) :: expansion_flag
  INTEGER(I4B) :: stage_flag
  INTEGER(I4B) :: oner_flag
!    auxiliary variables used (SAVEd) in Ar2HF potential.f90
  REAL(DP) :: r1q1,r1q2,r2q1,r2q2,r3q3,r4q3
END MODULE
