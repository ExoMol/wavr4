!  potential function for acetylene based on
!  Halonen, Child, Carter, Mol. Phys. 1982 v47 p1097
!  which is an update for potetial1.f90 based on
!  Carter, Mills, Murrell, Mol. Phys. 1980 v41 p191
!
!  INPUT: diatom-diatom orthogonal vectors, valence-like
!  OUTPUT: energy (in cm-1)
!
!  Intermediate step: interatomic distances are computed
!                     using formulae of MM supplement
!
!  H(1)                 H(4)
!     \ q1             /
!      \      q3      / q2
!       \------------/
!        C(2)       C(3)
!
!  21/11/02  I.K.
!
!  experimental minimum is
!  R(CC)= 1.203A; R(CH)= 1.061A; De= 17.55
!  diat-diat: q3=0; q1=1.203; q2=3.325; theta1=any theta=any phi=any
!  orth-sat:  q3=1.203; q1=q2=1.6625; theta1=Pi; theta0; phi=0 (in fact any)
!  d-d valence: q1=q2=1.061; q3=1.366; theta1=Pi; theta0; phi= any
!                            1.203+0.163 (because of CM shifts)

MODULE potential

INTERFACE v
  FUNCTION v(q1,q2,q3,theta1,theta,phi)
  USE types
  USE param
  REAL(DP), INTENT(IN) :: q1,q2,q3,theta1,theta,phi
  REAL(DP) :: v
  END FUNCTION v
END INTERFACE

END MODULE

FUNCTION v(q1,q2,q3,theta1,theta,phi)
  USE types
IMPLICIT NONE
  REAL(DP) :: v
  REAL(DP), INTENT(IN) ::  q1,q2,q3,theta1,theta,phi
REAL(DP) :: q1x,    q1z, q2x,q2y,q2z,         q3z  ! q1y, q3x,q3y,
! REAL(DP) :: tmp1,tmp2, stheta1,ctheta1,stheta,ctheta,sphi,cphi
REAL(DP) :: r1x,r1y,r1z, r2x,r2y,r2z, r3x,r3y,r3z, r4x,r4y,r4z
REAL(DP) :: r12,r13,r14,r23,r24,r34
! REAL(DP), SAVE :: a(4,3)
REAL(DP) :: a(4,3)
REAL(DP) :: ev2cm = 8065.540928841530_dp
INTEGER(I4B) :: aflag = 0
!$OMP THREADPRIVATE(aflag,a)
SAVE a, aflag

if(aflag == 0) then
  call ac(a)
  aflag = 1
end if

q1x = q1*sin(theta1)
! q1y = 0.0_dp
q1z = q1*cos(theta1)

q2x = q2*sin(theta)*cos(phi)
q2y = q2*sin(theta)*sin(phi)
q2z = q2*cos(theta)

! q3x = 0.0_dp
! q3y = 0.0_dp
q3z = q3

! this makes rather sparse matrix
! q1x  q2x   0
!  0   q2y   0
! q1z  q2z  q3z

! some terms removed because a(1,2)=0

r1x = a(1,1)*q1x
r1y =                                   0.0_dp
r1z = a(1,1)*q1z           +a(1,3)*q3z

! some terms removed because a(2,2)=0

r2x = a(2,1)*q1x
r2y =                                   0.0_dp
r2z = a(2,1)*q1z           +a(2,3)*q3z

! some terms removed because a(3,1)=0

r3x =            a(3,2)*q2x
r3y =            a(3,2)*q2y
r3z =            a(3,2)*q2z+a(3,3)*q3z

! some terms removed because a(4,1)=0

r4x =            a(4,2)*q2x
r4y =            a(4,2)*q2y
r4z =            a(4,2)*q2z+a(4,3)*q3z

r12 = sqrt((r1x-r2x)**2+(r1y-r2y)**2+(r1z-r2z)**2)
r13 = sqrt((r1x-r3x)**2+(r1y-r3y)**2+(r1z-r3z)**2)
r14 = sqrt((r1x-r4x)**2+(r1y-r4y)**2+(r1z-r4z)**2)
r23 = sqrt((r2x-r3x)**2+(r2y-r3y)**2+(r2z-r3z)**2)
r24 = sqrt((r2x-r4x)**2+(r2y-r4y)**2+(r2z-r4z)**2)
r34 = sqrt((r3x-r4x)**2+(r3y-r4y)**2+(r3z-r4z)**2)

!===================================================================
! 21/11/02:  Let's take  H1 -> h, H4 -> h1, C2 -> c, C3 -> c1, THEN
! the internuclear distances from S. Carter et al paper are
!      H1       r_hh: r6   r_ch: r1   r4   r5   r2   r_cc: r3
! r1  /| \ r4         r14        r12  r13  r24  r34        r23
!   /  |   \
! C2---+-r3-C3        r12        r13  r14  r23  r24        r34  (first implem)
!   \  r6  /
! r5  \| / r2
!      H4
v = v_hh(r14)+v_cc(r23)+v_ch(r12)+v_ch(r13)+v_ch(r24)+v_ch(r34) + &
    v_cch(r12,r23,r13)+v_cch(r24,r23,r34)+ & ! both C to the same H
    v_chh(r12,r24,r14)+v_chh(r13,r34,r14)+ & ! both H to the same C
    v_cchh(r12,r34,r23,r13,r24,r14) + 17.62922453833992_dp

! convert eV to cm-1
v = v* ev2cm

! write(6,'(6f10.6)') r12,r34,r13,r14,r23,r24

CONTAINS


  FUNCTION v_ch(rch)
  USE types
  REAL(DP) :: v_ch, x, rch
  REAL(DP), DIMENSION(3) :: a = (/ 3.726_dp, 3.086_dp, 1.877_dp /)
  REAL(DP) :: de = 3.648_dp
  REAL(DP) :: Re = 1.120_dp

  x = rch-re
  v_ch = -de*(one+a(1)*x+a(2)*x*x+a(3)*x**3)*exp(-a(1)*x)

  END FUNCTION v_ch

  FUNCTION v_cc(rcc)
  USE types
  REAL(DP) :: v_cc, x, rcc
  REAL(DP), DIMENSION(3) :: a = (/ 5.089_dp, 6.898_dp, 4.059_dp /)
  REAL(DP) :: de = 6.271_dp
  REAL(DP) :: Re = 1.243_dp

  x = rcc-re
  v_cc = -de*(one+a(1)*x+a(2)*x*x+a(3)*x**3)*exp(-a(1)*x)

  END FUNCTION v_cc

  FUNCTION v_hh(rhh)
  USE types
  REAL(DP) :: v_hh, x, rhh
  REAL(DP), DIMENSION(3) :: a = (/ 3.880_dp, 3.741_dp, 3.255_dp /)
  REAL(DP) :: de = 4.750_dp
  REAL(DP) :: Re = 0.742_dp

  x = rhh-re
  v_hh = -de*(one+a(1)*x+a(2)*x*x+a(3)*x**3)*exp(-a(1)*x)

  END FUNCTION v_hh

  FUNCTION v_cch(rhc, rcc1, rhc1)
  USE types
  REAL(DP) :: v_cch, rhc, rcc1, rhc1
  REAL(DP) :: x1, x2, x3
  REAL(DP) :: a1 = -0.494_dp
  REAL(DP) :: a2 =  0.842_dp
!  REAL(DP) :: a3 = -0.494_dp
  REAL(DP) :: a11 = -0.781_dp
  REAL(DP) :: a22 = -2.799_dp
!  REAL(DP) :: a33 = -0.781_dp
  REAL(DP) :: a12 = -0.181_dp
  REAL(DP) :: a13 =  1.937_dp
!  REAL(DP) :: a23 = -0.181_dp
  REAL(DP) :: re = 1.511_dp
  REAL(DP) :: gamma1 = 2.0_dp
  REAL(DP) :: vi = 1.884_dp

  x1 = rhc -re
  x2 = rcc1-re
  x3 = rhc1-re
  v_cch = vi* &
  (one+a1*(x1+x3)+a2*x2+ &
   a11*(x1*x1+x3*x3)+a22*x2*x2+  &
   a12*(x1*x2+x2*x3)+a13*x1*x3)*  &
   (one-tanh(gamma1*(x1+x2+x3)*0.5_dp/sqrt3))

  END FUNCTION v_cch

  FUNCTION v_chh(rch, rch1, rhh1)
  USE types
  REAL(DP) :: v_chh, rch, rch1, rhh1
  REAL(DP) :: x1, x2, x3
  REAL(DP) :: a1 = -0.854_dp
!  REAL(DP) :: a2 = -0.854_dp
  REAL(DP) :: a3 =  3.964_dp
  REAL(DP) :: a11 = -20.613_dp
!  REAL(DP) :: a22 = -20.613_dp
  REAL(DP) :: a33 = -25.692_dp
  REAL(DP) :: a12 = -41.970_dp
  REAL(DP) :: a13 =  46.835_dp
!  REAL(DP) :: a23 =  46.835_dp
  REAL(DP) :: rech = 1.088_dp
  REAL(DP) :: rehh = 2.018_dp
  REAL(DP) :: gamma1 = 3.0_dp
  REAL(DP) :: gamma3 = 1.0_dp
  REAL(DP) :: vi = -0.424_dp

  x1 = rch -rech
  x2 = rch1-rech
  x3 = rhh1-rehh
  v_chh = vi* &
  (one+a1*(x1+x2)+a3*x3+ &
   a11*(x1*x1+x2*x2)+a33*x3*x3+  &
   a12*x1*x2+a13*(x1*x3+x2*x3))* &
   (one-tanh(gamma1*x1*0.5_dp))* &
   (one-tanh(gamma1*x2*0.5_dp))* &
   (one-tanh(gamma3*x3*0.5_dp))

  END FUNCTION v_chh

  FUNCTION v_cchh(rc1h1, rc2h2, rc1c2, rc2h1, rc1h2, rh1h2)
  USE types
  REAL(DP) :: v_cchh, rc1h1, rc2h2, rc1c2, rc2h1, rc1h2, rh1h2
!  REAL(DP) :: r1, r2, r3, r4, r5, r6
  REAL(DP) :: s1, s2, s3, s4, s5, s6
  REAL(DP) :: a1  =  0.11107_dp
  REAL(DP) :: a2  = -1.31967_dp
  REAL(DP) :: a3  = -1.75405_dp
  REAL(DP) :: a11 =   0.40750_dp
  REAL(DP) :: a22 =  12.16518_dp
  REAL(DP) :: a33 =  -0.55689_dp
  REAL(DP) :: a44 =  -0.15247_dp
  REAL(DP) :: a55 =   0.04666_dp
  REAL(DP) :: a66 =  -0.50762_dp
  REAL(DP) :: a12 =  -4.25660_dp
  REAL(DP) :: a13 =  -0.35050_dp
  REAL(DP) :: a23 =   3.32452_dp
  REAL(DP) :: a166 =  -0.01873_dp
  REAL(DP) :: a266 =  -0.08141_dp
  REAL(DP) :: a344 =   0.33917_dp
  REAL(DP) :: a355 =  -0.11011_dp
  REAL(DP) :: a366 =   0.09814_dp
  REAL(DP) :: a6666 =  0.05775_dp
  REAL(DP) :: rech = 1.5_dp
  REAL(DP) :: recc = 1.2033_dp
  REAL(DP) :: rehh = 2.7481_dp
  REAL(DP) :: gamma1 =  2.9_dp
  REAL(DP) :: gamma2 = -9.0_dp
  REAL(DP) :: gamma3 = -1.5_dp
  REAL(DP) :: vi = 4.41238_dp

!  r1 = rc1h1 -rech
!  r2 = rc2h2 -rech
!  r3 = rc1c2 -recc
!  r4 = rc2h1 -rech
!  r2 = rc1h2 -rech
!  r6 = rh1h2 -rehh

  s1 = rc1h1+rc2h2+rc2h1+rc1h2-4.0_dp*rech  ! r1+r2+r3+r4
  s2 = rc1c2 -recc                          ! r3
  s3 = rh1h2 -rehh                          ! r6
  s4 = rc1h1+rc2h2-rc2h1-rc1h2              ! r1+r2-r4-r5
  s5 = rc1h1-rc2h2+rc2h1-rc1h2              ! r1-r2+r4-r5
  s6 = rc1h1-rc2h2-rc2h1+rc1h2              ! r1-r2-r4+r5

  v_cchh = vi* &
  (one+a1*s1+a2*s2+a3*s3+ &
   a11*s1*s1+a22*s2*s2+a33*s3*s3+a44*s4*s4+a55*s5*s5+a66*s6*s6+ &
   a12*s1*s2+a13*s1*s3+a23*s2*s3+ &
   a166*s1*s6*s6+a266*s2*s6*s6+a344*s3*s4*s4+  &
   a355*s3*s5*s5+a366*s3*s6*s6+a6666*s6**6 )*  &
   (one-tanh(gamma1*s1*0.25_dp))* &
   (one-tanh(gamma2*s2*0.5_dp))*(one-tanh(gamma3*s3*0.5_dp))

  END FUNCTION v_cchh


  SUBROUTINE ac(a)
  USE types
  USE param
  IMPLICIT NONE
  REAL(DP) :: a(4,3)
  REAL(DP) :: m1,m2,m3,m4, bigm, m12,m34

  ! write(6,*) 'warning: ac has been executed once.'

  m1 = mass1
  m2 = mass2
  m3 = mass3
  m4 = mass4

  m12= m1+m2
  m34= m3+m4
  bigm= m12+m34

  a(1,1)=  m2/m12
  a(1,2)=  0.0_dp
  a(1,3)= -m34/bigm
  a(2,1)= -m1/m12
  a(2,2)=  0.0_dp
  a(2,3)= -m34/bigm
  a(3,1)=  0.0_dp
  a(3,2)= -m4/m34
  a(3,3)=  m12/bigm
  a(4,1)=  0.0_dp
  a(4,2)=  m3/m34
  a(4,3)=  m12/bigm

!  write(6,*) a(1,:)
!  write(6,*) a(2,:)
!  write(6,*) a(3,:)
!  write(6,*) a(4,:)

  END SUBROUTINE ac


END FUNCTION
