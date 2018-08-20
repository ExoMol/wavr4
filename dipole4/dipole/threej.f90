! Simple code to compute 3J symbols in the form
! [ J'    1    J]
! [ M' (-M'-M) M]
! i.e. it's only suitable for dipole calculations.
!
! INPUT: my3j( J', M', J, M)
!
! Probably need to make sure that the formulae are valid for M < 0 too.
!
! 3J symbol properties:
!
! [j1 j2 j3]   [j2 j3 j1]   [j3 j1 j2]
! [m1 m2 m3] = [m2 m3 m1] = [m3 m1 m2]
!
! [j1 j2 j3]       (j1+j2+j3) [j2 j1 j3]
! [m1 m2 m3] = (-1)           [m2 m1 m3]
!
! [ j1  j2  j3]       (j1+j2+j3) [j1 j2 j3]
! [-m1 -m2 -m3] = (-1)           [m1 m2 m3]
!
! m1+m2+m3=0  => m2= -m1-m3
!
! The firmulae below only work for J' = J or J' = J+1
! Hence we need to swap J' and J if J'<J
!
! [j1 j2 j3]       (j1+j2+j3) [j3 j2 j1]
! [m1 m2 m3] = (-1)           [m3 m2 m1]
!

MODULE threej

INTERFACE my3j
  FUNCTION my3j(jj1,mm1,jj,mm)
  USE types
  INTEGER(I4B), INTENT(IN) :: jj1,mm1,jj,mm
  REAL(DP) :: my3j
  END FUNCTION my3j
END INTERFACE

END MODULE threej


FUNCTION my3j(jj1,mm1,jj,mm)
  USE types
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: jj1,mm1,jj,mm
  INTEGER(I4B) :: j1,m1,m2,j,m  ! working copies of input
  REAL(DP) :: my3j, factor

factor = 1.D0

if (abs(mm1+mm) > 1) then
  write(6,*) ' 3J: delta M > 1'
  stop
elseif (abs(jj1-jj) > 1) then
  write(6,*) ' 3J: delta J > 1'
  stop
elseif (jj1 < 0 .OR. jj < 0) then
  write(6,*) ' 3J: negative J'
  stop 
elseif (abs(mm1) > jj1 .OR. abs(mm) > jj) then
  write(6,*) ' 3J: M > J'
  stop 
endif

  j1=jj1
  m1=mm1
  j=jj
  m=mm

! swap (j1,m1) and (j,m) if j1 < j
if (j1 < j) then
  j1=jj
  m1=mm
  j=jj1
  m=mm1
  if(mod(j1+1+j,2) /=0)  factor = -factor
endif

m2 = -m1-m

if (j1==0 .AND. j==0) then
  my3j = 0.D0

elseif (j1==j) then
  if (m2==0) then
    my3j = m/(sqrt(j*(j+1.D0)*(2*j+1.D0)))
    if(mod(j-m,2) /=0) factor= -factor
  elseif (m2==  1) then
    my3j = sqrt((j-m)*(j+m+1.D0)/(2.D0*j*(j+1.D0)*(2*j+1.D0)))
    if(mod(j-m,2) /=0) factor= -factor
  elseif (m2== -1) then
    my3j = sqrt((j+m)*(j-m+1.D0)/(2.D0*j*(j+1.D0)*(2*j+1.D0)))
    if(mod(j-m+1,2) /=0) factor= -factor
  endif
else ! j1 > j
  if (m2==0) then
    my3j = sqrt(2.D0*(j+m+1.D0)*(j-m+1.D0)/((2.D0*j+1.D0)*(2.D0*j+2.D0)*(2.D0*j+3.D0)))
    if(mod(j-m+1,2) /=0) factor= -factor
  elseif (m2==  1) then
    my3j = sqrt((j+m+1.D0)*(j+m+2.D0)/((2.D0*j+1.D0)*(2.D0*j+2.D0)*(2.D0*j+3.D0)))
    if(mod(j-m,2) /=0) factor= -factor
  elseif (m2== -1) then
    my3j = sqrt((j-m+1.D0)*(j-m+2.D0)/((2.D0*j+1.D0)*(2.D0*j+2.D0)*(2.D0*j+3.D0)))
    if(mod(j-m,2) /=0) factor= -factor
  endif
endif

my3j = factor* my3j

! write(6,*) jj1,mm1,jj,mm,my3j

END FUNCTION
