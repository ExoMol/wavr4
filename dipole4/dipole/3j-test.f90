! 3j-test  - test 3J formulae and 3j.f90 code against drc3jj.f code
!
! how to compile:
! ifort -c -fixed drc3jj.f
! ifort -o 3j-test.exe  3j-test.f90 3j.f90 types.f90 drc3jj.o
! ./3j-test.exe 3

PROGRAM threejtest
  IMPLICIT NONE
  INTEGER(4) :: j1,m1,j,m, narg, jmax, js, je
  INTEGER(4) :: ier
  REAL(8) :: my3j, tmp, ref, l1min, l1max
  REAL(8) :: thrcof(10)
  CHARACTER*10 :: arg

! input processing
narg=iargc()

call getarg(1,arg)
read(arg,'(i10)') jmax

! get input
!if (narg > 3) then
!  call getarg(1,arg)
!  read(arg,'(i10)') j2
!  call getarg(2,arg)
!  read(arg,'(i10)') m2
!  call getarg(3,arg)
!  read(arg,'(i10)') j
!  call getarg(4,arg)
!  read(arg,'(i10)') m
!end if
!  end of user input

write(6,*) "(j1,m1, 1,m2, j, m)"

do j1=0,jmax
do m1=-j1,j1

js = max(j1-1,0)
je = max(j1+1,0)
do j= js, je
do m=-j,j

if (abs(m1+m) > 1) cycle

tmp = my3j(j1,m1,j,m)

! (j1,m1, 1,m2, j, m) = (-1)^(j1+1+j) ( 1,m2,j1,m1, j, m)

!     DRC3JJ (L2, L3, M2, M3, L1MIN, L1MAX, THRCOF, NDIM, IER)
call  DRC3JJ (dble(j1),dble(j), dble(m1),dble(m), l1min, l1max, THRCOF,   10, IER)

ref = thrcof(2-int(l1min))
if (mod(j1+1+j,2) > 0 )  ref = -ref

write(6,'(" (",i2,",",i2,", 1,",i2,",",i2,",",i2,")  ",3f20.16)') j1,m1,(-m1-m),j,m, tmp, ref, (ref-tmp)


enddo
enddo
enddo
enddo


END PROGRAM
