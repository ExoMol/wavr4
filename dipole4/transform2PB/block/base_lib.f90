! auxiliary programms library mainly based on Numerical Recipies

MODULE base_lib

INTERFACE print_stat
  SUBROUTINE print_stat(id)
  USE types
! IMPLICIT NONE
  INTEGER(I4B) :: id
  END SUBROUTINE
END INTERFACE

INTERFACE time_stamp
  SUBROUTINE time_stamp(output_unit)
  USE types
! IMPLICIT NONE
  INTEGER(I4B), OPTIONAL :: output_unit
  END SUBROUTINE
END INTERFACE

INTERFACE arth
  FUNCTION arth_d(first,increment,n)
  USE types
! IMPLICIT NONE
  REAL(DP), INTENT(IN) :: first,increment
  INTEGER(I4B), INTENT(IN) :: n
  REAL(DP), DIMENSION(n) :: arth_d
  END FUNCTION arth_d

  FUNCTION arth_i(first,increment,n)
  USE types
! IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: first,increment,n
  INTEGER(I4B), DIMENSION(n) :: arth_i
  END FUNCTION arth_i
END INTERFACE

INTERFACE lagpoly
  FUNCTION lagpoly(n,alf,x)
  USE types
! IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: n
  REAL(DP), INTENT(IN) :: x, alf
  REAL(DP) :: lagpoly
  END FUNCTION lagpoly
END INTERFACE

INTERFACE gammln
  FUNCTION gammln_d(xx)
  USE types
! IMPLICIT NONE
  REAL(DP),INTENT(IN)::xx
  REAL(DP)::gammln_d
  END FUNCTION gammln_d

  FUNCTION gammln_i(xx)
  USE types
! IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: xx
  REAL(DP)::gammln_i
  END FUNCTION gammln_i
END INTERFACE

INTERFACE gauleg
  SUBROUTINE gauleg_d(x1,x2,x,w)
  USE types
! IMPLICIT NONE
  REAL(DP),INTENT(IN)::x1,x2
  REAL(DP),DIMENSION(:),INTENT(OUT)::x,w
  END SUBROUTINE gauleg_d
END INTERFACE

INTERFACE plgndr
  FUNCTION plgndr_d(l,m,x)
  USE types
! IMPLICIT NONE
  INTEGER(I4B),INTENT(IN)::l,m
  REAL(DP),INTENT(IN)::x
  REAL(DP)::plgndr_d
  END FUNCTION plgndr_d
END INTERFACE

INTERFACE xn
  FUNCTION xn(i,j)
  USE types
! IMPLICIT NONE
  INTEGER(I4B) :: xn,i,j
  END FUNCTION xn
END INTERFACE

END MODULE


FUNCTION lagpoly(n,alf,x)
  USE types
IMPLICIT NONE
  REAL(DP) :: lagpoly
  REAL(DP), INTENT(IN) :: x, alf
  INTEGER(I4B), INTENT(IN) :: n
! Computes the normalized generalized Laguerre polynomial L^alf_n (x).
INTEGER(I4B) :: j
REAL(DP) :: p1,p2,p3
! recurrence relation to get the Laguerre polynomials evaluated at x
if (n<0) then
  p1=1.0_dp
  write(u6, *) ' lagpoly: n is negative!'
  stop
else if (n==0) then
  p1=1.0_dp
else
  p1=1.0_dp
  p2=0.0_dp
  do j=1,n
    p3=p2
    p2=p1
    p1=((2.0_dp*j-1.0_dp+alf-x)*p2-(j-1.0_dp+alf)*p3)/j
  end do
end if

lagpoly=p1

! devide by the norm factor
! lagpoly_d = p1 /sqrt(exp(gammln(alf+n+1)-gammln(real(n+1,dp))))

END FUNCTION lagpoly

! Num Rec routine
! changed to double prec
! NOTE that this is actually not LOG((xx)!) but LOG(GAMMA(xx)); N!=GAMMA(N+1)
!
FUNCTION gammln_d(xx)
  USE types
  USE base_lib, ONLY : arth
IMPLICIT NONE
  REAL(DP),INTENT(IN) :: xx
  REAL(DP) :: gammln_d
! Returns the value ln[Gamma(xx)] for xx > 0.
REAL(DP) :: tmp,x
! Internal arithmetic will be done in double precision, a nicety that you can
! omit if five-figure accuracy is good enough.
REAL(DP) :: stp =2.5066282746310005_dp
REAL(DP), DIMENSION(6) :: coef =(/76.18009172947146_dp,&
  -86.50532032941677_dp,24.01409824083091_dp,-1.231739572450155_dp,&
  0.1208650973866179e-2_dp,-0.5395239384953e-5_dp/)
if (xx<=0.0) then
  write(u6,*) ' negative or 0 argument in gammln_d', xx
  stop
end if
x=xx
tmp=x+5.5_dp
tmp=(x+0.5_dp)*log(tmp)-tmp
gammln_d=tmp+log(stp*(1.000000000190015_dp+&
  sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
END FUNCTION gammln_d

FUNCTION gammln_i(n)
  USE types
IMPLICIT NONE
  INTEGER(I4B),INTENT(IN) :: n
  REAL(DP) :: gammln_i
! Returns the value ln[Gamma(n)] for n >= 0.
INTEGER(I4B) :: i
if (n<1) then
  write(u6,*) ' negative argument in gammln_i', n
  stop
else if (n<3) then
  gammln_i=0.0_dp
else
  gammln_i=0.0_dp
  do i=2,n-1
    gammln_i=gammln_i+log(real(i,dp))
  end do
endif
END FUNCTION gammln_i

! Num Rec routine
! changed to double prec
!
SUBROUTINE gauleg_d(x1,x2,x,w)
  USE types
  USE base_lib, ONLY : arth
IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x1,x2
  REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
REAL(DP), PARAMETER :: EPS=1.0e-14_dp
! Given the lower and upper limits of integration x1 and x2, this routine returns
! arrays x and w of length N containing the abscissas and weights of the
! Gauss-Legendre N-point quadrature formula. The parameter EPS is the relative
! precision. Note that internal computations are done in double recision.
INTEGER(I4B) :: its,j,m,n
INTEGER(I4B), PARAMETER :: MAXIT=10
REAL(DP) :: xl,xm
REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished

if (size(x) == size(w)) then
  n=size(x)
else
  write (u6,*)' gauleg_d: array sizes are different'
  STOP
end if

m=(n+1)/2                                   ! The roots are symmetric in the interval,
xm=0.5_dp*(x2+x1)                           !    so we only have to find half of them.
xl=0.5_dp*(x2-x1)
z=cos(PI*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))   ! Initial approximations to the roots.
unfinished=.true.
do its=1,MAXIT                              ! Newton's method carried out simultane-
  where (unfinished)                        !    ously on the roots.
    p1=1.0_dp
    p2=0.0_dp
  end where
  do j=1,n                                  ! Loop up the recurrence relation to get
    where (unfinished)                      !   the Legendre polynomials evaluated at z.
      p3=p2
      p2=p1
      p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
    end where
  end do
  ! p1 now contains the desired Legendre polynomials. We next compute pp,
  ! the derivatives, by a standard relation involving also p2, the polynomials 
  ! of one lower order.
  where (unfinished)
    pp=n*(z*p1-p2)/(z*z-1.0_dp)
    z1=z
    z=z1-p1/pp                              ! Newton s method.
    unfinished=(abs(z-z1)>EPS)
  end where
  if (.not. any(unfinished)) exit
end do
if (its ==MAXIT+1) then
  write(u6,*) 'gauleg_d: too many iterations in gauleg'
  STOP
end if
x(1:m)=xm-xl*z                              ! Scale the root to the desired interval,
x(n:n-m+1:-1)=xm+xl*z                       ! and put in its symmetric counterpart.
w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)      ! Compute the weight
w(n:n-m+1:-1)=w(1:m)                        ! and its symmetric counterpart.
END SUBROUTINE gauleg_d

! Num Rec routine
! changed to double prec
!
FUNCTION plgndr_d(l,m,x)
  USE types
  USE base_lib, ONLY : arth
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: l,m
  REAL(DP), INTENT(IN) :: x
REAL(DP) :: plgndr_d
! Computes the associated Legendre polynomial P^m_l (x). Here m and l are
! integers satisfying 0=< m =< l, while x lies in the range -1=< x =< 1
INTEGER(I4B) :: ll
REAL(DP) :: pll,pmm,pmmp1,somx2
! call assert(m >=0,m <=l,abs(x)<=1.0,'plgndr_d args ')
if (m < 0 .or. m > l .or. abs(x) > 1.0) then
  write(u6,*) ' check arguments in PLGNDR ', l, m, x
  stop
end if
pmm=1.0_dp                        ! Compute P^m_m
if (m > 0) then
  somx2=sqrt((1.0_dp-x)*(1.0_dp+x))
  pmm=product(arth(1.0_dp,2.0_dp,m))*somx2**m
  ! if (mod(m,2)==1) pmm= -pmm      ! I.K: we don't need sign alternation here
end if
if (l == m) then
  plgndr_d= pmm
else
  pmmp1=x*(2*m+1)*pmm             ! Compute P^m_(m+1)
  if (l ==m+1)then
    plgndr_d=pmmp1
  else                            ! Compute P^m_l, l > m+1
  do ll=m+2,l
    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
    pmm=pmmp1
    pmmp1=pll
  end do
  plgndr_d=pll
  end if
end if
END FUNCTION plgndr_d

! =============================================================================
! Returns an array of length n containing an arithmetic progression whose first
! value is "first" and whose increment is "increment". If "first" and "increment"
! have rank greater than zero, returns an array of one larger rank, with the
! last subscript having size "n" and indexing the progressions. Note that
! the following reference implementation (for the scalar case) is definitional
! only, and neither parallelized nor optimized for roundoff error. See 22.2
! and Appendix C1 for implementation by subvector scaling.


FUNCTION arth_d(first,increment,n)
  USE types
IMPLICIT NONE
  REAL(DP), INTENT(IN) :: first,increment
  INTEGER(I4B), INTENT(IN) :: n
  REAL(DP), DIMENSION(n) :: arth_d
INTEGER(I4B), PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8  ! Each NPAR2 must be =< then
                                                       ! corresponding NPAR.
INTEGER(I4B) :: k,k2
REAL(DP) :: temp
if (n >0)arth_d(1)=first
if (n <=NPAR_ARTH)then
  do k=2,n
    arth_d(k)=arth_d(k-1)+increment
  end do
else
  do k=2,NPAR2_ARTH
    arth_d(k)=arth_d(k-1)+increment
  end do
  temp=increment*NPAR2_ARTH
  k=NPAR2_ARTH
  do
    if (k >=n)exit
    k2=k+k
    arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
    temp=temp+temp
    k=k2
  end do
end if
END FUNCTION arth_d


FUNCTION arth_i(first,increment,n)
  USE types
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: first,increment,n
  INTEGER(I4B), DIMENSION(n) :: arth_i
INTEGER(I4B), PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8  ! Each NPAR2 must be =< then
                                                       ! corresponding NPAR.
INTEGER(I4B) :: k,k2,temp
if (n >0)arth_i(1)=first
if (n <=NPAR_ARTH) then
  do k=2,n
    arth_i(k)=arth_i(k-1)+increment
  end do
else
  do k=2,NPAR2_ARTH
    arth_i(k)=arth_i(k-1)+increment
  end do
  temp=increment*NPAR2_ARTH
  k=NPAR2_ARTH
  do
    if (k >=n) exit
    k2=k+k
    arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
    temp=temp+temp
    k=k2
  end do
end if
END FUNCTION arth_i

! include 'timer_unix.inc'
! include 'timer_pc.inc'

!  unix timer, it works as a dummy with MSF portlib.lib
SUBROUTINE time_stamp(output_unit)
  !EXTERNAL dtime,etime
  USE types
IMPLICIT NONE
  INTEGER(I4B), OPTIONAL :: output_unit
REAL(SP) :: tarray(2),etime,dtime

  if(.NOT.present(output_unit)) then
    output_unit = 6
  end if
  write(output_unit,'(" Last step:",F10.2,"s; total:",F10.2,"s.")')  &
    dtime(tarray), etime(tarray)

END SUBROUTINE


SUBROUTINE print_stat(id)
  USE types
IMPLICIT NONE
  INTEGER(I4B) :: id
  if(id /= 0) then
    write(u6,'("  allocation problem:",i5)') id
    stop
  end if
END SUBROUTINE

!  function xn(i,j) maps double index i,j to single index  xn
!  assuming UPPER triangle matrix
!
FUNCTION xn(i,j)
  USE types
IMPLICIT NONE
INTEGER(I4B) :: xn,i,j
if(1<=i .AND. i<=j) then
  xn = i+(j-1)*j/2
else if(1<=j .AND. j<i) then  ! swap i and j
  xn = j+(i-1)*i/2
else
  write(*,*) '  xn: wrong indices:', i,j
end if
END FUNCTION
