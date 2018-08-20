! eps has been changed to 1.0d-14
!
! this is modified Stroud & Secrest program in Fortran 90
! alf  can be real
!
! Notes: 
! 1. JT devided weights by gamma(nn+alf+1) (also see comments below)
!    I used more standard definition of Laguerre polynomial but still not
!    quite standard (minus sign) and devided weights by gamma(nn+alf+1)
! 2. the accuracy of weights as tested against the sum is worse when
!    computed using gammln than as a loop product at least for small n.
!    coordinates are fine (relative accuracy 1.E-16) but (tsa-csa)/tsa
!    was only about 1E-12, alf=51 with eps=1E-12. However when eps
!    increased so increased the accuracy (itmax must be increased also)
!    but still can not be decreased down to 1E-15 even at itmax=100.
!    Perhaps new algorithm to compute w_i is needed. As a result eps 
!    was set to 1E-14 (i.e. slightly smaller then original 1E-12)
! 3. At n=160 the initial guess about x1 becomes bad

MODULE laguer_grid

INTERFACE laguer
  SUBROUTINE laguer(nn,x,a,alf,csx,csa,tsx,tsa)
  USE types
  USE base_lib, ONLY: lagpoly, gammln
  INTEGER(I4B), INTENT(IN) :: nn
  REAL(DP), INTENT(IN) :: alf
  REAL(DP), INTENT(OUT) :: x(:),a(:), csx,csa,tsx,tsa
  END SUBROUTINE laguer
END INTERFACE

INTERFACE lgroot
  SUBROUTINE lgroot(x,nn,alf,dpn,pn1)
  USE types
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: nn
  REAL(DP), INTENT(IN) :: alf
  REAL(DP), INTENT(OUT) :: x,dpn,pn1
  END SUBROUTINE lgroot
END INTERFACE

INTERFACE lgrecr
  SUBROUTINE lgrecr(pn,dpn,pn1,x,nn,alf)
  USE types
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: nn
  REAL(DP), INTENT(IN) :: x, alf
  REAL(DP), INTENT(OUT) :: pn,dpn,pn1
  END SUBROUTINE lgrecr
END INTERFACE

END MODULE


SUBROUTINE laguer(nn,x,a,alf,csx,csa,tsx,tsa)
!
!     calculates points & weights for gauss-laguerre integration see:
!     "Gaussian quadrature formulas" by A.H.Stroud & D.Secrest
!      1966, Prentice-Hall, p.32.
!
!
! for all N less or equal to the highest degree nn.
!
!   csx = calc sum x(i)     tsx = true sum x(i)
!   csa = calc sum a(i)     tsa = true sum a(i)
!
USE types
USE base_lib, ONLY: lagpoly, gammln
USE laguer_grid, ONLY: lgroot, lgrecr
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: nn
REAL(DP), INTENT(IN) :: alf
REAL(DP), INTENT(OUT) :: x(:),a(:), csx,csa,tsx,tsa
REAL(DP) :: fn, r1,r2,ratio, xt,dpn,pn1,fi
INTEGER(I4B) :: i,j
REAL(DP) :: step,xt1,xt2,pt,pt2

fn = dble(nn)
csx = 0.0d0
csa = 0.0d0
tsx = fn*(fn+alf)
tsa=1.d0
do i=2,nn
  tsa=tsa*i/(alf+i-1.d0)
end do

xt=(1.0d0+alf)*(2.0d0+alf)/(1.0d0+3.0d0*fn+2.0d0*alf)
step=3.0d0*(1.0d0+alf)/(1.0d0+3.0d0*fn+alf)
call lgrecr(pt,dpn,pn1,xt,nn,alf)

do i=1,nn

!**************************************************************
! J.T.: formulas for initial point & step chosen because they work!
!
  if ( i==1.or.i==2 ) then
    ! smallest two zeros: found by "brute force" search
    iter_loop: do j=1,100
      xt2 = xt + step
      call lgrecr(pt2,dpn,pn1,xt2,nn,alf)
      if (sign(1.d0,pt)*sign(1.d0,pt2) > 0.0d0) then
        pt = pt2
        xt = xt2
      else
        pt = pt2
        xt = 0.5d0 * (xt + xt2)
        exit
      end if
    end do iter_loop

    if (j>=100) then
      write(u6,*) ' laguer: iteration number reached the limit'
      stop
    end if

  else
    !  all other zeros
    fi = real((i-2),dp)
    r1 = (1.0d0+2.55d0*fi)/(1.9d0*fi)
    r2 = 1.26d0*fi*alf/(1.0d0+3.5d0*fi)
    ratio = (r1+r2)/(1.0d0+0.3d0*alf)
    xt = xt + ratio*(xt-xt2)
  end if

  call lgroot(xt,nn,alf,dpn,pn1)

  xt2=xt1
  xt1=xt
  x(i) = xt
  a(i) = 1.d0/dpn/pn1
  csx = csx + xt
  csa = csa + a(i)
  ! write(u6, '(2x,i4,f26.12,e26.12)') i, x(i), a(i)

end do

END SUBROUTINE laguer


SUBROUTINE lgroot(x,nn,alf,dpn,pn1)
!
!     improves the approximate root x; in addition obtains
!          dpn = derivative of p(n) at x
!          pn1 = value of p(n-1) at x
!     this routine is due to Stroud & Secrest (see subroutine laguer)
!
USE types
USE laguer_grid, ONLY: lgrecr
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: nn
REAL(DP), INTENT(IN) :: alf
REAL(DP), INTENT(INOUT) :: x
REAL(DP), INTENT(OUT) :: dpn,pn1
REAL(DP) :: eps,dpp,d,p
INTEGER(I4B) :: iter, itmax

!*** set accuracy of finding the roots:
! eps = 1.0d-12
! itmax = 10
eps = 5.0d-15
! eps = 1.0d-14
itmax = 5000

do iter=1,itmax
  call lgrecr(p,dpp,pn1,x,nn,alf)
  d = p/dpp
  x = x-d
  if ((dabs(d/x)-eps)<0) exit
end do
if ( (dabs(d/x)-eps)>0 ) then
  write(u6,*) ' insufficient iter=itmax, delta= ', dabs(d/x)
  stop
end if
dpn = dpp

END SUBROUTINE lgroot


SUBROUTINE lgrecr(pn,dpn,pn1,x,nn,alf)
!
!  uses recurrence relations to set up polynomials
!  this routine is due to stroud & secrest (see subroutine laguer)
!  modifications have been made to use standard Laguerre polynomials
!  i.e. the original recursion relation 
!     p(n) = (x-b(n))*p(n-1) - c(n)*p(n-2)
!  has been replaced 
!     p(n) = [(x-B(n))*p(n-1) - C(n)*p(n-2)]/n
!  where
!     B(n) = (alf+2*n-1)
!     C(n) = (alf+n-1)
!  it is still not standard since there should be -(x-B(n))
!     L(n) [standard] = (-1)^n /n!*L(n) [Stroud] 
!  Standard recurrence relation
!     p(n) = {[(2*n+a-1)-x]*p(n-1)-(n+a-1)*p(n-2)}/n
!  and the STANDARD derivative is
!     p'(n) = {[(2*n+a-1)-x]*p'(n-1)-(n+a-1)*p'(n-2)-p(n-1)}/n
!  but Stroud has
!     p'(n) = {[x-(2*n+a-1)]*p'(n-1)-(n+a-1)*p'(n-2)+p(n-1)}/n
!
!  have not changed this since that might have something to do
!  with how derivatives behave near zero and I don't want to study
!  this. So the corrections are: c changed and P and P' devided by n
!
!  no need anymore to prepare in advance the coefficients
!     b(n) = (alf+2*n-1)
!     c(n) = (n-1)*(alf+n-1)
!
!  note: 
!  since lgroot uses only the ratio p(n)/p'(n)  no changes are required
!
USE types
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: nn
REAL(DP), INTENT(IN) :: x, alf
REAL(DP), INTENT(OUT) :: pn,dpn,pn1
REAL(DP) :: p,p1,dpp,dp1,q,dq, b,c
INTEGER(I4B) :: j
!  initiate:
p1 = 1.0d0             ! p(n=0)
p = x - alf - 1.0d0    ! p(n=1); standard = alf+1+x
dp1 = 0.0d0            ! p'(n=0)
dpp = 1.0d0            ! p'(n=1)
do j=2,nn
  b = alf+2.d0*real(j,dp)-1.d0
  c = alf+j-1.d0
  q  = ((x-b)*p-c*p1)/real(j,dp)       ! p(j)
  dq = ((x-b)*dpp-c*dp1+p)/real(j,dp)  ! p'(j)
  p1 = p                               ! p(j-2) for the next loop
  p  = q                               ! p(j-1) for the next loop
  dp1= dpp                             ! p'(j-2) for the next loop
  dpp = dq                             ! p'(j-1) for the next loop
end do
pn = p                 ! p(nn)
dpn= dpp               ! p'(nn)
pn1= p1                ! p(nn-1)

END SUBROUTINE lgrecr
