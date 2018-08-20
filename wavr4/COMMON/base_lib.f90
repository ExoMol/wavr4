! auxiliary programs

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

INTERFACE diag
  SUBROUTINE diag(ham,eng,nlev,nout,icut,encut,istage,idiag)
  USE types
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEV, LA_SYEVD, LA_SYEVX, LA_SYEVR
! IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)  :: nlev
  INTEGER(I4B), INTENT(OUT) :: nout,idiag
  INTEGER(I4B), INTENT(IN)  :: icut,istage
  REAL(DP), INTENT(INOUT) :: ham(nlev,nlev)
  REAL(DP), INTENT(OUT)   :: eng(nlev)
  REAL(DP), INTENT(IN)    :: encut
  END SUBROUTINE
END INTERFACE

INTERFACE diag_p
SUBROUTINE diag_p(ham,eng,vectors,nlev,nout,icut,encut,istage,idiag)
  USE types
  USE LA_PRECISION, ONLY: WP => DP
!  USE my_lapack95
  USE F95_LAPACK, ONLY: LA_SPEVX
! IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)  :: nlev
  INTEGER(I4B), INTENT(OUT) :: nout,idiag
  INTEGER(I4B), INTENT(IN)  :: icut,istage
  REAL(DP), INTENT(INOUT) :: ham((nlev*(nlev+1))/2)
  REAL(DP), INTENT(OUT)   :: eng(nlev)
  REAL(DP), INTENT(OUT)   :: vectors(nlev,icut)
  REAL(DP), INTENT(IN)    :: encut
  END SUBROUTINE
END INTERFACE


END MODULE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!================================================================
!
! LOG( GAMMA (n) ) for integer arguments
!
! NOTE that this is actually not LOG((xx)!) but LOG(GAMMA(xx)); N!=GAMMA(N+1)

FUNCTION gammln_i(n)
  USE types
IMPLICIT NONE
  INTEGER(I4B),INTENT(IN) :: n
  REAL(DP) :: gammln_i
! Computes value of ln[Gamma(n)] for n >= 0.
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

!================================================================
!
! LOG( GAMMA (n) ) for real*8 arguments

FUNCTION gammln_d(XX)
USE types
IMPLICIT NONE
REAL(DP)::gammln_d,gamma_log,XX
gammln_d= gamma_log(XX)
RETURN
END FUNCTION gammln_d

!================================================================

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

! ===========================================================================
! Gauss-Legendre N-point quadrature (call Jacobi quandrature)
!
! Originally using lower and upper limits of integration x1 and x2
! which are presently ignored assuming -1 .. 1 interval
! arrays x(1:N), w(1:N) return the abscissas and weights respectively

SUBROUTINE gauleg_d(x1,x2,x,w)
  USE types
IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x1,x2
  REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
INTEGER(I4B) :: n, nn, nb, i, j, itmp
REAL(DP),ALLOCATABLE :: x_temp(:),w_temp(:),basis(:,:)
REAL(DP) :: tmp

if (size(x) == size(w)) then
  n=size(x)
else
  write (u6,*)' gauleg_d: array sizes are different'
  STOP
end if

nn=n-1
nb=nn
ALLOCATE(x_temp(0:nn),w_temp(0:nn),basis(0:nb,0:nn))

CALL jac_basis(nn,nb,0.0D0,0.0D0,x_temp,w_temp,basis)

x(1:n)=x_temp(0:nn)
w(1:n)=w_temp(0:nn)

DEALLOCATE(x_temp,w_temp,basis)

! sort the grid in increasing order

do i=1,n
  itmp=i
  do j=i+1,n
    if(x(itmp) > x(j)) itmp = j
  end do
  ! swap
  tmp=x(itmp)
  x(itmp)=x(i)
  x(i)=tmp

  tmp=w(itmp)
  w(itmp)=w(i)
  w(i)=tmp

end do

END SUBROUTINE gauleg_d

! ===========================================================================
!
! subroutines and functions associated with computation of
! Jacobi quadrature. Used here to compute Gauss-Legendre
! quadrature by setting alpha=beta=0

! PROGRAM jacobi
! USE types
! IMPLICIT NONE
! INTEGER :: nb, nn, alpha, beta
! REAL(DP) :: alf, bet
! REAL(DP),ALLOCATABLE :: x(:),w(:),basis(:,:)
! print *, ' enter: alpha, beta, Nmax, Npoints'
! read *, alpha, beta, nn, nb
!   nb=40
!   nn=40
! alf=1.0D0*alpha
! bet=1.0D0*beta
! ALLOCATE(x(0:nn),w(0:nn),basis(0:nb,0:nn))
! CALL jac_basis(nn,nb,alf,bet,x,w,basis)
! DEALLOCATE(x,w,basis)
! END PROGRAM jacobi

SUBROUTINE jac_basis(nn,nb,alf,bet,x,w,basis)
USE types
IMPLICIT NONE
INTEGER :: n,i,nb,nn
REAL(DP) :: alf, bet,sum
REAL(DP) :: x(0:nn),w(0:nn),basis(0:nb,0:nn)
REAL(DP), ALLOCATABLE :: bass(:,:),norm(:) 
ALLOCATE(bass(0:nb,0:nn),norm(0:nb))
CALL norms(norm,nn,alf,bet)
CALL jacgt(x,w,bass,alf,bet,nn,nb)
DO i=0,nb
   DO n=0,nn
      basis(i,n)=bass(i,n)*norm(i)
   END DO
END DO
! do i=0,nb/2
!    print '(4(2x,f18.16))', x(i), x(nb-i), w(i), (x(i)+x(nb-i))
! end do
DO i=0,nb
      sum=0.0
      DO n=0,nn
         sum=sum+w(n)*basis(i,n)*basis(i,n)
      END DO
!  print '(" norm= ",f18.16)',sum
END DO
DEALLOCATE(bass,norm)
END SUBROUTINE jac_basis

SUBROUTINE jacgt(x,w,bass,alf,bet,nn,nb)
USE types
USE base_lib, ONLY: gammln
IMPLICIT NONE
INTEGER :: n,I,nn,nb
! REAL(DP),EXTERNAL :: gammln
REAL(DP) :: alf, bet,lmd, x1,x2
REAL(DP) :: x(0:nn),w(0:nn),bass(0:nb,0:nn)
REAL(DP), ALLOCATABLE :: A1n(:),A2n(:),A3n(:),A4n(:)
DATA x1,x2/1.0d0,2.0d0/
lmd=alf+bet+x1
ALLOCATE(A1n(nb),A2n(nb),A3n(nb),A4n(nb))
DO n=1,nb
   A1n(n)=x2*(n+x1)*(n+lmd)*(x2*n+lmd-x1)
   A2n(n)=(x2*n+lmd)*(alf*alf-bet*bet)
   A3n(n)=(x2*n+lmd-x1)*(x2*n+lmd)*(x2*n+lmd+x1)
   A4n(n)=x2*(n+alf)*(n+bet)*(x2*n+lmd+x1)
END DO   
bass=0.0D0
DO I=0,nn
  CALL  gaujac(x,w,nn+1,alf,bet)
  !  print '(" i= ",i3,2(2x,D22.16))',i,x(i),w(i)
  bass(0,I)=x1
  IF(nb.LT.1) THEN
  ELSE
    bass(1,I)=(alf-bet+(lmd+x1)*x(I))/x2
    DO n=2,nb
       bass(n,I)=((A2n(n-1)+A3n(n-1)*x(I))*bass(n-1,I)&
          -A4n(n-1)*bass(n-2,I))/A1n(n-1)
    END DO
  END IF
END DO
DEALLOCATE(A1n,A2n,A3n,A4n)
END SUBROUTINE jacgt

SUBROUTINE gaujac(x,w,n,alf,bet)
USE types
USE base_lib, ONLY: gammln
IMPLICIT NONE
INTEGER :: n
REAL(DP):: alf,bet,x1,x2,x3
REAL(DP):: w(n),x(n)
REAL(DP),PARAMETER :: EPS=1.0D-14
INTEGER,PARAMETER :: MAXIT=10
INTEGER :: i,its,j
REAL(DP)::alfbet,an,bn,r1,r2,r3
REAL(DP)::c1,c2,c3,p1,p2,p3,pp,temp,z,z1
! REAL(DP),EXTERNAL :: gammln
DATA x1,x2,x3/1.0d0,2.0d0,3.0d0/
DO i=1,n
  IF(i==1)THEN
     an=alf/DBLE(n)
     bn=bet/DBLE(n)
     r1=(x1+alf)*(2.78D0/(4.0D0+DBLE(n*n))+0.768D0*an/DBLE(n))
     r2=x1+1.48D0*an+0.96D0*bn+0.452D0*an*an+0.83D0*an*bn
     z =x1-r1/r2
  ELSE IF(i==2)THEN
     r1=(4.1D0+alf)/((x1+alf)*(x1+0.156D0*alf))
     r2=x1+0.06D0*(DBLE(n)-8.0D0)*(1.0D0+0.12D0*alf)/DBLE(n)
     r3=x1+0.012*bet*(x1+0.25D0*DABS(alf))/DBLE(n)
     z=z-(x1-z)*r1*r2*r3
  ELSE IF(i==3)THEN
     r1=(1.67D0+0.28D0*alf)/(x1+0.37D0*alf)
     r2=x1+0.22D0*(DBLE(n)-8.0D0)/DBLE(n)
     r3=x1+8.0D0*bet/((6.28D0+bet)*DBLE(n*n))
     z=z-(x(1)-z)*r1*r2*r3
  ELSE IF(i==n-1)THEN
     r1=(x1+0.235D0*bet)/(0.766D0+0.119D0*bet)
     r2=x1/(x1+0.639D0*(DBLE(n)-4.0D0)/&
         (x1+0.71D0*(DBLE(n)-4.0D0)))
     r3=x1/(x1+20.0D0*alf/((7.5D0+alf)*DBLE(n*n)))
     z=z+(z-x(n-3))*r1*r2*r3
  ELSE IF(i==n)THEN
     r1=(x1+0.37D0*bet)/(1.67D0+0.28D0*bet)
     r2=x1/(x1+0.22D0*DBLE(n-8)/DBLE(n))
     r3=x1/(x1+8.0D0*alf/((6.28D0+alf)*DBLE(n*n)))
     z=z+(z-x(n-2))*r1*r2*r3
  ELSE
     z=x3*x(i-1)-x3*x(i-2)+x(i-3)
  ENDIF
  alfbet=alf+bet
  DO its=1,MAXIT
    temp=x2+alfbet
    p1=(alf-bet+temp*z)/x2
    p2=x1
    DO j=2,n
       p3=p2
       p2=p1
       temp=x2*DBLE(j)+alfbet
       c1=x2*DBLE(j)*(DBLE(j)+alfbet)*(temp-x2)
       c2=(temp-x1)*(alf*alf-bet*bet+temp*(temp-x2)*z)
       c3=x2*(DBLE(j-1)+alf)*(DBLE(j-1)+bet)*temp
       p1=(c2*p2-c3*p3)/c1
    END DO
    pp=(DBLE(n)*(alf-bet-temp*z)*p1+x2*(DBLE(n)+alf)*(DBLE(n)+bet)*p2)/(temp*(x1-z*z))
    z1=z
    z=z1-p1/pp
    IF(ABS(z-z1).LE.EPS) EXIT
  END DO
  x(i)=z 
  ! special case alpha=beta=0 has been introduced to increase
  ! the accuracy of weights for Legendre quadrature
  if(alf /= 0.0D0 .AND. bet /= 0.0D0) then
    w(i)=DEXP(gammln(alf+DBLE(n))+gammln(bet+DBLE(n)) &
     -gammln(DBLE(n+1))-gammln(DBLE(n)+alfbet+x1))  &
     *temp*x2**alfbet/(pp*p2)
  else
    w(i)= temp*x2**alfbet/(pp*p2)/(DBLE(n))**2
  end if
END DO
RETURN
END SUBROUTINE gaujac

! FUNCTION gammln_d(XX)
! USE types
! IMPLICIT NONE
! INTEGER :: j
! REAL(DP)::GAMMLN,XX
! REAL(DP)::SER,STP,TMP,X,COF(6)
! REAL(DP)::FPF
! DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
!    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
! DATA FPF/5.5D0/
! X=XX-ONE
! TMP=X+FPF
! TMP=(X+HALF)*LOG(TMP)-TMP
! SER=ONE
! DO J=1,6
!   X=X+ONE
!   SER=SER+COF(J)/X
! END DO
! gammln_d=TMP+LOG(STP*SER)
! RETURN
! END FUNCTION gammln_d

SUBROUTINE norms(norm,nn,alf,bet)
USE types
USE base_lib, ONLY: gammln
IMPLICIT NONE
INTEGER :: n,nn
REAL(DP) :: alf, bet,lmd,x1,x2,norm1
REAL(DP) :: norm(0:nn)
REAL(DP) :: a1,a2,a3,a4
! REAL(DP), EXTERNAL :: gammln
DATA x1,x2/1.0d0,2.0d0/
lmd=alf+bet+x1
do n=0,nn
   a1=gammln(DBLE(n+1))
   a2=gammln(DBLE(n)+lmd)
   a3=gammln(DBLE(n)+alf+x1)
   a4=gammln(DBLE(n)+bet+x1)
   norm1=2**(-lmd)*(x2*DBLE(n)+lmd)*exp(a1+a2-a3-a4)
   norm(n)=SQRT(norm1)
end do
END SUBROUTINE norms


! ===========================================================================
! Evaluates associated Legendre polynomial P^m_l (x) 
! m and l are integers for which 0=< m =< l
! x is in range -1=< x =< 1
!
! this is just a wrapper around call to legendre_associated
! not be used frequently because of allocate/delallocate
!
FUNCTION plgndr_d(l,m,x)
  USE types
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: l,m
  REAL(DP), INTENT(IN) :: x
REAL(DP) :: plgndr_d
  REAL(DP), ALLOCATABLE :: cx(:)

allocate(cx(0:l))
call legendre_associated ( l, m, x, cx )
plgndr_d = cx(l)
deallocate(cx)

! do not use sign alternation
if (mod(m,2)==1) plgndr_d = -plgndr_d

END FUNCTION plgndr_d


! ===========================================================================
!
! subroutine and function by   John Burkardt
! http://www.csit.fsu.edu/~burkardt/f_src/polpak/polpak.html
!
subroutine legendre_associated ( n, m, x, cx )
!
!*******************************************************************************
!
!! LEGENDRE_ASSOCIATED evaluates the associated Legendre functions.
!
!
!  Differential equation:
!
!    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
!
!  First terms:
!
!    M = 0  ( = Legendre polynomials of first kind P(N)(X) )
!
!    P00 =    1
!    P10 =    1 X
!    P20 = (  3 X**2 -   1)/2
!    P30 = (  5 X**3 -   3 X)/2
!    P40 = ( 35 X**4 -  30 X**2 +   3)/8
!    P50 = ( 63 X**5 -  70 X**3 +  15 X)/8
!    P60 = (231 X**6 - 315 X**4 + 105 X**2 -  5)/16
!    P70 = (429 X**7 - 693 X**5 + 315 X**3 - 35 X)/16
!
!    M = 1
!
!    P01 =   0
!    P11 =   1 * SQRT(1-X*X)
!    P21 =   3 * SQRT(1-X*X) * X
!    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
!    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
!
!    M = 2
!
!    P02 =   0
!    P12 =   0
!    P22 =   3 * (1-X*X)
!    P32 =  15 * (1-X*X) * X
!    P42 = 7.5 * (1-X*X) * (7*X*X-1)
!
!    M = 3
!
!    P03 =   0
!    P13 =   0
!    P23 =   0
!    P33 =  15 * (1-X*X)**1.5
!    P43 = 105 * (1-X*X)**1.5 * X
!
!    M = 4
!
!    P04 =   0
!    P14 =   0
!    P24 =   0
!    P34 =   0
!    P44 = 105 * (1-X*X)**2
!
!  Recursion:
!
!    if N < M:
!      P(N,M) = 0
!    if N = M:
!      P(N,M) = (2*M-1)!! * (1-X*X)**(M/2) where N!! means the product of
!      all the odd integers less than or equal to N.
!    if N = M+1:
!      P(N,M) = X*(2*M+1)*P(M,M)
!    if M+1 < N:
!      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
!
!  Restrictions:
!
!    -1 <= X <= 1
!     0 <= M <= N
!
!  Special values:
!
!    P(N,0)(X) = P(N)(X), that is, for M=0, the associated Legendre
!    function of the first kind equals the Legendre polynomial of the
!    first kind.
!
!  Modified:
!
!    14 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz and Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964.
!
!  Parameters:
!
!    Input, integer N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X, the point at which the function is to be
!    evaluated.  X must satisfy -1 <= X <= 1.
!
!    Output, real ( kind = 8 ) CX(0:N), the values of the first N+1 functions.
!
  implicit none
!
  integer n
!
  real ( kind = 8 ) cx(0:n)
  real ( kind = 8 ) fact
  integer i
  integer m
  real ( kind = 8 ) somx2
  real ( kind = 8 ) x
!
  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,i6)' ) '  Input value of M is ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop
  end if
 
  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,i6)' ) '  Input value of M = ', m
    write ( *, '(a,i6)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but M must be less than or equal to N.'
    stop
  end if
 
  if ( x < -1.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no less than -1.'
    stop
  end if
 
  if ( 1.0D+00 < x ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED - Fatal error!'
    write ( *, '(a,g14.6)' ) '  Input value of X = ', x
    write ( *, '(a)' ) '  but X must be no more than 1.'
    stop
  end if
  
  cx(0:m-1) = 0.0D+00

  cx(m) = 1.0D+00
  somx2 = sqrt ( 1.0D+00 - x * x )
 
  fact = 1.0D+00
  do i = 1, m
    cx(m) = - cx(m) * fact * somx2
    fact = fact + 2.0D+00
  end do
 
  if ( m == n ) then
    return
  end if

  cx(m+1) = x * real ( 2 * m + 1, kind = 8 ) * cx(m)

  do i = m+2, n
    cx(i) = ( real ( 2 * i     - 1, kind = 8 ) * x * cx(i-1) &
            + real (   - i - m + 1, kind = 8 ) *     cx(i-2) ) &  
            / real (     i - m,     kind = 8 )
  end do

  return
end

! ===========================================================================
!
function gamma_log ( x )
!
!*******************************************************************************
!
!! GAMMA_LOG calculates the natural logarithm of GAMMA ( X ) for positive X.
!
!
!  Discussion:
!
!    Computation is based on an algorithm outlined in references 1 and 2.
!    The program uses rational functions that theoretically approximate
!    log ( GAMMA(X) ) to at least 18 significant decimal digits.  The
!    approximation for 12 < X is from reference 3, while approximations
!    for X < 12.0 are similar to those in reference 1, but are unpublished.
!    The accuracy achieved depends on the arithmetic system, the compiler,
!    intrinsic functions, and proper selection of the machine-dependent
!    constants.
!
!  Modified:
!
!    16 June 1999
!
!  Authors:
!
!    W. J. Cody and L. Stoltz
!    Argonne National Laboratory
!
!  Reference:
!
!    # 1)
!    W. J. Cody and K. E. Hillstrom,
!    Chebyshev Approximations for the Natural Logarithm of the Gamma Function,
!    Mathematics of Computation,
!    Volume 21, 1967, pages 198-203.
!
!    # 2)
!    K. E. Hillstrom,
!    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
!    May 1969.
!
!    # 3)
!    Hart, Cheney, Lawson, Maehly, Mesztenyi, Rice, Thacher, Witzgall,
!    Computer Approximations,
!    Wiley, 1968.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the Gamma function. 
!    X must be positive.
!
!    Output, real ( kind = 8 ) GAMMA_LOG, the logarithm of the Gamma 
!    function of X.  If X <= 0.0, or if overflow would occur, the
!    program returns the value HUGE().
!
!  Machine-dependent constants:
!
!    BETA   - radix for the floating-point representation.
!
!    MAXEXP - the smallest positive power of BETA that overflows.
!
!    XBIG   - largest argument for which LN(GAMMA(X)) is representable
!             in the machine, i.e., the solution to the equation
!             LN(GAMMA(XBIG)) = BETA**MAXEXP.
!
!    XINF   - largest machine representable floating-point number;
!             approximately BETA**MAXEXP.
!
!    FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!    Approximate values for some important machines are:
!
!                              BETA      MAXEXP         XBIG
!
!    CRAY-1        (S.P.)        2        8191       9.62D+2461
!    Cyber 180/855
!      under NOS   (S.P.)        2        1070       1.72D+319
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)        2         128       4.08D+36
!    IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)        2        1024       2.55D+305
!    IBM 3033      (D.P.)       16          63       4.29D+73
!    VAX D-Format  (D.P.)        2         127       2.05D+36
!    VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                            FRTBIG
!
!    CRAY-1        (S.P.)   3.13D+615
!    Cyber 180/855
!      under NOS   (S.P.)   6.44D+79
!    IEEE (IBM/XT,
!      SUN, etc.)  (S.P.)   1.42D+9
!    IEEE (IBM/XT,
!      SUN, etc.)  (D.P.)   2.25D+76
!    IBM 3033      (D.P.)   2.56D+18
!    VAX D-Format  (D.P.)   1.20D+9
!    VAX G-Format  (D.P.)   1.89D+76
!
  implicit none
!
  real ( kind = 8 ), parameter, dimension ( 7 ) :: c = (/ &
    -1.910444077728D-03, &
     8.4171387781295D-04, &
    -5.952379913043012D-04, &
     7.93650793500350248D-04, &
    -2.777777777777681622553D-03, &
     8.333333333333333331554247D-02, &
     5.7083835261D-03 /)
  real ( kind = 8 ) corr
  real ( kind = 8 ), parameter :: d1 = - 5.772156649015328605195174D-01
  real ( kind = 8 ), parameter :: d2 =   4.227843350984671393993777D-01
  real ( kind = 8 ), parameter :: d4 =   1.791759469228055000094023D+00
  real ( kind = 8 ) eps
  real ( kind = 8 ), parameter :: frtbig = 1.42D+09
  integer i
  real ( kind = 8 ) gamma_log
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p1 = (/ &
    4.945235359296727046734888D+00, &
    2.018112620856775083915565D+02, &
    2.290838373831346393026739D+03, &
    1.131967205903380828685045D+04, &
    2.855724635671635335736389D+04, &
    3.848496228443793359990269D+04, &
    2.637748787624195437963534D+04, &
    7.225813979700288197698961D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p2 = (/ &
    4.974607845568932035012064D+00, &
    5.424138599891070494101986D+02, &
    1.550693864978364947665077D+04, &
    1.847932904445632425417223D+05, &
    1.088204769468828767498470D+06, &
    3.338152967987029735917223D+06, &
    5.106661678927352456275255D+06, &
    3.074109054850539556250927D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: p4 = (/ &
    1.474502166059939948905062D+04, &
    2.426813369486704502836312D+06, &
    1.214755574045093227939592D+08, &
    2.663432449630976949898078D+09, &
    2.940378956634553899906876D+10, &
    1.702665737765398868392998D+11, &
    4.926125793377430887588120D+11, &
    5.606251856223951465078242D+11 /)
  real ( kind = 8 ), parameter :: pnt68 = 0.6796875D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q1 = (/ &
    6.748212550303777196073036D+01, &
    1.113332393857199323513008D+03, &
    7.738757056935398733233834D+03, &
    2.763987074403340708898585D+04, &
    5.499310206226157329794414D+04, &
    6.161122180066002127833352D+04, &
    3.635127591501940507276287D+04, &
    8.785536302431013170870835D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q2 = (/ &
    1.830328399370592604055942D+02, &
    7.765049321445005871323047D+03, &
    1.331903827966074194402448D+05, &
    1.136705821321969608938755D+06, &
    5.267964117437946917577538D+06, &
    1.346701454311101692290052D+07, &
    1.782736530353274213975932D+07, &
    9.533095591844353613395747D+06 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: q4 = (/ &
    2.690530175870899333379843D+03, &
    6.393885654300092398984238D+05, &
    4.135599930241388052042842D+07, &
    1.120872109616147941376570D+09, &
    1.488613728678813811542398D+10, &
    1.016803586272438228077304D+11, &
    3.417476345507377132798597D+11, &
    4.463158187419713286462081D+11 /)
  real ( kind = 8 ) res
  real ( kind = 8 ), parameter :: sqrtpi = 0.9189385332046727417803297D+00
  real ( kind = 8 ) x
  real ( kind = 8 ), parameter :: xbig = 4.08D+36
  real ( kind = 8 ) xden
  real ( kind = 8 ) xm1
  real ( kind = 8 ) xm2
  real ( kind = 8 ) xm4
  real ( kind = 8 ) xnum
  real ( kind = 8 ) xsq
!
!  Return immediately if the argument is out of range.
!
  if ( x <= 0.0D+00 .or. xbig < x ) then
    gamma_log = huge ( gamma_log )
    return
  end if

  eps = epsilon ( eps )

  if ( x <= eps ) then

    res = - log ( x )

  else if ( x <= 1.5D+00 ) then

    if ( x < pnt68 ) then
      corr = - log ( x )
      xm1 = x
    else
      corr = 0.0D+00
      xm1 = ( x - 0.5D+00 ) - 0.5D+00
    end if

    if ( x <= 0.5D+00 .or. pnt68 <= x ) then

      xden = 1.0D+00
      xnum = 0.0D+00

      do i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      end do

      res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

    else

      xm2 = ( x - 0.5D+00 ) - 0.5D+00
      xden = 1.0D+00
      xnum = 0.0D+00
      do i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      end do

      res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

    end if

  else if ( x <= 4.0D+00 ) then

    xm2 = x - 2.0D+00
    xden = 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    end do

    res = xm2 * ( d2 + xm2 * ( xnum / xden ) )

  else if ( x <= 12.0D+00 ) then

    xm4 = x - 4.0D+00
    xden = - 1.0D+00
    xnum = 0.0D+00
    do i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    end do

    res = d4 + xm4 * ( xnum / xden )

  else

    res = 0.0D+00

    if ( x <= frtbig ) then

      res = c(7)
      xsq = x * x

      do i = 1, 6
        res = res / xsq + c(i)
      end do

    end if

    res = res / x
    corr = log ( x )
    res = res + sqrtpi - 0.5D+00 * corr
    res = res + x * ( corr - 1.0D+00 )

  end if

  gamma_log = res

  return
end

! ===========================================================================
!
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


! ===========================================================================
! interfaces to LAPACK diagonalizers
!
! regular diagonalisers

SUBROUTINE diag(ham,eng,nlev,nout,icut,encut,istage,idiag)
  USE types
!  USE param, ONLY: test_flag
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEV, LA_SYEVD, LA_SYEVX, LA_SYEVR
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)  :: nlev
  INTEGER(I4B), INTENT(OUT) :: nout,idiag
  INTEGER(I4B), INTENT(IN)  :: icut,istage
  REAL(DP), INTENT(INOUT) :: ham(nlev,nlev)
  REAL(DP), INTENT(OUT)   :: eng(nlev)
  REAL(DP), INTENT(IN)    :: encut
! INTEGER(I4B) idiag

! if(test_flag >2) then
!  write(u6,'(" main: stage",i2,": diagonalizing...")') istage
! end if

nout = nlev

select case (icut)
case(0)
  ! find all levels
  !  write(u6,'("  finding all levels")')
  call la_syev(ham,eng,JOBZ='V',UPLO='U',INFO=idiag)
case(1:)
  ! find icut1 levels
  !  write(u6,'("  finding ",i4," lowest levels")')  min(icut,nlev)
  ! call la_syevx(ham,eng,JOBZ='V',UPLO='U',IU=min(icut,nlev),M=nout,INFO=idiag)
  call la_syevr(ham,eng,JOBZ='V',UPLO='U',IU=min(icut,nlev),M=nout,INFO=idiag)
case(:-1)
  ! find levels E < encut
  !  write(u6,'("  finding lowest levels with E<",f12.3)')  encut
  ! call la_syevx(ham,eng,JOBZ='V',UPLO='U',VU=encut,M=nout,INFO=idiag)
  call la_syevr(ham,eng,JOBZ='V',UPLO='U',VU=encut,M=nout,INFO=idiag)
end select

if(idiag /= 0) then
  write(u6,'(" at stage",i4,"  idiag=",i4)') istage,idiag
  stop
end if

! if(test_flag >2) then
!   write(u6,'(" main: stage",i2,": diagonalization done")') istage
! end if

END SUBROUTINE

! (almost) same as previous but ham is in packed form
!
!  these are regular diagonalisers
! USE F95_LAPACK, ONLY: LA_SYEV, LA_SYEVD, LA_SYEVX, LA_SYEVR
!  these are packed diagonalisers
! USE F95_LAPACK, ONLY: LA_SPEV, LA_SPEVD, LA_SPEVX
!  this is my "fixed" SPEVX

SUBROUTINE diag_p(ham,eng,vectors,nlev,nout,icut,encut,istage,idiag)
  USE types
!  USE param, ONLY: test_flag
  USE LA_PRECISION, ONLY: WP => DP
!  USE my_lapack95
  USE F95_LAPACK, ONLY: LA_SPEVX
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)  :: nlev
  INTEGER(I4B), INTENT(OUT) :: nout,idiag
  INTEGER(I4B), INTENT(IN)  :: icut,istage
  REAL(DP), INTENT(INOUT) :: ham((nlev*(nlev+1))/2)
  REAL(DP), INTENT(OUT)   :: eng(nlev)
  REAL(DP), INTENT(OUT)   :: vectors(nlev,icut)
  REAL(DP), INTENT(IN)    :: encut
! INTEGER(I4B) idiag
INTEGER(I4B) ii

! if(test_flag >2) then
!   write(u6,'(" main: stage",i2,": diagonalizing...")') istage
! end if

ii = min(icut,nlev)

select case (icut)
case(0)
  ! find all levels
  write(u6,'("  option 0: no eigenvectors")')
  call la_spevx(ham,eng,UPLO='U',IU=ii,M=nout,INFO=idiag)
case(1:)
  ! find icut levels
  !  write(u6,'("  finding ",i4," lowest levels")')  min(icut,nlev)
  call la_spevx(ham,eng,UPLO='U',Z=vectors(:,1:ii),IU=ii,M=nout,INFO=idiag)
case(:-1)
  ! find levels E < encut
  !  write(u6,'("  finding lowest levels with E<",f12.3)')  encut
  call la_spevx(ham,eng,UPLO='U',Z=vectors(:,1:ii),VU=encut,M=nout,INFO=idiag)
end select

if(idiag /= 0) then
  write(u6,'(" at stage",i4,"  idiag=",i4)') istage,idiag
  stop
end if

! if(test_flag >2) then
!   write(u6,'(" main: stage",i2,": diagonalization done")') istage
! end if

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
