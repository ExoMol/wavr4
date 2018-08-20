! 22/04/2004  toler was increased from 1E-14 to 5E-14
! stretch -- morse-laguerre basis
! stretch2 -- spherical oscillator basis
!   20/06/02: corrections in <n1|H|n>
!   07/07/02: analytic matrix elements for 1/R^2

MODULE laguer_dvr

INTERFACE stretch
  SUBROUTINE stretch(nq,r,hdvr,oner)
  USE types
  USE param
  USE base_lib, ONLY: lagpoly, gammln
  USE laguer_grid, ONLY: laguer
!  USE potential
!  USE LA_PRECISION, ONLY: WP => DP
!  USE F95_LAPACK, ONLY: LA_SYEV
  INTEGER(I4B), INTENT(IN) :: nq,hdvr(nn(nq),nn(nq)),oner(nn(nq),nn(nq))
  REAL(DP), INTENT(OUT) :: r(nn(nq))

  END SUBROUTINE
END INTERFACE

INTERFACE stretch2
  SUBROUTINE stretch2(nq,r,hdvr,oner)
  USE types
  USE param
  USE base_lib, ONLY: lagpoly, gammln
  USE laguer_grid, ONLY: laguer
!  USE potential
!  USE LA_PRECISION, ONLY: WP => DP
!  USE F95_LAPACK, ONLY: LA_SYEV
  INTEGER(I4B), INTENT(IN) :: nq,hdvr(nn(nq),nn(nq)),oner(nn(nq),nn(nq))
  REAL(DP), INTENT(OUT) :: r(nn(nq))

  END SUBROUTINE
END INTERFACE

END MODULE

SUBROUTINE stretch(nq,r,hdvr,oner)
  USE types
  USE param
  USE base_lib, ONLY: lagpoly, gammln
  USE laguer_grid, ONLY: laguer
!  USE potential
!  USE LA_PRECISION, ONLY: WP => DP
!  USE F95_LAPACK, ONLY: LA_SYEV
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: nq,hdvr(nn(nq),nn(nq)),oner(nn(nq),nn(nq)) ! hdvr and oner are not actually used
  REAL(DP), INTENT(OUT) :: r(nn(nq))

REAL(DP), ALLOCATABLE :: tmp(:,:),y(:),w(:),t(:,:),e1(:)
REAL(DP) :: csx,csa,tsx,tsa, wt, q(6),temp
REAL(DP) :: toler = 5.0d-14
REAL(DP) :: alf, beta, ec
INTEGER(I4B) :: i,j,nx,n

if(test_flag > 0) write(u7,'(" stretch: processing r",i1,"...")') nq

! beta in Angstroms Eqn. (38)
! beta = we*sqrt(c*mu*Pi/De/hbar)
!   where we - omega_e in cm-1; c - speed of light in cm/sec; 
!   mu - reduced mass in kg; De - dissoc. energy in cm-1;
!   hbar = h/2/Pi in J*s. After feeding all the constants, for
!   mu in au and beta in 1/Angstrom:

if(test_flag > 1) &
  write(u6,'("  DVR parameters:  re=",f10.6,";  we=",f12.3,";  De=",f12.3)') &
  re(nq),we(nq),De(nq)
if(test_flag > 1) &
  write(u6,'("  Equilibr./Ref.: Qe",i1,"=",f10.6)') nq,qe(nq)

beta = 0.1217788108405819_dp *we(nq)*sqrt(muq(nq)/De(nq))
ec = we(nq)*we(nq)/De(nq)/16.d0
alf = 4.d0*De(nq)/we(nq)

if(test_flag > 1) &
  write(u7,'(" beta=",f8.2, " alpha=",f8.2)') beta,alf

nx = nn(nq)

allocate(t(nx,nx))
allocate(y(nx),w(nx))

call laguer(nx,y,w,alf,csx,csa,tsx,tsa)

! n runs over basis functions from 1 to nx (in fact from 0 to nx-1)
! j over the grid points      from 1 to nx

if(test_flag > 1) &
  write(u7,*) '    j          y(j)              r(j)               w(j)'
do j=1,nx
  r(j)=log(alf/y(j))/beta+re(nq)
  if(test_flag > 1)  write(u7,'(2x,i4,2f18.12,e20.12)') j, y(j), r(j), w(j)
end do

do n=1,nx
  wt=exp(0.5d0*(gammln(nx+alf)-gammln(n+alf)+gammln(real(n,dp))-gammln(nx+1.d0)))
  do j=1,nx
    t(n,j)= sqrt(w(j))*wt*lagpoly(n-1,alf,y(j))
  end do
end do

if(test_flag > 1) then
  write(u7,'(" grid sum    should be         is               difference")')
  write(u7,'("  r:",3d20.12)') tsx, csx, tsx-csx
  write(u7,'("  w:",3d20.12)') tsa, csa, tsa-csa

  write(u7, *) ' points  relative accuracy: ', dabs((csx-tsx)/tsx)
  write(u7, *) ' weights relative accuracy: ', dabs((csa-tsa)/tsa)
endif

if (abs((csx-tsx)/tsx) > toler) then 
  write(u7, *) ' laguer: points in error, adjust algorithm '
  stop
end if

if (abs((csa-tsa)/tsa) > toler) then 
  write(u7, *) ' laguer: weights in error, adjust algorithm '
  stop 
end if

deallocate(w,y)

! hdvr  calculations removed

if(test_flag > 0) write(u7,'(" stretch: done.",/)')

END SUBROUTINE


!  spherical oscillator
!  input: referencing index 1, 2, or 3
!  output: r     grid,
!          hdvr  kinetic energy matrix in DVR repr
!          oner  hbar^2/2/mu(nq)/R^2 in DVR repr
!                either exact or in DVR approx

SUBROUTINE stretch2(nq,r,hdvr,oner)
  USE types
  USE param
  USE base_lib, ONLY: lagpoly, gammln
  USE laguer_grid, ONLY: laguer
!  USE potential
!  USE LA_PRECISION, ONLY: WP => DP
!  USE F95_LAPACK, ONLY: LA_SYEV
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: nq,hdvr(nn(nq),nn(nq)),oner(nn(nq),nn(nq)) ! hdvr and oner are not actually used
  REAL(DP), INTENT(OUT) :: r(nn(nq))

REAL(DP), ALLOCATABLE :: tmp(:,:),y(:),w(:),t(:,:),e1(:)
REAL(DP), ALLOCATABLE :: y1(:),w1(:)
REAL(DP) :: csx,csa,tsx,tsa, wt, q(6),temp,tempr
REAL(DP) :: toler = 1.0d-13
REAL(DP) :: alf, beta, ec, eta, norm
INTEGER(I4B) :: i,i1,nx,n,n1, j

if(test_flag > 0) write(u7,'(" stretch: processing r",i1,"...")') nq

if(test_flag > 1) &
  write(u6,'("  DVR parameters:  we=",f12.3)') we(nq)
if(test_flag > 1) &
  write(u6,'("  Equilibr./Ref.: Qe",i1,"=",f10.6)') nq,qe(nq)

! eta=zero
! eta=one
 eta= re(nq)
write(u6,*) '  eta=', eta
alf = eta+half
beta = muq(nq)*we(nq)* atomic_mass*TWOPI*speed_of_light/(hbar) *angstrom**2
ec = we(nq)*half
! g = we(nq)**2*PI*speed_of_light*muq(nq)*atomic_mass/(hbar) *angstrom**2
! g2 = hbar/(4*PI*atomic_mass*speed_of_light*angstrom**2) /muq(nq)
! g2 = g2*eta*(eta+one)
! print *, '    g=', g, '    g2=', g2

if(test_flag > 1) &
  write(u7,'(" beta=",f8.2, " alpha=",f8.2)') beta,alf

nx = nn(nq)

allocate(t(nx,nx))
allocate(y(nx),w(nx),y1(nx),w1(nx))

call laguer(nx,y1,w1,alf,csx,csa,tsx,tsa)

! reverse the order of the points - consistency with Morse-Laguerre
do j=1,nx
  y(j)= y1(nx+1-j)
  w(j)= w1(nx+1-j)
end do

! n runs over basis functions from 1 to nx (in fact from 0 to nx-1)
! j over the grid points      from 1 to nx

if(test_flag > 1) &
  write(u7,*) '    j          y(j)              r(j)               w(j)'
do j=1,nx
  r(j)=sqrt(y(j)/beta)
  if(test_flag > 1)  write(u7,'(2x,i4,2f18.12,e20.12)') j, y(j), r(j), w(j)
end do

do n=1,nx
  wt=exp(0.5d0*(gammln(nx+alf)-gammln(n+alf)+gammln(real(n,dp))-gammln(nx+1.d0)))
  do j=1,nx
    t(n,j)= sqrt(w(j))*wt*lagpoly(n-1,alf,y(j))
  end do
end do

if(test_flag > 1) then
  write(u7,'(" grid sum    should be         is               difference")')
  write(u7,'("  r:",3d20.12)') tsx, csx, tsx-csx
  write(u7,'("  w:",3d20.12)') tsa, csa, tsa-csa

  write(u7, *) ' points  relative accuracy: ', dabs((csx-tsx)/tsx)
  write(u7, *) ' weights relative accuracy: ', dabs((csa-tsa)/tsa)
endif

if (abs((csx-tsx)/tsx) > toler) then
  write(u7, *) ' laguer: points in error, adjust algorithm '
  stop
end if

if (abs((csa-tsa)/tsa) > toler) then
  write(u7, *) ' laguer: weights in error, adjust algorithm '
  stop
end if

deallocate(w,y,w1,y1)

! hdvr  calculations removed

if(test_flag > 0) write(u7,'(" stretch: done.",/)')

END SUBROUTINE
