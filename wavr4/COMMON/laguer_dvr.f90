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
  USE potential
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEV
  INTEGER(I4B), INTENT(IN) :: nq
  REAL(DP), INTENT(OUT) :: r(nn(nq)),hdvr(nn(nq),nn(nq)), &
                           oner(nn(nq),nn(nq))
  END SUBROUTINE
END INTERFACE

INTERFACE stretch2
  SUBROUTINE stretch2(nq,r,hdvr,oner)
  USE types
  USE param
  USE base_lib, ONLY: lagpoly, gammln
  USE laguer_grid, ONLY: laguer
  USE potential
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEV
  INTEGER(I4B), INTENT(IN) :: nq
  REAL(DP), INTENT(OUT) :: r(nn(nq)),hdvr(nn(nq),nn(nq)), &
                           oner(nn(nq),nn(nq))
  END SUBROUTINE
END INTERFACE

END MODULE

SUBROUTINE stretch(nq,r,hdvr,oner)
  USE types
  USE param
  USE base_lib, ONLY: lagpoly, gammln
  USE laguer_grid, ONLY: laguer
  USE potential
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEV
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: nq
  ! REAL(DP), INTENT(OUT) :: r(:),hdvr(:,:)
  REAL(DP), INTENT(OUT) :: r(nn(nq)),hdvr(nn(nq),nn(nq)), &
                           oner(nn(nq),nn(nq))
REAL(DP), ALLOCATABLE :: tmp(:,:),y(:),w(:),t(:,:),e1(:)
REAL(DP) :: csx,csa,tsx,tsa, wt, q(6),temp
REAL(DP) :: toler = 1.0d-14
REAL(DP) :: alf, beta, ec
INTEGER(I4B) :: i,j,nx,n

write(u7,'(" stretch: processing r",i1,"...")') nq

! beta in Angstroms Eqn. (38)
! beta = we*sqrt(c*mu*Pi/De/hbar)
!   where we - omega_e in cm-1; c - speed of light in cm/sec; 
!   mu - reduced mass in kg; De - dissoc. energy in cm-1;
!   hbar = h/2/Pi in J*s. After feeding all the constants, for
!   mu in au and beta in 1/Angstrom:

write(u6,'("  DVR parameters:  re=",f10.6,";  we=",f12.3,";  De=",f12.3)') &
  re(nq),we(nq),De(nq)
write(u6,'("  Equilibr./Ref.: Qe",i1,"=",f10.6)') nq,qe(nq)

beta = 0.1217788108405819_dp *we(nq)*sqrt(muq(nq)/De(nq))
ec = we(nq)*we(nq)/De(nq)/16.d0
alf = 4.d0*De(nq)/we(nq)

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

write(u7,'(" grid sum    should be         is               difference")')
write(u7,'("  r:",3d20.12)') tsx, csx, tsx-csx
write(u7,'("  w:",3d20.12)') tsa, csa, tsa-csa

write(u7, *) ' points  relative accuracy: ', dabs((csx-tsx)/tsx)
write(u7, *) ' weights relative accuracy: ', dabs((csa-tsa)/tsa)

if (abs((csx-tsx)/tsx) > toler) then 
  write(u7, *) ' laguer: points in error, adjust algorithm '
  stop
end if

if (abs((csa-tsa)/tsa) > toler) then 
  write(u7, *) ' laguer: weights in error, adjust algorithm '
  stop 
end if

deallocate(w,y)

write(u7,'(" stretch: grid & weights done.")')

allocate(tmp(nx,nx))
hdvr = 0.0_dp

! make hfbr
!
! implementation of formular by JT & BTS 77, 4065 JCP 1982
! Eqn (45) note the sign error!
!
do i=1,nx
  !
  ! n=n'=i-1 (because n starts from 0)
  !
  n=i-1
  hdvr(i,i)= ec*(2.d0*n*(alf+n+1.d0)+alf+1.d0)
  !
  ! n=n'+2=i+1
  !
  if((i+2)<=nx) then 
    n=i+1
    hdvr(i+2,i)=-ec*sqrt((alf+n)*(alf+n-1.d0)*n*(n-1.d0))
    hdvr(i,i+2)=hdvr(i+2,i)
  end if
end do

tmp  = matmul(transpose(t),hdvr)
hdvr = matmul(tmp,t)

deallocate(tmp)

allocate(e1(nx))
t=hdvr
q=qe
do i=1,nx
  q(nq)=r(i)
  temp=v(q(1),q(2),q(3),q(4),q(5),q(6))
  ! write(u7,'(1x,i3,f10.6,f14.3)') i,r(i),temp
  t(i,i)=t(i,i)+temp
end do
write(u7,'(" stretch: solving 1D...")')
call la_syev(t,e1,JOBZ='N',UPLO='U')
write(u7,'(" stretch: eigenvalues:",/)')
write(u7,'("   #        E            r            v")')
do i=1,nx
  q(nq)=r(i)
  temp=v(q(1),q(2),q(3),q(4),q(5),q(6))
  write(u7,'(1x,i3,f12.3,f14.6,f14.3)') i,e1(i),r(i),temp
end do
if(nx>1)   write(u7,'(/,"   we=",f12.3)') e1(2)-e1(1)
write(6,*) ' trace=',sum(e1)
deallocate(e1)
deallocate(t)

write(u7,'("  compare exact and DVR 1/R^2")')
!  assuming:     q3 -> R,       q1 -> d1,       q2 -> d2
!  resp.      mu_q3 -> mu_R, mu_q1 -> mu_d1, mu_q2 -> mu_d2

temp  = convf/two/muq(nq)

oner = zero

do i=1,nx
  write(u7,'(1x,i3,2f14.6)') i,oner(i,i),temp/r(i)/r(i)
  oner(i,i) = temp/r(i)/r(i)
end do

write(u7,'(" stretch: done.",/)')

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
  USE potential
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEV
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: nq
  ! REAL(DP), INTENT(OUT) :: r(:),hdvr(:,:)
  REAL(DP), INTENT(OUT) :: r(nn(nq)),hdvr(nn(nq),nn(nq)), &
                           oner(nn(nq),nn(nq))
REAL(DP), ALLOCATABLE :: tmp(:,:),y(:),w(:),t(:,:),e1(:)
REAL(DP), ALLOCATABLE :: y1(:),w1(:)
REAL(DP) :: csx,csa,tsx,tsa, wt, q(6),temp,tempr
REAL(DP) :: toler = 1.0d-13
REAL(DP) :: alf, beta, ec, eta, norm
INTEGER(I4B) :: i,i1,nx,n,n1, j

write(u7,'(" stretch: processing r",i1,"...")') nq

write(u6,'("  DVR parameters:  we=",f12.3)') we(nq)
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

write(u7,'(" grid sum    should be         is               difference")')
write(u7,'("  r:",3d20.12)') tsx, csx, tsx-csx
write(u7,'("  w:",3d20.12)') tsa, csa, tsa-csa

write(u7, *) ' points  relative accuracy: ', dabs((csx-tsx)/tsx)
write(u7, *) ' weights relative accuracy: ', dabs((csa-tsa)/tsa)

if (abs((csx-tsx)/tsx) > toler) then
  write(u7, *) ' laguer: points in error, adjust algorithm '
  stop
end if

if (abs((csa-tsa)/tsa) > toler) then
  write(u7, *) ' laguer: weights in error, adjust algorithm '
  stop
end if

deallocate(w,y,w1,y1)

write(u7,'(" stretch: grid & weights done.")')

hdvr = 0.0_dp
oner = 0.0_dp

! matrix elements of d^2/dr^2  in FBR reps
! can include plus extra 1/r^2 term
! Note: n runs from 0 to nn-1
!
! implementation of Eqs (21),(22) by JT & BTS 101, p71 JMS 1983
! Note:
! N{m,eta} are taken to be Laguerre norms;
! Eqs (21),(22) are multiplied by sqrt(2*sqrt(zeta)) to get the right norm
! zeta=beta here
! then the factor hbar^2/(4*mu)*sqrt(zeta)  -> hbar^2/(2*mu)*zeta
! which turns out to be simply we/2 where we is in cm-1
! where zeta=beta= mu*we/hbar
!
! <n1|H|n>
do i=1,nx
  n=i-1      ! right index

  do i1=1,i
    n1=i1-1  ! left index

    norm= sqrt(exp(gammln(n1+one)-gammln(alf+n1+one)+ &
                   gammln(n+one)-gammln(alf+n+one)))

    temp= 0.0D0
    do j=0,n1
      temp= temp+exp(gammln(j+eta+half)-gammln(j+one))
    end do
    temp= -eta*(eta+one)*temp

    if(n==n1) then
      temp= temp+ (2*n+eta+1.5D0)*exp(gammln(n+eta+1.5D0)-gammln(n+one))
    end if

    if(n1==n-1) then
      temp= temp+ exp(gammln(n1+eta+2.5D0)-gammln(n1+one))
    end if

    hdvr(i1,i) = norm*temp*ec

!========  1/R^2  term  ========
    tempr= 0.0D0
    do j=0,n1
      tempr= tempr+exp(gammln(j+eta+half)-gammln(j+one))
    end do
    oner(i1,i) = norm*ec*tempr
!===============================
!   the exact spherical oscillator problem has the term
!   eta*(eta+one)/R^2 and therefore the conribution is
!    oner(i1,i) = norm*ec*eta*(eta+one)*tempr

    if(i1 /= i) then
      hdvr(i,i1) = hdvr(i1,i)
      oner(i,i1) = oner(i1,i)
    end if

  end do

end do

allocate(tmp(nx,nx))

! transform K.E. matrix to DVR repr

tmp  = matmul(transpose(t),hdvr)
hdvr = matmul(tmp,t)

! transform 1/R^2 matrix to DVR repr

tmp  = matmul(transpose(t),oner)
oner = matmul(tmp,t)

deallocate(tmp)

! write(u7,'("   #       r          v")')
allocate(e1(nx))
t=hdvr
q=qe
do i=1,nx
  q(nq)=r(i)
  temp=v(q(1),q(2),q(3),q(4),q(5),q(6))
  ! write(u7,'(1x,i3,f10.6,f14.3)') i,r(i),temp
  t(i,i)=t(i,i)+temp
end do
write(u7,'(" stretch: solving 1D...")')
call la_syev(t,e1,JOBZ='N',UPLO='U')
write(u7,'(" stretch: eigenvalues:",/)')
write(u7,'("   #        E            r            v")')
do i=1,nx
  q(nq)=r(i)
  temp=v(q(1),q(2),q(3),q(4),q(5),q(6))
  write(u7,'(1x,i3,f12.3,f14.6,f14.3)') i,e1(i),r(i),temp
end do
if(nx>1)   write(u7,'(/,"   we=",f12.3)') e1(2)-e1(1)
write(6,*) ' trace=',sum(e1)
deallocate(e1)
deallocate(t)

write(u7,'("  compare exact and DVR 1/R^2")')
!  assuming:     q3 -> R,       q1 -> d1,       q2 -> d2
!  resp.      mu_q3 -> mu_R, mu_q1 -> mu_d1, mu_q2 -> mu_d2

temp  = convf/two/muq(nq)

!  the IF structure below to be used in conjunction with removal of
!  off-diagonal contributions of 1/R^2 (see main.f90)
if(oner_flag==0) then
  do i=1,nx
    oner(i,i) = temp/r(i)/r(i)
    do i1=1,i-1
      oner(i,i1)=zero
      oner(i1,i)=zero
    end do
  end do
end if


do i=1,nx
  write(u7,'(1x,i3,2f14.6)') i,oner(i,i),temp/r(i)/r(i)
end do

write(u7,'(" stretch: done.",/)')

END SUBROUTINE
