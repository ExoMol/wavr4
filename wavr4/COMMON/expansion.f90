! 19/04/02
! purpose: expand angular potential for every radial grid point
! input:   allocated array vex and radial grid points
! input (implicit): expansion parameters are provided via USE param
!   full grid sizes:             theta1 - fg1, theta - fg2, phi - fg3
!   symmetry reduced grid sizes: theta1 - rg1, theta - rg2, phi - rg3
!   jmax - ne1, lmax - ne2, kmax - ne3
!   note: kmax = min(jmax,lmax,kmax)
! output:  expansion coefficients vex() computed at the grid points

MODULE expansion

INTERFACE vexpand
SUBROUTINE vexpand(vex,q1,q2,q3,testpot)
  USE types
  USE param
  USE base_lib, ONLY: plgndr, gammln, gauleg
!  USE potential
  USE workarrays, ONLY: v3
! IMPLICIT NONE
  REAL(DP), INTENT(IN) :: q1(1:nn(1)), q2(1:nn(2)), q3(1:nn(3))
  REAL(DP), INTENT(OUT) :: vex(0:ne1,0:ne2,0:ne3,nn(1),nn(2),nn(3))
  LOGICAL(LGT), INTENT(OUT) :: testpot(nn(1),nn(2),nn(3))
  END SUBROUTINE
END INTERFACE

INTERFACE threej
FUNCTION threej(j1,j3,m1,m3,j2min,j2max,j2size,j2sb)
  USE types
! IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: j1,j3,m1,m3,j2min,j2max,j2size,j2sb
  REAL(DP) :: threej(j2size)
  END FUNCTION threej
END INTERFACE

INTERFACE setfac
  SUBROUTINE setfac()
  USE types
  USE param
  USE workarrays, ONLY: threej0
! IMPLICIT NONE
  END SUBROUTINE setfac
END INTERFACE

INTERFACE pot_sym
SUBROUTINE pot_sym(r1,r2,r3, theta1,theta,phi)
  USE types
  USE param
  USE potential
  USE workarrays, ONLY: v3
! IMPLICIT NONE
  REAL(DP), INTENT(IN) :: r1,r2,r3
  REAL(DP), INTENT(IN) :: theta1(nt1)
  REAL(DP), INTENT(IN) :: theta(nt)
  REAL(DP), INTENT(IN) :: phi(nphi)
  END SUBROUTINE pot_sym
END INTERFACE

INTERFACE test_potential
FUNCTION test_potential(r1,r2,r3, theta1,theta,phi)
  USE types
  USE param
  USE potential
! IMPLICIT NONE
  REAL(DP), INTENT(IN) :: r1,r2,r3
  REAL(DP), INTENT(IN) :: theta1(nt1)
  REAL(DP), INTENT(IN) :: theta(nt)
  REAL(DP), INTENT(IN) :: phi(nphi)
  LOGICAL(LGT) :: test_potential
  END FUNCTION test_potential
END INTERFACE

END MODULE

SUBROUTINE vexpand(vex,q1,q2,q3,testpot)
  USE types
  USE param
  USE base_lib, ONLY: plgndr, gammln, gauleg
!  USE potential
  USE workarrays, ONLY: v3
!  this is not required (recursion?)
! USE expansion, ONLY: test_potential, pot_sym
IMPLICIT NONE
  REAL(DP), INTENT(IN) :: q1(1:nn(1)), q2(1:nn(2)), q3(1:nn(3))
  REAL(DP), INTENT(OUT) :: vex(0:ne1,0:ne2,0:ne3,nn(1),nn(2),nn(3))
  LOGICAL(LGT), INTENT(OUT) :: testpot(nn(1),nn(2),nn(3))
REAL(DP), ALLOCATABLE :: x(:),wx(:),phi(:),tphi(:,:)
REAL(DP), ALLOCATABLE :: theta1(:),ttheta1(:,:,:),theta(:),ttheta(:,:,:)
REAL(DP) :: norm, norm0, sum1, sum2, sum3
INTEGER(I4B) :: n1,n2,n3, i1,i2,i3
INTEGER(I4B) :: i,j,l,k
INTEGER(I4B) :: fg1,fg2,fg3, rg1,rg2,rg3
INTEGER(I4B) :: is1,is2,is3
LOGICAL(LGT) :: test_potential
!  reset to zero threshold when doing expansion
REAL(DP), PARAMETER :: thresh=1.0e-6_dp

write(u6, '(/," expansion: start...")')

!  initialize expansion terms
vex = 0.0_dp

!  determine the grid sizes
!  reduced grid
rg1= nt1
rg2= nt
rg3= nphi
!  full grid
fg1= 2*rg1
fg2= 2*rg2
fg3= 2*rg3

!  then for i=1..rg1  x<0; for i=rg1+1..fg1  x>0 (0 < theta < Pi/2)
!  then for i=1..rg2  x<0; for i=rg2+1..fg2  x>0
!  then for i=1..rg3  x>0; for i=rg3+1..fg3  x<0

write(u6, '(" computing theta1 grid...")')

allocate(x(fg1),wx(fg1), &
  theta1(rg1),ttheta1(0:ne1,0:ne3,rg1))

call gauleg(-1.0_dp, 1.0_dp, x,wx)

theta1 = acos(x(rg1+1:fg1))
do l=0, ne1
  do k=0, min(ne3,l)
    norm = sqrt((2._dp*l+1._dp)/2._dp*exp(gammln(l-k+1._dp)-gammln(l+k+1._dp)))
    do i=1,rg1
      ttheta1(l,k,i)= norm*plgndr(l,k,x(i+rg1))*wx(i+rg1)
    end do
  end do
end do

deallocate(x,wx)

write(u6, '(" done.")')
write(u6, '(" computing theta grid...")')

allocate(x(fg2),wx(fg2), &
  theta(rg2),ttheta(0:ne2,0:ne3,rg2))

call gauleg(-1.0_dp, 1.0_dp, x,wx)

theta = acos(x(rg2+1:fg2))
do l=0, ne2
  do k=0, min(ne3,l)
    norm = sqrt((2._dp*l+1._dp)/2._dp*exp(gammln(l-k+1._dp)-gammln(l+k+1._dp)))
    do i=1,rg2
      ttheta(l,k,i)= norm*plgndr(l,k,x(i+rg2))*wx(i+rg2)
    end do
  end do
end do

deallocate(x,wx)

write(u6, '(" done.")')
write(u6, '(" computing phi grid...")')

allocate(phi(rg3),tphi(0:ne3,rg3))

do i=1,rg3
  phi(i)=PI*(i-0.5_dp)/real(fg3,dp)
end do

! use this if expnasion in normalised series is required
norm0 = sqrt(2.0_dp*PI)/real(fg3,dp)  !  norm INCLUDING uniform weights for k=0
norm  = 2.0_dp*sqrt(PI)/real(fg3,dp)  !  norm INCLUDING uniform weights for k>0

! use this if expnasion in non-normalised series is required
! i.e. natural expansion in cos(x)
! norm0 = 1.0_dp/real(nsz,dp)  !  norm INCLUDING uniform weights for k=0
! norm  = 2.0_dp/real(nsz,dp)  !  norm INCLUDING uniform weights for k>0

  tphi(0,:)= norm0
do k= 1, ne3
  if(mod(k,2)==0) then
    tphi( k,:)= norm*cos(k*phi(:))
  else
    tphi( k,:)= -norm*cos(k*phi(:))
  end if
end do

! above (-1)**k phase change has been included

write(u6, '(" done.")')

write(u6, '(/," computing expansion...")')

testpot = .FALSE.

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP& PRIVATE(n1,n2,n3,is1,is2,is3,sum1,sum2,sum3,i1,i2,i3)

allocate(v3(nt1,nt,nphi,0:1,0:1,0:1))

!$OMP DO SCHEDULE(DYNAMIC)

do n1=1,nn(1); do n2=1,nn(2); do n3=1,nn(3)

  if(test_potential(q1(n1),q2(n2),q3(n3), theta1,theta,phi)) then
    write(6,'(" Vmin > Encut1 for grid point",3i4)') n1,n2,n3
    testpot(n1,n2,n3) = .TRUE.
    cycle
  end if

  call pot_sym(q1(n1),q2(n2),q3(n3), theta1,theta,phi)

  do j=0,ne1, j_parity_max+1
    do l=0,ne2, l_parity_max+1
      if(jl_parity_max==1 .AND. mod(j+l,2) /= 0) cycle
      do k=0, min(l,j,ne3)

        is1=mod(j+k,2)
        is2=mod(l+k,2)
        is3=mod(k,2)

        sum3=0.0_dp
        do i3=1,rg3
    
          sum2=0.0_dp
          do i2=1,rg2

            sum1=0.0_dp
            do i1=1,rg1
              sum1=sum1+ v3(i1,i2,i3,is1,is2,is3)*ttheta1(j,k,i1)
            end do
            sum2=sum2+ ttheta(l,k,i2)*sum1

          end do
          sum3=sum3+ tphi(k,i3)*sum2

        end do

        if(abs(sum1) > thresh)  vex(j,l,k,n1,n2,n3)=sum3
        if(test_flag >1) write(u8,'(2x,3i3,f20.5)') j,l,k, sum3

      end do
    end do
  end do

end do; end do; end do

!$OMP END DO

deallocate(v3)

!$OMP END PARALLEL

deallocate(phi,tphi,theta,ttheta,theta1,ttheta1)

write(u6, '(" expansion: done.",/)')

END SUBROUTINE


! threej  is a wrapper for drc3jj which computes an
! array of 3-J Symbols
! j2min..j2max  is the desired range of the length 
!   j2size = j2max-j2min+1
! the full length is from abs(j1-j3)..(j1+j3) with length
!   j2sizebig = (j1+j3)+abs(j1-j3)+1 = 2*min(j1,j3)+1
!
FUNCTION threej(j1,j3,m1,m3,j2min,j2max,j2size,j2sb)
  USE types
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: j1,j3,m1,m3,j2min,j2max,j2size,j2sb
  REAL(DP) :: threej(j2size)
INTEGER(I4B) :: ier, j2min_allowed,j2max_allowed
REAL(DP) :: tp(j2sb)

j2min_allowed = max(abs(j1-j3),abs(m1+m3))
j2max_allowed = j1+j3

! write(6,'(5i5)') j1,j3, j2min, j2max, j2size

! all these IF's slow things down, uncomment if desired
! should be alright w/o them
!
! there are also some checks in drc3jj.f
! could remove it also

! if(j2size /= j2max-j2min+1) then
!   write(6,*) ' j2size is incorrect:', j2size
!   stop
! end if

! if(j2sb /= j2max_allowed-j2min_allowed+1) then
!   write(6,*) ' j2sb is incorrect:', j2sb
!   stop
! end if

! if(j2min<j2min_allowed) then
!   write(6,*) ' j2min < j2min_allowed', j2min,j2min_allowed
!   stop
! end if

! if(j2max>j2max_allowed) then
!   write(6,*) ' j2max > j2max_allowed', j2max,j2max_allowed
!   stop
! end if

call DRC3JJ(dble(j1), dble(j3), dble(m1), dble(m3), &
  dble(j2min_allowed), dble(j2max_allowed), tp, j2sb,IER)

if(ier>0) then
  write(6,*) '  ERROR in DRC3JJ: IER=', ier
  stop
end if

!do ier=1,j2sb
!  if(abs(tp(ier))>1.D9) then
!    write(6,*) ' tp is too big:', tp(ier)
!    write(6,*) j1,j3,m1,m3, j2min,j2max,j2sb
!    stop
!  end if
!end do

threej(1:j2size) = &
  tp(j2min-j2min_allowed+1:j2max-j2min_allowed+1)

!do ier=1,j2size
!  if(abs(threej(ier))>1.D9) then
!    write(6,*) ' threej is too big:', threej(ier)
!    stop
!  end if
!end do

END FUNCTION threej


!  setfac initialises threej0 array
!  which contains 3J symbols with M's = 0
!  multiplied by sqrt((2*i1+1)*(2*i2+1)*(2*i3+1)/2)
!
SUBROUTINE setfac()
  USE types
  USE param
  USE workarrays, ONLY: threej0
  !$OMP THREADPRIVATE(threej0)
IMPLICIT NONE

INTERFACE threej
FUNCTION threej(j1,j3,m1,m3,j2min,j2max,j2size,j2sizebig)
  USE types
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: j1,j3,m1,m3,j2min,j2max, &
    j2size, j2sizebig
  REAL(DP) :: threej(j2size)
  END FUNCTION threej
END INTERFACE

INTEGER(I4B) :: i1,i2,i2min,i2max,i3,i2size,i2sizebig
! REAL(DP) :: threej(nbin2+1)

write(6,*) '  setfac: starting...'

! initialise
threej0 = zero

do i1=0,nbin1
  do i3=0,nbin1
    ! need not to check whether i2min >= abs(m1+m3) = 0
    i2max = min(i3+i1,nbin2)
    i2min = abs(i3-i1)
    if(i2min>i2max) cycle
    i2size = i2max-i2min+1
!    i2sizebig = 2*min(i1,i3)+1
    i2sizebig = i3+i1-abs(i3-i1)+1

    threej0(i2min:i2max,i1,i3) = &
      threej(i1,i3,0,0,i2min,i2max,i2size,i2sizebig)

!    write(6,*) i1,i3
!    do i2=0,nbin2
!      write(6,*) i2,i1,i3, threej0(i2,i1,i3)
!    end do

    ! make sure all even sums are zero
    do i2=i2min,i2max
      ! if(mod(i1+i2+i3,2)==1)  threej0(i2,i1,i3) = zero
      threej0(i2,i1,i3) = sqrt((two*i1+1)*(i2+half)*(two*i3+1))* &
                          threej0(i2,i1,i3)
    end do

  end do
end do

write(6,*) '  setfac: done.'

END SUBROUTINE setfac

!  pot_sym returns symmetrical & asymmetrical parts of the potential
!  on reduced grid for theta1, theta and phi
!  it is used in numerical quadrature
!
SUBROUTINE pot_sym(r1,r2,r3, theta1,theta,phi)
  USE types
  USE param
  USE potential
  USE workarrays, ONLY: v3
  !$OMP THREADPRIVATE(v3)
IMPLICIT NONE
!  REAL(DP), INTENT(OUT) :: v3(nt1,nt,phi,0:1,0:1,0:1)
  REAL(DP), INTENT(IN) :: r1,r2,r3
  REAL(DP), INTENT(IN) :: theta1(nt1)
  REAL(DP), INTENT(IN) :: theta(nt)
  REAL(DP), INTENT(IN) :: phi(nphi)
REAL(DP) :: ppp,mpp,pmp,mmp,ppm,mpm,pmm,mmm
REAL(DP) :: pps,mps,pms,mms,ppa,mpa,pma,mma
REAL(DP) :: pss,mss,pas,mas,psa,msa,paa,maa
INTEGER(I4B) :: i1,i2,i3

!write(u8,'("  pot_sym: start...")')
  do i1=1, nt1; do i2=1, nt; do i3=1, nphi

    !  compute V on the full grid
    ppp = v(r1,r2,r3,theta1(i1),    theta(i2),    phi(i3))
    mpp = v(r1,r2,r3,PI-theta1(i1), theta(i2),    phi(i3))
    pmp = v(r1,r2,r3,theta1(i1),    PI-theta(i2), phi(i3))
    mmp = v(r1,r2,r3,PI-theta1(i1), PI-theta(i2), phi(i3))
    ppm = v(r1,r2,r3,theta1(i1),    theta(i2),    PI-phi(i3))
    mpm = v(r1,r2,r3,PI-theta1(i1), theta(i2),    PI-phi(i3))
    pmm = v(r1,r2,r3,theta1(i1),    PI-theta(i2), PI-phi(i3))
    mmm = v(r1,r2,r3,PI-theta1(i1), PI-theta(i2), PI-phi(i3))

    !  find V sym/asym in the respect to phi
    pps = ppp+ppm;   ppa = ppp-ppm
    pms = pmp+pmm;   pma = pmp-pmm
    mps = mpp+mpm;   mpa = mpp-mpm
    mms = mmp+mmm;   mma = mmp-mmm

    !  find V sym/asym in the respect to theta
    pss = pps+pms;  pas = pps-pms;  psa = ppa+pma;  paa = ppa-pma
    mss = mps+mms;  mas = mps-mms;  msa = mpa+mma;  maa = mpa-mma

    !  find V sym/asym in the respect to theta1 and store it
    v3(i1,i2,i3,0,0,0) = pss+mss;  v3(i1,i2,i3,1,0,0) = pss-mss
    v3(i1,i2,i3,0,1,0) = pas+mas;  v3(i1,i2,i3,1,1,0) = pas-mas
    v3(i1,i2,i3,0,0,1) = psa+msa;  v3(i1,i2,i3,1,0,1) = psa-msa
    v3(i1,i2,i3,0,1,1) = paa+maa;  v3(i1,i2,i3,1,1,1) = paa-maa

    !  symmetry index: 0 means sym; 1 means asym
    !  runs through theta1, theta and phi respectively

  end do; end do; end do

!write(u8,'("  pot_sym: done")')

END SUBROUTINE pot_sym


! function test_potential returns
! .TRUE.  if all angular grid points have E > Ecut1
! .FALSE. otherwise
!

FUNCTION test_potential(r1,r2,r3, theta1,theta,phi)
  USE types
  USE param
  USE potential
IMPLICIT NONE
  REAL(DP), INTENT(IN) :: r1,r2,r3
  REAL(DP), INTENT(IN) :: theta1(nt1)
  REAL(DP), INTENT(IN) :: theta(nt)
  REAL(DP), INTENT(IN) :: phi(nphi)
  LOGICAL(LGT) :: test_potential
REAL(DP) :: ppp,mpp,pmp,mmp,ppm,mpm,pmm,mmm
INTEGER(I4B) :: i1,i2,i3

! initialize
test_potential = .TRUE.

do i1=1, nt1; do i2=1, nt; do i3=1, nphi

    ppp = v(r1,r2,r3,theta1(i1),    theta(i2),    phi(i3))
    mpp = v(r1,r2,r3,PI-theta1(i1), theta(i2),    phi(i3))
    pmp = v(r1,r2,r3,theta1(i1),    PI-theta(i2), phi(i3))
    mmp = v(r1,r2,r3,PI-theta1(i1), PI-theta(i2), phi(i3))
    ppm = v(r1,r2,r3,theta1(i1),    theta(i2),    PI-phi(i3))
    mpm = v(r1,r2,r3,PI-theta1(i1), theta(i2),    PI-phi(i3))
    pmm = v(r1,r2,r3,theta1(i1),    PI-theta(i2), PI-phi(i3))
    mmm = v(r1,r2,r3,PI-theta1(i1), PI-theta(i2), PI-phi(i3))

    ! test if at least one grid point is less than Ecut1

    ! if(min(ppp,mpp,pmp,mmp,ppm,mpm,pmm,mmm) <= encut1) then
    if(min(ppp,mpp,pmp,mmp,ppm,mpm,pmm,mmm) <= margin1) then
      test_potential = .FALSE.
      return
    end if

end do; end do; end do

END FUNCTION
