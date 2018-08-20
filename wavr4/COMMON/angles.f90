! 10/09/02 correction to kp loopin angular_states_qnumbers
!          (moved up) to make better use of cache
!
! 11/08/02 cache matrix added in computing integrals in angle3d
!
! 01/08/02 phase change corrections in angleT3D
!
! 01/08/02 call to angleT3D in angle3d and angle3de is moved to main.f90
!          so that angular kinetic matrix can be saved and re-used
! 
! 30/07/02 adapted to rv57; minor corrections (sqrt(2.0_dp) replaced by 
!          SQRT2; norm factors IF construct)
!
! 06/07/02 corrections
!
! 3D angular problem:
! compute eigenvalues for a set of fixed radial coordinates
!
! quantum numbers convention:
! p -> p            parity (in the respect to inversion)
! J -> jr  K -> kr  rotational quantum numbers
! l -> l   k -> k   spherical harmonic q. numbers
! j -> j            assoc Legendre polynomial
!                   (m being determined from k and K)
! kappa -> kp       sign
!
! when parity of j & l is important (i.e. when xx_parity_max >0)
! j_parity          whether j is even or odd
! l_parity          whether l is even or odd
! jl_parity         whether j+l is even or odd
!
! r1,r2,r3 stretching grid point where the angular problem is
!          being computed
! ha, en   matrix and energy arrays (inout)
! e0       energy of the stretch, to be added to the diagonal
! nlev     number of levels returned with E < encut

MODULE angles

INTERFACE angular_states_max
SUBROUTINE angular_states_max()
  USE types
  USE param
  END SUBROUTINE
END INTERFACE

INTERFACE angular_states_qnumbers
SUBROUTINE angular_states_qnumbers(indx,krn)
  USE types
  USE param
  INTEGER(I4B), INTENT(OUT) :: indx(namax,6)
  TYPE(triple), INTENT(OUT) :: krn(0:krmax)
  END SUBROUTINE
END INTERFACE

INTERFACE angle3d
SUBROUTINE angle3d(r1,r2,r3,b1,b2,b3,ha,en,ts,nlev, &
    theta1,ttheta1, theta,ttheta, phi,tphi, indx)
  USE types
  USE param
  USE workarrays, ONLY: tp1,tp2,tp, v3
  USE expansion, ONLY: pot_sym
  USE potential
! IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: ha(na,na), en(na)
  REAL(DP), INTENT(IN)  :: r1,r2,r3,b1,b2,b3, ts
  REAL(DP), INTENT(IN) :: theta1(nt1),ttheta1(0:jmax,0:mmax,nt1)
  REAL(DP), INTENT(IN) :: theta(nt),   ttheta(0:lmax,0:kmax,nt)
  REAL(DP), INTENT(IN) :: phi(nphi),     tphi(0:2*kmax,nphi)
  INTEGER(I4B), INTENT(IN) :: indx(namax,6)
  INTEGER(I4B), INTENT(OUT) :: nlev
  END SUBROUTINE
END INTERFACE

INTERFACE angle3de
SUBROUTINE angle3de(r1,r2,r3,b1,b2,b3,ha,en,ts,nlev, vex3, indx)
  USE types
  USE param
  USE potential
! IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: ha(na,na), en(na)
  REAL(DP), INTENT(IN) :: r1,r2,r3,b1,b2,b3,ts, vex3(0:ne1,0:ne2,0:ne3)
  INTEGER(I4B), INTENT(IN) :: indx(namax,6)
  INTEGER(I4B), INTENT(OUT) :: nlev
  END SUBROUTINE
END INTERFACE

INTERFACE angle3d_diag
SUBROUTINE angle3d_diag(ha,en, r1,r2,r3, fd1r,fd2r,frr,nlev)
  USE types
  USE param
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEVR, LA_SYEV, LA_SYEVX
! IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: ha(na,na), en(na)
  REAL(DP), INTENT(IN) :: r1,r2,r3, fd1r,fd2r,frr
  INTEGER(I4B), INTENT(OUT) :: nlev
  END SUBROUTINE
END INTERFACE


INTERFACE angleT3D
SUBROUTINE angleT3D(ha,indx,i1start,i1size,istart,isize)
  USE types
  USE param
  INTEGER(I4B), INTENT(IN) :: indx(namax,6),i1start,i1size,istart,isize
  REAL(DP), INTENT(INOUT) :: ha(i1size,isize)
  END SUBROUTINE
END INTERFACE

END MODULE

SUBROUTINE angular_states_max()
  USE types
  USE param
IMPLICIT NONE
INTEGER(I4B) :: l,j,kp,kr,k, ktop, n, krtop

!if(test_flag >2) write(u8,'(/," angular_states_max",/,  &
! "    # kap: jr kr:  j  m:  l  k")')

krtop=min(jrmax,krmax)

namax=0
mmax=0

jr_loop: do jr=jrmax,jrmax
  p_loop: do p=0,1
    if(test_flag >1) write(u8,'(" p=",i3)') p

  j_parity_loop: do j_parity =0,j_parity_max
  l_parity_loop: do l_parity =0,l_parity_max
  jl_parity_loop: do jl_parity =0,jl_parity_max

n=0
do l=0,lmax                    ! l - orbital momentum
  do j=0,jmax                  ! j - another momentum

      if(j_parity_max==1 .AND. mod(j,2) /= j_parity) cycle
      if(l_parity_max==1 .AND. mod(l,2) /= l_parity) cycle
      if(jl_parity_max==1 .AND. mod(j+l,2) /= jl_parity) cycle


    do kp=-1,1,2               ! kappa =-1, 1 (i.e. B & A functions
      do kr=0,krtop            ! K - rotational number
        if (kr==0 .AND. kp==1) cycle
        ktop=min(l,kmax)
        do k=0,ktop            ! k - its projection

          ! if (abs(j-k) > 20 .OR. abs(l-k) > 20) cycle
          ! if (abs(j+l-2*k) > 20) cycle

          if (k==0 .AND. kp==1) cycle
          if (kr==0 .AND. k==0 .AND. (-1)**(jr+p) /= 1) cycle
          if (abs(k-kp*kr)>j) cycle
          n=n+1
          ! if(test_flag >2) write(u8,'(i5,i4,3(":",2i3))')  &
          !   n,kp,jr,kr,j,abs(k-kp*kr),l,k
          if(mmax<abs(k-kp*kr)) mmax=abs(k-kp*kr)   !  find mmax
        end do
      end do
    end do

  end do
end do

if(n>namax) namax = n


  end do jl_parity_loop
  end do l_parity_loop
  end do j_parity_loop

  end do p_loop
end do jr_loop

write(u6,*) ' angular_states_max: namax=',namax
write(u6,*) ' angular_states_max: mmax =',mmax

!  compute how many expansion terms are there for given ne1,ne2,ne3
!  Kr=0 => kp=-1

n=0
do j=0,ne1

  if(j_parity_max==1 .AND. mod(j,2) /= 0) cycle

  do l=0,ne2

    if(l_parity_max==1 .AND. mod(l,2) /= 0) cycle
    if(jl_parity_max==1 .AND. mod(j+l,2) /= 0) cycle

    ktop=min(l,j,ne3)
    n=n+(ktop+1)

  end do
end do

nexp = n
write(u6,*) ' angular_states_max: nexp =',nexp

END SUBROUTINE  angular_states_max


SUBROUTINE angular_states_qnumbers(indx,krn)
  USE types
  USE param
IMPLICIT NONE
  INTEGER(I4B), INTENT(OUT) :: indx(namax,6)
  TYPE(triple), INTENT(OUT) :: krn(0:krmax)
INTEGER(I4B) :: l,j,kp,kr,k, n ! ktop

if(test_flag >2) write(u8,'(/," angular_states:",/,  &
 & "    # kap: jr kr:  j  m:  l  k")')

indx(:,:) = 0
krn(:)%size = 0
n=0
! K goes first to make the basis in K-blocks
! k goes next to make use cache when computing integrals
do kr=0,krmax                    ! K - rotational number
  do k=0,kmax                    ! k - projection of l
    do kp=-1,1,2                 ! kappa =-1, 1 (i.e. B & A functions
      if (kr==0 .AND. kp==1) cycle
      if (k==0 .AND. kp==1) cycle

      do l=0,lmax                  ! l - orbital momentum assoc. with q2
        if(k > l) cycle
        do j=0,jmax                ! j - another momentum assoc. with q1

          if(j_parity_max==1 .AND. mod(j,2) /= j_parity) cycle
          if(l_parity_max==1 .AND. mod(l,2) /= l_parity) cycle
          if(jl_parity_max==1 .AND. mod(j+l,2) /= jl_parity) cycle

          ! if (abs(j-k) > 20 .OR. abs(l-k) > 20) cycle
          ! if (abs(j+l-2*k) > 20) cycle

!  the change of the line below was required because only jp is known
!  when the subroutine is being called
!          if (kr==0 .AND. k==0 .AND. (-1)**(jr+p) /= 1) cycle
          if (kr==0 .AND. k==0 .AND. jp == 1) cycle
          if (abs(k-kp*kr)>j) cycle
          n=n+1

          indx(n,1)=l
          indx(n,2)=j
          indx(n,3)=kp
          indx(n,4)=kr
          indx(n,5)=k
          indx(n,6)=abs(k-kp*kr)

          krn(kr)%size = krn(kr)%size +1

          if(test_flag >1) write(u8,'(i5,i4,3(":",2i3))')   &
             n, kp, jr,kr,j,abs(k-kp*kr),l,k

        end do

      end do
    end do
  end do
end do

na = n
write(u6,*) ' angular_states_qnumbers: na=',na
if(na>namax) then
  write(u6,*) ' angular_states_qnumbers: na>namax: stop'
  stop
end if

write(u7,'("  number of basis functions for various K",/,"    K    # start  end")')
do kr=0,krmax
  if(kr==0) then
    krn(kr)%start = 1
    krn(kr)%end = krn(kr)%size
  else
    krn(kr)%start = krn(kr-1)%end +1
    krn(kr)%end = krn(kr-1)%end +krn(kr)%size
  end if
  write(u7,'(4i5)') kr,krn(kr)
end do

END SUBROUTINE  angular_states_qnumbers

!  angle3d: standard version
!
! For historical reasons jr plays the role of kr here.
! Think of it as the sub-block of J matrix with K=J.
! kr value can be also determined from indx
! but kr1 and kr must always be the same
! Note: that jr is a global variable accessible via param
!       thereas kr is local

SUBROUTINE angle3d(r1,r2,r3,b1,b2,b3,ha,en,ts,nlev,  &
    theta1,ttheta1, theta,ttheta, phi,tphi, indx)
  USE types
  USE param
  USE workarrays, ONLY: tp1,tp2,tp, v3, krn
  USE expansion, ONLY: pot_sym
  !$OMP THREADPRIVATE(tp1,tp2,tp,v3)
IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: ha(na,na), en(na)
  REAL(DP), INTENT(IN) :: r1,r2,r3,b1,b2,b3,ts
  REAL(DP), INTENT(IN) :: theta1(nt1),ttheta1(0:jmax,0:mmax,nt1)
  REAL(DP), INTENT(IN) :: theta(nt),   ttheta(0:lmax,0:kmax,nt)
  REAL(DP), INTENT(IN) :: phi(nphi),     tphi(0:2*kmax,nphi)
  INTEGER(I4B), INTENT(IN) :: indx(namax,6)
  INTEGER(I4B), INTENT(OUT) :: nlev
REAL(DP) ::  sum3  ! sum1, sum2
REAL(DP) ::  norm1, norm, fct
REAL(DP), ALLOCATABLE :: v2(:,:,:,:), v1(:)
INTEGER(I4B) :: i1,i2,i3, dkp,dkm
INTEGER(I4B) :: n, j, l, k, kp, kr, m
INTEGER(I4B) :: n1,j1,l1,k1,kp1,kr1,m1
INTEGER(I4B) :: i, is1,is2,is3, ii, ii1
TYPE fournumbers
  INTEGER(I4B) :: k,k1,kp,kp1
END TYPE
TYPE(fournumbers) :: cache(0:1,0:1)

if(test_flag >2) then
  write(u8,'(" angle3d:...")')
end if

! this is sub-block with J=K, therefore
kr1= jr
kr = jr

allocate(v2(nt1,nt,0:1,0:1),v1(nt1))

!  stretch dependent coefficients
!  are now given in the input: b1, b2, b3
! call fdr(r1,r2,r3,frr,fd1r,fd2r)

!  kinetic matrix elements proportional to frr
!  01/08/02  it should be pre-computed in main.f90
! call angleT3D(ha,indx,krn(jr)%start,krn(jr)%size,krn(jr)%start,krn(jr)%size)

!  radial factor for angular kinetic matrix
ha = b3*ha

do ii=1,krn(jr)%size

  n =ii+krn(jr)%start-1

  l =indx(n,1)
  j =indx(n,2)

  ! Eqn 48 (the rest of)
  ! ha(ii,ii) = ha(ii,ii) +ts +0.5_dp*( fd1r*j*(j+1._dp)+fd2r*l*(l+1._dp) )
  ha(ii,ii) = ha(ii,ii) +ts +(b1+b3)*j*(j+one)+(b2+b3)*l*(l+one)
end do

call pot_sym(r1,r2,r3, theta1,theta,phi)

! write(u8,'("  angle3D: computing integrals...")')

do ii1=1,krn(jr)%size

  n1 =ii1+krn(jr)%start-1

  l1 =indx(n1,1)
  j1 =indx(n1,2)
  kp1=indx(n1,3)
  ! don't need kr1 here anymore because kr1=kr=jr
  ! kr1=indx(n1,4)
  k1 =indx(n1,5)
  m1 =indx(n1,6)

  if(kr1==0 .AND. k1==0) then
    norm1 = 0.5_dp
  else
    norm1 = 1._dp/SQRT2
  end if

  ! make sure cache is empty
  cache(:,:)%k1 = -999
  cache(:,:)%k  = -999
  cache(:,:)%kp = -999
  cache(:,:)%kp1= -999

  !  only Upper part (n=n1,na) as in the calling angle3D subroutine
  do ii =ii1,krn(jr)%size

    n =ii+krn(jr)%start-1

    l =indx(n,1)
    j =indx(n,2)
    kp=indx(n,3)
    ! don't need kr here anymore because kr1=kr=jr
    ! kr=indx(n,4)
    k =indx(n,5)
    m =indx(n,6)

    !  if K' =/= K then the integrals of V are zero
    ! if (kr1 /= kr) cycle

    if(kr==0 .AND. k==0) then
      norm = 0.5_dp
    else
      norm = 1._dp/SQRT2
    end if

    fct = norm1*norm
    if (mod(k1+k, 2) /= 0) fct = -fct

    is1=mod(j1+m1+j+m,2)
    is2=mod(l1+k1+l+k,2)
    is3=mod(k1+k,2)

    tp1(:) = ttheta1(j1,m1,:)*ttheta1(j,m,:)
    tp2(:) = ttheta(l1,k1,:)*ttheta(l,k,:)

    if(k1  == cache(is1,is2)%k1 .AND. &
       k   == cache(is1,is2)%k  .AND. &
       kp  == cache(is1,is2)%kp .AND. &
       kp1 == cache(is1,is2)%kp1) then

      ! then previous v2 can be re-used
    else
      cache(is1,is2)%k1 =k1
      cache(is1,is2)%k  =k
      cache(is1,is2)%kp =kp
      cache(is1,is2)%kp1=kp1

      !  this implements cases I & II and special via temporary array tp
      !  01/08/02: phase corrections: (-1)^(J+K+p) removed

      if (kr == 0) then                       ! special case K=0
        dkp = k1+k
        dkm = abs(k1-k)
        if (jp /= 0) then     !  (-1)**(J+p) factor
          do i=1,nphi; tp(i) = tphi(dkm,i)-tphi(dkp,i); end do
        else
          do i=1,nphi; tp(i) = tphi(dkm,i)+tphi(dkp,i); end do
        end if
      else if (kr /= 0 .AND. kp1 == kp) then  ! case I
        dkm = abs(k1-k)
          do i=1,nphi; tp(i) = tphi(dkm,i); end do
      else                                    ! case II
        dkp = k1+k
        ! if (mod(jp+kr, 2) /= 0) then     !  (-1)**(J+K+p) factor removed
        !   do i=1,nphi; tp(i) = -tphi(dkp,i); end do
        ! else
          do i=1,nphi; tp(i) =  tphi(dkp,i); end do
        ! end if
      end if

      v2(:,:,is1,is2)= zero
      do i3=1,nphi
        do i2=1,nt
          do i1=1,nt1
             v2(i1,i2,is1,is2)= &
               v2(i1,i2,is1,is2)+tp(i3)*v3(i1,i2,i3,is1,is2,is3)
          end do
        end do
      end do
    end if

    ! matmul(v2,tp2)
    call dgemv('N',nt1,nt,one,v2(:,:,is1,is2),nt1,tp2,1,zero,v1(:),1)

    ! dot_product(v1(:),tp1(:))
    sum3= dot_product(v1(:),tp1(:))

    sum3 = 2._dp*fct*sum3
    !  check the integrals
    if(test_flag >2) write(u8,'(2x,4i4,f20.6,8i3)') &
                     ii1,ii,kr,kp1*kp, sum3,j1,m1,l1,k1,j,m,l,k
    ha(ii1,ii) = ha(ii1,ii) + sum3

    ! write(u8,'(2x,2i5,f20.10)') n1,n,ha(ii1,ii)

  end do
end do

deallocate(v2,v1)

call angle3d_diag(ha,en, r1,r2,r3, b1,b2,b3,nlev)

if(test_flag >2) then
  write(u8,'(" angle3d: done")')
end if

END SUBROUTINE angle3d


!  angle3de: angle3d version for expanded potential
!
! For historical reasons jr plays the role of kr here.
! Think of it as the sub-block of J matrix with K=J.
! kr value can be also determined from indx
! but kr1 and kr must always be the same
! Note: that jr is a global variable accessible via param
!       thereas kr is local

SUBROUTINE angle3de(r1,r2,r3,b1,b2,b3,ha,en,ts,nlev, vex3, indx)
  USE types
  USE param
  USE workarrays, ONLY: threej0, tp1, tp2, krn
  USE expansion, ONLY: threej
  !$OMP THREADPRIVATE(tp1,tp2,threej0)
IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: ha(na,na), en(na)
  REAL(DP), INTENT(IN) :: r1,r2,r3,b1,b2,b3,ts, vex3(0:ne1,0:ne2,0:ne3)
  INTEGER(I4B), INTENT(IN) :: indx(namax,6)
  INTEGER(I4B), INTENT(OUT) :: nlev
REAL(DP) :: norm1,norm, summ,sum1,sum2,sfac, one_over_sqrt2
INTEGER(I4B) :: j2,l2,k2, ii,ii1 !  i
INTEGER(I4B) :: n, j, l, k, kp, kr, m, m_abs
INTEGER(I4B) :: n1,j1,l1,k1,kp1,kr1,m1,m1_abs
INTEGER(I4B) :: j2min,j2max,l2min,l2max, j2sb,l2sb

one_over_sqrt2 = one/sqrt2

if(test_flag >2) then
  write(u8,'(" angle3de:...")')
end if

! this is sub-block with J=K, therefore
kr1= jr
kr = jr

! this is just a test
!do j=0,nbin1
!  do j1=0,nbin1
!    do j2=0,nbin2
!      if(threej0(j2,j1,j) ==zero) 
!write(6,*) j2,j1,j,threej0(j2,j1,j)
!    end do
!  end do
!end do

!  diagnostic block
! print *, ' ttheta size', size(ttheta,1), lbound(ttheta,1), ubound(ttheta,1)
!  end diagnostic block

!  stretch dependent coefficients
!  are now given in the input: b1, b2, b3
! call fdr(r1,r2,r3,frr,fd1r,fd2r)

!  kinetic matrix elements proportional to frr
!  01/08/02  it should be pre-computed in main.f90
! call angleT3D(ha,indx,krn(jr)%start,krn(jr)%size,krn(jr)%start,krn(jr)%size)

!  stretching dependence in the kinetic matrix
ha = b3*ha

do ii=1,krn(jr)%size

  n =ii+krn(jr)%start-1

  l =indx(n,1)
  j =indx(n,2)

  ! Eqn 48 (the rest of)
  ! ha(ii,ii) = ha(ii,ii) +ts +0.5_dp*( fd1r*j*(j+1._dp)+fd2r*l*(l+1._dp) )
  ha(ii,ii) = ha(ii,ii) +ts +(b1+b3)*j*(j+1._dp)+(b2+b3)*l*(l+1._dp)
end do

! write(u8,'("  angle3de: computing integrals...")')

do ii1=1,krn(jr)%size

  n1 =ii1+krn(jr)%start-1

  l1 =indx(n1,1)
  j1 =indx(n1,2)
  kp1=indx(n1,3)
  ! don't need kr1 here anymore because kr1=kr=jr
  ! kr1=indx(n1,4)
  k1 =indx(n1,5)
  m1_abs =indx(n1,6)

  if(kr1==0 .AND. k1==0) then
    norm1 = half
  else
    norm1 = one_over_SQRT2
  end if

  !  only Upper part (n=n1,na) as in the calling angle3d subroutine
  do ii =ii1,krn(jr)%size

    n =ii+krn(jr)%start-1

    l =indx(n,1)
    j =indx(n,2)
    kp=indx(n,3)
    ! don't need kr here anymore because kr1=kr=jr
    ! kr=indx(n,4)
    k =indx(n,5)
    m_abs =indx(n,6)

    !  if K' =/= K then the integrals of V are zero
    ! if (kr1 /= kr) cycle

    if(kr==0 .AND. k==0) then
      norm = half
    else
      norm = one_over_SQRT2
    end if

    !  compute integral: PART A & B
    !  the case K=0 is a sum of A + B

    j2max = j1+j
    call sameparity_minus(j2max,ne1)
    ! if(abs(j1-j)>j2max) cycle

    l2max = l1+l
    call sameparity_minus(l2max,ne2)
    ! if(abs(l1-l)>l2max) cycle

    summ=0.0_dp

    ! we do a lot of checks below (which slows down everything)
    ! to comply with 3J rules and be within the expansion limits
    ! j2sb & l2sb are full lengths of possible j2 and l2 respecitively
    ! without the restrictions implied by the expansion limits

    !  PART A
    !

    k2 = abs(k1-k)
    j2min = abs(j1-j)
    l2min = abs(l1-l)
    call sameparity_plus(j2min,k2)
    call sameparity_plus(l2min,k2)

    if(abs(k2)<=ne3 .AND. j2min<=j2max .AND. l2min<=l2max .AND. &
       (kr == 0 .OR. (kr /=0 .AND. kp1==kp)) ) then

        j2sb  = j1+j-max(abs(j1-j),k2)+1
        l2sb  = l1+l-max(abs(l1-l),k2)+1

        select case(k1-k)
        case(0) ! k1=k
          !  k2 = 0

          m1 = -(k1-kp*kr)
          m  =  (k -kp*kr)

          sfac = two*(-1)**(k +(m1-m1_abs+m-m_abs)/2)*one_over_sqrt2pi
          tp1(j2min:j2max)=threej(j1,j,m1,m,j2min,j2max,j2max-j2min+1,j2sb)
          tp2(l2min:l2max)=threej(l1,l,-k,k,l2min,l2max,l2max-l2min+1,l2sb)

        case(1:) ! k1>k

          m1 = -(k1-kp*kr)
          m  =  (k -kp*kr)

          sfac = (-1)**(k1 +(m1-m1_abs+m-m_abs)/2)*one_over_sqrtpi
          tp1(j2min:j2max)=threej(j1,j, m1,m,j2min,j2max,j2max-j2min+1,j2sb)
          tp2(l2min:l2max)=threej(l1,l,-k1,k,l2min,l2max,l2max-l2min+1,l2sb)

        case(:-1) ! k>k1
          ! k2 = -k2   ! k2 is positive anyway

          m1 =  (k1-kp*kr)
          m  = -(k -kp*kr)

          sfac = (-1)**(k +(m1-m1_abs+m-m_abs)/2)*one_over_sqrtpi
          tp1(j2min:j2max)=threej(j1,j,m1, m,j2min,j2max,j2max-j2min+1,j2sb)
          tp2(l2min:l2max)=threej(l1,l,k1,-k,l2min,l2max,l2max-l2min+1,l2sb)

        end select

        ! write(6,'(" A:",7i5)') j2min,j2max, l2min,l2max, m1,m, k2

        sum1 = zero
        do l2 =l2min,l2max,2
          sum2 = zero
          do j2 =j2min,j2max,2
            ! if(vex3(j2,l2,k2)==0.0_dp) cycle
            ! write(6,*) j2,l2,k2,vex3(j2,l2,k2),threej0(j2,j1,j),tp1(j2)
            sum2 = sum2+vex3(j2,l2,k2)*threej0(j2,j1,j)*tp1(j2)
          end do
          sum1 = sum1+ sum2*threej0(l2,l1,l)*tp2(l2)
        end do

        !if(abs(sum1)==0.D0) then
        !  write(6,'(" A:",7i5)') j2min,j2max, l2min,l2max, m1,m, k2
        !  write(6,'("  sum1",f20.9)') sum1
        !  write(6,*) ' 3j0 j:',threej0(j2min:j2max,j1,j)
        !  write(6,*) ' tp1:',tp1(j2min:j2max)
        !  write(6,*) ' 3j0 l:',threej0(l2min:l2max,l1,l)
        !  write(6,*) ' tp2:',tp2(l2min:l2max)
        !  write(6,*) ' ',threej(j1,j,m1, m,j2min,j2max,j2max-j2min+1,j2sb)
        !  ! stop
        !end if

        summ = summ+sfac*sum1

    end if

    !  PART B
    !

    k2 = k1+k
    j2min = abs(j1-j)
    l2min = abs(l1-l)
    call sameparity_plus(j2min,k2)
    call sameparity_plus(l2min,k2)

    if(abs(k2)<=ne3 .AND. j2min<=j2max .AND. l2min<=l2max .AND. &
       (kr == 0 .OR. (kr /=0 .AND. kp1==-kp)) ) then

        j2sb  = j1+j-max(abs(j1-j),k2)+1
        l2sb  = l1+l-max(abs(l1-l),k2)+1

        if(k2==0) then
          !  k2 = 0, k1 = 0 & k = 0
          !  and we could use threej0 instead of computing tp2

          m1 = -kp*kr
          m  =  kp*kr

          ! 02/08/02; 21/02/03 phase correction change
          ! sfac = two*(-1)**(jp+kr +(m1-m1_abs+m-m_abs)/2)*one_over_sqrt2pi
          sfac = two*(-1)**((m1-m1_abs+m-m_abs)/2)*one_over_sqrt2pi
          if(kr == 0) sfac = sfac*(-1)**(jp+kr)
          tp1(j2min:j2max)=threej(j1,j,m1,m,j2min,j2max,j2max-j2min+1,j2sb)
          tp2(l2min:l2max)=threej(l1,l, 0,0,l2min,l2max,l2max-l2min+1,l2sb)

        else ! if (k2 /= 0) then
          !  k2 = k1+k

          m1 = -(k1+kp*kr)
          m  = -(k -kp*kr)

          ! 02/08/02; 21/02/03 phase correction change
          ! sfac = (-1)**(jp+kr+k1+k +(m1-m1_abs+m-m_abs)/2)*one_over_sqrtpi
          sfac = (-1)**(k1+k +(m1-m1_abs+m-m_abs)/2)*one_over_sqrtpi
          if(kr == 0) sfac = sfac*(-1)**(jp+kr)
          tp1(j2min:j2max)=threej(j1,j, m1, m,j2min,j2max,j2max-j2min+1,j2sb)
          tp2(l2min:l2max)=threej(l1,l,-k1,-k,l2min,l2max,l2max-l2min+1,l2sb)

        end if

        ! write(6,'(" B:",7i5)') j2min,j2max, l2min,l2max, m1,m, k2

        sum1 = zero
        do l2 =l2min,l2max,2
          sum2 = zero
          do j2 =j2min,j2max,2
            ! if(vex3(j2,l2,k2)==0.0_dp) cycle
            ! write(6,*) j2,l2,k2,vex3(j2,l2,k2),threej0(j2,j1,j),tp1(j2)
            sum2 = sum2+vex3(j2,l2,k2)*threej0(j2,j1,j)*tp1(j2)
          end do
          sum1 = sum1+ sum2*threej0(l2,l1,l)*tp2(l2)
        end do

        !if(abs(sum1)==0.D0) then
        !  write(6,'(" B:",7i5)') j2min,j2max, l2min,l2max, m1,m, k2
        !  write(6,'("  %",f20.9)') sum1
        !  write(6,*) ' 3j0 j:',threej0(j2min:j2max,j1,j)
        !  write(6,*) ' tp1:',tp1(j2min:j2max)
        !  write(6,*) ' 3j0 l:',threej0(l2min:l2max,l1,l)
        !  write(6,*) ' tp2:',tp2(l2min:l2max)
        !  write(6,*) ' ',threej(j1,j,m1, m,j2min,j2max,j2max-j2min+1,j2sb)
        !  stop
        !end if

        summ = summ+sfac*sum1

    end if

    summ = norm1*norm*summ
    !  check the integrals
    if(test_flag >2) write(u8,'(2x,4i4,f20.6,8i3)') &
      n1,n,kr,kp1*kp, summ,j1,m1_abs,l1,k1,j,m_abs,l,k
    ha(ii1,ii) = ha(ii1,ii) + summ

    ! write(u8,'(2x,2i5,f20.10)') n1,n,ha(ii1,ii)

  end do
end do

call angle3d_diag(ha,en, r1,r2,r3, b1,b2,b3,nlev)

if(test_flag >2) then
  write(u8,'(" angle3de: done")')
end if

CONTAINS

! sameparity_plus  returns ii not less than kk and ii
! with the same parity (even or odd) as in the input
!
SUBROUTINE sameparity_plus(ii,kk)
  USE types
IMPLICIT NONE
  INTEGER(I4B), INTENT(INOUT) :: ii
  INTEGER(I4B), INTENT(IN) :: kk
if(ii<kk) then
  if(mod(ii+kk,2)==0) then
    ii=kk
  else
    ii=kk+1
  end if
end if

END SUBROUTINE

! sameparity_minus returns ii not greater than kk and ii
! with the same parity (even or odd) as in the input
!
SUBROUTINE sameparity_minus(ii,kk)
  USE types
IMPLICIT NONE
  INTEGER(I4B), INTENT(INOUT) :: ii
  INTEGER(I4B), INTENT(IN) :: kk
if(ii>kk) then
  if(mod(ii+kk,2)==0) then
    ii=kk
  else
    ii=kk-1
  end if
end if

END SUBROUTINE

END SUBROUTINE angle3de


!  diagonalizer common for angle3d & angle3de
!  returns eigenvalues and eigenvectors

SUBROUTINE angle3d_diag(ha,en, r1,r2,r3, fd1r,fd2r,frr,nlev)
  USE types
  USE param
  USE LA_PRECISION, ONLY: WP => DP
  USE F95_LAPACK, ONLY: LA_SYEVR, LA_SYEV, LA_SYEVX
IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: ha(na,na), en(na)
  REAL(DP), INTENT(IN) :: r1,r2,r3, fd1r,fd2r,frr
  INTEGER(I4B), INTENT(OUT) :: nlev
INTEGER(I4B) i, idiag

!  if all levels will be computed, return nlev=na
nlev = na

if(test_flag >2) write(u8,'(" angle3d_diag: diagonalizing...")')

  select case (icut1)
  case(0)
    ! find all levels
    call la_syev(ha,en,JOBZ='V',UPLO='U',INFO=idiag)
  case(1:)
    ! find icut1 levels
    call la_syevr(ha,en,JOBZ='V',UPLO='U',IU=min(icut1,nlev),M=nlev,INFO=idiag)
  case(:-1)
    ! find levels E < encut1
    ! Somehow syevx made crashes on win32 and wrong results on Tru64
    ! with certain input parameters. The reason is unclear but syevr works OK
    ! call la_syevx(ha,en,JOBZ='V',UPLO='U',VU=encut1,M=nlev,INFO=idiag)
    call la_syevr(ha,en,JOBZ='V',UPLO='U',VU=encut1,M=nlev,INFO=idiag)
  end select

if(angular_problem_only_flag == 1) then
  write(u8,'(//,"  computing J=",i2,"  p=",i2)') jr,p
  !   ALL LEVELS
  ! call la_syev(ha,en,JOBZ='N',UPLO='U',INFO=idiag)
  !
  !   A FEW lowest by energy
  ! call la_syevx(ha,en,JOBZ='N',UPLO='U',VU=encut1,M=nlev,INFO=idiag)
  !     ABSTOL=2*LA_LAMCH(1.0_dp,'Safe minimum'))
  !  A FEW lowest by number
  ! call la_syevr(ha,en,JOBZ='N',UPLO='U',IU=min(icut1,na),M=nlev,INFO=idiag)
  ! call la_syevx(ha,en,JOBZ='N',UPLO='U',IU=min(icut1,na),M=nlev,INFO=idiag)
  !   ABSTOL=LA_LAMCH(1.0_dp,'Safe minimum'))

  !  test print out
  write(u8,'("   r1=",f16.10,"   r2=",f16.10,"  r3=",f16.10)') r1,r2,r3
  write(u8,'(" fd1r=",f16.10," fd2r=",f16.10," frr=",f16.10)') fd1r,fd2r,frr
  write(u8,'(" 2B1 =",f16.10," 2B2 =",f16.10," 2B3=",f16.10)') fd1r-2*frr,fd2r-2*frr,2*frr
  write(u8,*) ' ha energies:'
  write(u8,'("   E0=",f26.6)') en(1)
  do i=1,min(na,nlev)
!    write(u8,'(2x,i4,f26.5)') i,en(i)-en(1)
    write(u8,'(2x,i4,f26.5)') i,en(i)
  end do
end if

if(idiag /= 0) then
  write(u6,'("   angle3d_diag: failed to diagonalise, idiag=",i4)') idiag
  stop
end if

if(test_flag >2) write(u8,'(" angle3d_diag: diagonalization done")')

END SUBROUTINE angle3d_diag

! radial- independent angular kinetic matrix
! generalized to return a spesified block
! (i1start:i1start+n1size)x(istart:istart+nsize)
!
! 01/08/02: phase change corrections:
! we introduce the new basis = c(-1)^(J+K+p)*[old basis kappa=1]. This 
! does affect the matrix elements where kappa and K do not change:
! Eqs: 50, 51, 52-55 and potential matrix.
! All these matrix element are multiplied by (-1)**(jp+kr).
! If it was already there then it was removed. If it was not then 
! the necessary factor inserted. The reason to change the phase is that
! this allows to re-use the eigenvectors for K>0 p=0 for p=1

SUBROUTINE angleT3D(ha,indx,i1start,i1size,istart,isize)
  USE types
  USE param
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: indx(namax,6),i1start,i1size,istart,isize
  REAL(DP), INTENT(INOUT) :: ha(i1size,isize)
  REAL(DP) :: fct
INTEGER(I4B) :: kp ,k, j, l, kr
INTEGER(I4B) :: kp1,k1,j1,l1,kr1
INTEGER(I4B) :: i,i1, n,n1

ha(1:i1size,1:isize) = 0.0_dp

do i1=1,i1size

  n1 =i1+i1start-1

  ! print *, '  n1=',n1

  l1 =indx(n1,1)
  j1 =indx(n1,2)
  kp1=indx(n1,3)
  kr1=indx(n1,4)
  k1 =indx(n1,5)

!  only Upper (i=i1,na) or Lower(i=1,i1) triangular part is sufficient
  do i=1,isize

    n =i+istart-1

    ! diagonal contribution when quantum numbers are the same, ie
    if(n1 == n) then
      ! Eqn 48
      ! ha(i1,i1) = jr*(jr+1._dp)-2._dp*(kr1*kr1+k1*k1-kp1*kr1*k1)
      !             jr*(jr+1._dp) has been moved to main program as "eshift"
      ha(i1,i1) = -two*(kr1*kr1+k1*k1-kp1*kr1*k1)
    end if

    l =indx(n,1)
    j =indx(n,2)
    kp=indx(n,3)
    kr=indx(n,4)
    k =indx(n,5)

    !  kintetic matrix is block diagonal in l,j
    if(l1==l .AND. j1==j) then

      if(kp1==kp) then

        ! Eqs 49-51

        if (kr1==kr .AND. k1==k+1) then
!          write(u8,'(" Eq. 49:L")')
          if(kr==0 .AND. k==0) then
            fct= sqrt2
          else
            fct= one
          end if
          ha(i1,i)=ha(i1,i)+ &
           sign(1,k-kp*kr)*fct*lf(j,k-kp*kr,+1)*lf(l,k,+1)
        end if

        if (kr1==kr .AND. k1==k-1) then
!          write(u8,'(" Eq. 49:U")')
          if(kr==0 .AND. k==1) then
            fct= sqrt2
          else
            fct= one
          end if
          ha(i1,i)=ha(i1,i)+ &
!           sign(1,k-kp*kr)*fct*lf(j,k-kp*kr,-1)*lf(l,k,-1)
! this formular was obtained from Eq 49:L
            sign(1,k1-kp*kr)*fct*lf(j,k-kp*kr,-1)*lf(l,k,-1)
        end if

        if (kr1==kr-1 .AND. k1==k-kp) then
!          write(u8,'(" Eq. 50:Ua")')
          fct= one;  if(kp == -1) fct=-fct
          ha(i1,i)=ha(i1,i) &
           -lf(jr,kr,-1)*lf(l,k,-sign(1,kp))*fct
        end if

        if (kr1==kr-1 .AND. k1==k) then
!          write(u8,'(" Eq. 50:Ub")')
          if(kr==1 .AND. k==0) then
            fct= sqrt2
          else
            fct= one
          end if
          if(kp == -1) fct=-fct
          ha(i1,i)=ha(i1,i) &
           -sign(1,k-kp*kr)*fct*lf(jr,kr,-1)*lf(j,k-kp*kr,sign(1,kp))
        end if

        if (kr1==kr+1 .AND. k1==k+kp) then
!          write(u8,'(" Eq. 51:La")')
          fct= one;  if(kp == -1) fct=-fct
          ha(i1,i)=ha(i1,i) &
           -lf(jr,kr,+1)*lf(l,k,sign(1,kp))*fct
        end if

        if (kr1==kr+1 .AND. k1==k) then
!          write(u8,'(" Eq. 51:Lb")')
          if(kr==0 .AND. k==0) then
            fct= sqrt2
          else
            fct= one
          end if
          if(kp == -1) fct=-fct
          ha(i1,i)=ha(i1,i) &
!           -sign(1,k-kp*kr)*fct*lf(jr,kr,+1)*lf(j,k-kp*kr,-sign(1,kp))
! this formular was obtained from Eq. 50b
           -sign(1,k-kp*kr1)*fct*lf(jr,kr,+1)*lf(j,k-kp*kr,-sign(1,kp))
        end if

      end if

      ! Eqs 52-55 (Lower Triangle since kp runs first

      if (kp1==1 .AND. kp==-1) then

        if (kr1==kr .AND. k1==1 .AND. k==0) then
!          write(u8,'(" Eq. 52:L")')
!          ha(i1,i)=ha(i1,i)- (-1)**(jp+kr)*lf(j,kr,-1)*lf(l,0,-1)
          ha(i1,i)=ha(i1,i)- lf(j,kr,-1)*lf(l,0,-1)
        end if

        if (kr1==1 .AND. kr==0 .AND. k1==k) then
!          write(u8,'(" Eq. 53:L")')
!          ha(i1,i)=ha(i1,i)- lf(jr,0,+1)*lf(j,k,-1)
          ha(i1,i)=ha(i1,i)- (-1)**(jp+kr)*lf(jr,0,+1)*lf(j,k,-1)
        end if

        if (kr1==1 .AND. kr==0 .AND. k1==k+1) then
!          write(u8,'(" Eq. 54:L")')
          if(k==0) then
            fct= sqrt2
          else
            fct= one
          end if
!          ha(i1,i)=ha(i1,i)- fct*lf(jr,0,+1)*lf(l,k,+1)
          ha(i1,i)=ha(i1,i)- (-1)**(jp+kr)*fct*lf(jr,0,+1)*lf(l,k,+1)
        end if

        if (kr1==kr+1 .AND. k1==1 .AND. k==0 .AND. kr /= 0) then
!          write(u8,'(" Eq. 55:L")')
          if(kr==0) then
            fct= sqrt2
          else
            fct= one
          end if
!          ha(i1,i)=ha(i1,i)- (-1)**(jp+kr)*fct*lf(jr,kr,+1)*lf(l,0,+1)
          ha(i1,i)=ha(i1,i)- fct*lf(jr,kr,+1)*lf(l,0,+1)
        end if
      end if

      ! Eqs 52-55 modified for Upper Triangle (since kp runs first)

      if (kp1==-1 .AND. kp==1) then

        if (kr1==kr .AND. k1==0 .AND. k==1) then
!          write(u8,'(" Eq. 52:U")')
!          ha(i1,i)=ha(i1,i)- (-1)**(jp+kr)*lf(j,kr,-1)*lf(l,0,-1)
          ha(i1,i)=ha(i1,i)- lf(j,kr,-1)*lf(l,0,-1)
        end if

        if (kr1==0 .AND. kr==1 .AND. k1==k) then
!          write(u8,'(" Eq. 53:U")')
!          ha(i1,i)=ha(i1,i)- lf(jr,0,+1)*lf(j,k,-1)
          ha(i1,i)=ha(i1,i)- (-1)**(jp+kr1)*lf(jr,0,+1)*lf(j,k,-1)
        end if

        if (kr1==0 .AND. kr==1 .AND. k1==k-1) then
!          write(u8,'(" Eq. 54:U")')
          if(k1==0) then   ! OR equivalently k==1
            fct= sqrt2
          else
            fct= one
          end if
!          ha(i1,i)=ha(i1,i)- fct*lf(jr,0,+1)*lf(l,k1,+1)
          ha(i1,i)=ha(i1,i)- (-1)**(jp+kr1)*fct*lf(jr,0,+1)*lf(l,k1,+1)
        end if

        if (kr1==kr-1 .AND. k1==0 .AND. k==1 .AND. kr1 /= 0) then
!          write(u8,'(" Eq. 55:U")')
          if(kr1==0) then   ! OR equivalently kr==1
            fct= sqrt2
          else
            fct= one
          end if
!          ha(i1,i)=ha(i1,i)- (-1)**(jp+kr1)*fct*lf(jr,kr1,+1)*lf(l,0,+1)
          ha(i1,i)=ha(i1,i)- fct*lf(jr,kr1,+1)*lf(l,0,+1)
        end if
      end if

    end if

  end do
end do

CONTAINS

FUNCTION lf(l,k,s)
USE types
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: l,k,s
REAL(DP) :: lf
lf = sqrt(l*(l+1._dp)-k*(k+s))
END FUNCTION

END SUBROUTINE angleT3D
