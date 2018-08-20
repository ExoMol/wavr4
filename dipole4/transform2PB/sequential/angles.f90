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
  USE workarrays, ONLY: angsym  
! INTEGER(I4B), INTENT(OUT) :: angsym(0:krmax,namax)
  INTEGER(I4B), INTENT(OUT) :: indx(namax,6)
  TYPE(triple), INTENT(OUT) :: krn(0:krmax)
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
  USE workarrays, ONLY: angsym  
! INTEGER(I4B), INTENT(OUT) :: angsym(0:krmax,namax)
IMPLICIT NONE
  INTEGER(I4B), INTENT(OUT) :: indx(namax,6)
  TYPE(triple), INTENT(OUT) :: krn(0:krmax)
INTEGER(I4B) :: l ,j, kp, kr ,k, n, m, i ! ktop
INTEGER(I4B) :: l1,j1,kp1,kr1,k1,n1,m1,i1

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

! now form the symmetry map for the angular basis 
! in the respect to H-H permutation

do n=1,na
  l =indx(n,1)
  j =indx(n,2)
  kp=indx(n,3)
  kr=indx(n,4)
  k =indx(n,5)
  m =indx(n,6)

  ! i = n - krn(kr)%start + 1

  ! restriction: K=0 => the rules are simple:  l <-> j
  ! find n1 such that l = j1 and j = l1
  do n1=1,na
    l1 =indx(n1,1)
    j1 =indx(n1,2)
    kp1=indx(n1,3)
    kr1=indx(n1,4)
    k1 =indx(n1,5)
    m1 =indx(n1,6)

    ! i1 = n1 - krn(kr1)%start + 1

    if(l1==j .and. j1==l .and. kp1==kp .and. kr1==kr .and. k1==k .and. m1==m) then
      angsym(kr,n) = n1
!      write(6,*) n,' - ',n1
      exit
    end if

  end do

end do

END SUBROUTINE  angular_states_qnumbers
