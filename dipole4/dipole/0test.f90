! a stand along test to check quadrature over Phi

PROGRAM threejtest
  USE types
  IMPLICIT NONE
  INTEGER(4) :: narg, kmax, k1,k2,n, nphi
  REAL(8), ALLOCATABLE :: phi(:),tphi_cos(:,:),tphi_sin(:,:)
  REAL(8) :: cc, ss, cs
  real(8), parameter :: eps= 0.00001D0
  CHARACTER*10 :: arg

! input processing
narg=iargc()

if (narg < 2 ) then
  write(6,*) " usage: 0test kmax nphi"
  write(6,*) "        the general rule for nphi=kmax/2+1"
  stop
end if

call getarg(1,arg)
read(arg,'(i10)') kmax

call getarg(2,arg)
read(arg,'(i10)') nphi

allocate(phi(nphi),tphi_cos(0:2*kmax,nphi),tphi_sin(0:2*kmax,nphi))

! choose either ... or ...
! the grid points normalised for computing f(x_i)*cos(k x_i)*cos(k' x_i)
! ie sum tphi_cos * tphi_cos = 1
!call grid_phi2(phi,tphi_cos,tphi_sin,kmax,nphi)
! the grid points normalised for computing f(x_i)*cos(k x_i)
! ie sum cos * tphi_cos = 1
call grid_phi(phi,tphi_cos,kmax,nphi)
call grid_phi_sin(phi,tphi_sin,kmax,nphi)

write(6,'("integrals greater than EPS=",f21.16)'), eps
write(6,*) '   k1   k2  cos(k1 x) cos(k2 x)  sin(k1 x) sin(k2 x)'
do k1=0,kmax
  !write(6,*) '   k1   k2  cos(k1 x) cos(k2 x)  sin(k1 x) sin(k2 x)'
  do k2=0,kmax
    cc=zero; ss=zero; cs=zero
    do n=1,nphi

      !cc = cc+tphi_cos(k1,n)*tphi_cos(k2,n)
      cc = cc+tphi_cos(k1,n)*cos(k2*phi(n))

      !ss = ss+tphi_sin(k1,n)*tphi_sin(k2,n)
      ss = ss+tphi_sin(k1,n)*sin(k2*phi(n))

    enddo
    cc=cc*4; ss=ss*4
    if (k2==0) then
      cc=cc/SQRT2; ss=ss/SQRT2
    end if
    if (k1==0) then
      cc=cc/SQRT2; ss=ss/SQRT2
    end if

    if(mod(k1+k2,2)==1) then
      ! write(6,'(2i5)') k1, k2
    else
      if(abs(cc) > eps .OR. abs(ss) > eps) then
        if(k1 /= k2) then
          write(6,'(2i5,2f21.15," *")') k1, k2, cc, ss
        else
          write(6,'(2i5,2f21.15)') k1, k2, cc, ss
        endif
      end if
    end if
  enddo
enddo

!do k1=0,kmax
!  write(6,*) '   k1   k2  cos(k1 x) cos(k2 x)  sin(k1 x) sin(k2 x)  sin(k1 x) cos(k2 x)'
!  do k2=0,kmax
!    cc=zero; ss=zero; cs=zero
!    do n=1,nphi
!      cc = cc+tphi_cos(k1,n)*tphi_cos(k2,n)
!      ss = ss+tphi_sin(k1,n)*tphi_sin(k2,n)
!      cs = cs+tphi_cos(k1,n)*tphi_sin(k2,n)
!    enddo
!    write(6,'(2i5,3f21.16)') k1, k2, 4*cc, 4*ss, 4*cs
!  enddo
!enddo

!write(6,*) '   k1  cos(k1 x)            sin(k1 x)'
!do k1=0,2*kmax
!    cc=zero; ss=zero
!    do n=1,nphi
!      cc = cc+tphi_cos(k1,n)
!      ss = ss+tphi_sin(k1,n)
!    enddo
!    write(6,'(i5,3f21.16)') k1, 4*cc, 4*ss
!enddo


END PROGRAM
