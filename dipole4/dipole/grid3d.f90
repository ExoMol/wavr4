! purpose:	preparation for Gaussian quadrature
! input:	allocated arrays and array boundaries
! returns:	grid points and basis functions computed at the grid points
! last modified: 03/04/02
!
! angular grids are modified so that only half of the usual grid is returned
! because we are using symmetry reduced grid size

MODULE grid3d

INTERFACE radial_grids
SUBROUTINE radial_grids(q1,q2,q3,t1,t2,t3,oner1,oner2,oner3)
  USE types
  USE param
!  USE potential
! IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: q1(nn(1)),t1(nn(1),nn(1)),  &
                           q2(nn(2)),t2(nn(2),nn(2)),  &
                           q3(nn(3)),t3(nn(3),nn(3)),  &
                                  oner1(nn(1),nn(1)),  &
                                  oner2(nn(2),nn(2)),  &
                                  oner3(nn(3),nn(3))
  END SUBROUTINE
END INTERFACE

INTERFACE angular_grids
SUBROUTINE angular_grids(theta1,theta,phi,ttheta1,ttheta,tphi)
  USE types
  USE param
!  USE potential
! IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: theta1(nt1),ttheta1(0:jmax,0:mmax,nt1)
  REAL(DP), INTENT(OUT) :: theta(nt),ttheta(0:lmax,0:kmax,nt)
  REAL(DP), INTENT(OUT) :: phi(nphi),tphi(0:2*kmax,nphi)
  END SUBROUTINE
END INTERFACE

INTERFACE grid_theta
SUBROUTINE grid_theta(theta,ttheta,lmax,kmax,nt)
  USE types
  USE base_lib, ONLY: plgndr, gammln, gauleg
  INTEGER(I4B), INTENT(IN) :: lmax,kmax,nt
  REAL(DP), INTENT(OUT) :: theta(nt),ttheta(0:lmax,0:kmax,nt)
  END SUBROUTINE
END INTERFACE

INTERFACE grid_phi
SUBROUTINE grid_phi(phi,tphi,kmax,nphi)
  USE types
  INTEGER(I4B), INTENT(IN) :: kmax,nphi
  REAL(DP), INTENT(OUT) :: phi(nphi),tphi(0:kmax,nphi)
  END SUBROUTINE
END INTERFACE

INTERFACE grid_phi_sin
SUBROUTINE grid_phi_sin(phi,tphi,kmax,nphi)
  USE types
  INTEGER(I4B), INTENT(IN) :: kmax,nphi
  REAL(DP), INTENT(OUT) :: phi(nphi),tphi(0:kmax,nphi)
  END SUBROUTINE
END INTERFACE

INTERFACE grid_phi2
SUBROUTINE grid_phi2(phi,tphi_cos,tphi_sin,kmax,nphi)
  USE types
  INTEGER(I4B), INTENT(IN) :: kmax,nphi
  REAL(DP), INTENT(OUT) :: phi(nphi),tphi_cos(0:kmax,nphi),tphi_sin(0:kmax,nphi)
  END SUBROUTINE
END INTERFACE


END MODULE

!  common grid calls

SUBROUTINE radial_grids(q1,q2,q3,t1,t2,t3,oner1,oner2,oner3)
  USE types
  USE param
!  USE potential
IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: q1(nn(1)),t1(nn(1),nn(1)),  &
                           q2(nn(2)),t2(nn(2),nn(2)),  &
                           q3(nn(3)),t3(nn(3),nn(3)),  &
                                  oner1(nn(1),nn(1)),  &
                                  oner2(nn(2),nn(2)),  &
                                  oner3(nn(3),nn(3))

nsmax = nn(1)*nn(2)*nn(3)

!  prepare stretching grids and Tkin in DVR basis
!  igq  controls the type of Gaussian Quadrature
!       presently 1 -- Morse-Laguerre
!                else  Spherical Oscillator

  if(nn(1) /= 1) then
    if(igq(1)==1) then 
      call stretch(1,q1,t1,oner1)
    else
      call stretch2(1,q1,t1,oner1)
    end if
  else
    q1(1)=qe(1)
  end if

  if(nn(2) /= 1) then
    if(igq(2)==1) then 
      call stretch(2,q2,t2,oner2)
    else
      call stretch2(2,q2,t2,oner2)
    end if
  else
    q2(1)=qe(2)
  end if

  if(nn(3) /= 1) then
    if(igq(3)==1) then 
      call stretch(3,q3,t3,oner3)
    else
      call stretch2(3,q3,t3,oner3)
    end if
  else
    q3(1)=qe(3)
  end if

END SUBROUTINE radial_grids

SUBROUTINE angular_grids(theta1,theta,phi,ttheta1,ttheta,tphi)
  USE types
  USE param
!  USE potential
IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: theta1(nt1),ttheta1(0:jmax,0:mmax,nt1)
  REAL(DP), INTENT(OUT) :: theta(nt),ttheta(0:lmax,0:kmax,nt)
  REAL(DP), INTENT(OUT) :: phi(nphi),tphi(0:2*kmax,nphi)
INTEGER(I4B) :: i

!  prepare angular grids

call grid_theta(theta1,ttheta1,jmax,mmax,nt1)
call grid_theta(theta,ttheta,lmax,kmax,nt)
call grid_phi(phi,tphi,kmax,nphi)

END SUBROUTINE angular_grids

! purpose:	Gauss-Legendre Quadrature
! input/output:	theta(i) - angular grid array in 0..Pi, cos(theta_i) = x_i
!		ttheta(l,k,i) - normalized associated Lengendre polynomial
!		P^k_l(x_i)*sqrt(weight_i)
!
! NOTE: the number of grid points nsz is determined by the theta array size.
! For the Gauss-Legendre quadrature to be successful nsz must be at least
! l_max+1 or greater. This is however for k's of the same parity.
! In our case k can be any => nsz must be 2*(l_max+1) or greater.
!
SUBROUTINE grid_theta(theta,ttheta,lmax,kmax,nt)
USE types
USE base_lib, ONLY: plgndr, gammln, gauleg
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: lmax,kmax,nt
REAL(DP), ALLOCATABLE :: x(:),wx(:)
REAL(DP), INTENT(OUT) :: theta(nt),ttheta(0:lmax,0:kmax,nt)
REAL(DP) :: norm
INTEGER(I4B) :: i,l,k, lkm, nsz, kmx
INTEGER(I4B) :: test_flag = 0

if(test_flag > 0) write(u6, '(" computing theta grid...")')

ttheta=0.0_dp
nsz = 2*size(theta)
kmx = ubound(ttheta,2)

allocate(x(nsz),wx(nsz))

call gauleg(-1.0_dp, 1.0_dp, x,wx)

!  reduced grid =  half of the full grid
theta = acos(x(nsz/2+1:nsz))
wx=sqrt(wx)
do l=0, ubound(ttheta,1)
  lkm = min( kmx,l)
  do k=0,lkm
!    norm = (-1)**k*sqrt((2._dp*l+1._dp)/2._dp*exp(gammln(l-k+1._dp)-gammln(l+k+1._dp)))
    norm = sqrt((2._dp*l+1._dp)/2._dp*exp(gammln(l-k+1._dp)-gammln(l+k+1._dp)))
    do i=1,nsz/2
      ttheta(l,k,i)= norm*plgndr(l,k,x(i+nsz/2))*wx(i+nsz/2)
    end do
  end do
end do

deallocate(x,wx)

!debug! tests of grid orthogonality
!
! this gives 1/2 due to halved integration region
!
!do l=0, ubound(ttheta,1)
!  lkm = min( kmx,l)
!  do k=0,lkm
!    norm = zero
!    do i=1,nsz/2
!      norm = norm+ ttheta(l,k,i)**2
!    end do
!    write(6,*) l,k,norm
!  end do
!end do

if(test_flag > 0) write(u6, '(" done.")')

END SUBROUTINE

! purpose:	Gauss (Fourier) Quadrature
! input/output:	phi(i) = angular grid array in 0..Pi
!		tphi(k,i) = cos(k*phi_i)/(2*Pi)*weight_i
!
! NOTE: the true integration region is 0..2*Pi but because of the inversion
! symmetry we reduce it to 0..Pi and double the result when integrating
! cos(?*phi)*V(?,?,phi) where V(?,?,phi)=V(?,?,2*Pi-phi).
! If there were no the symmetry, 0..2*Pi interval is required thus doubling
! the grid size.
!
! NOTE: original basis functions are exp(I*k*phi)/sqrt(2*Pi)
! so 1/(2*Pi) is norm^2. 
!
! NOTE: the number of grid points nsz is determined by the phi array size.
! For the Gauss quadrature to be successful nsz must be at least 
! k_max+1 or greater.
!
SUBROUTINE grid_phi(phi,tphi,kmax,nphi)
USE types
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: kmax,nphi
REAL(DP), INTENT(OUT) :: phi(nphi),tphi(0:2*kmax,nphi)
REAL(DP), ALLOCATABLE :: wphi(:)
REAL(DP) :: norm
INTEGER(I4B) :: i,k, nsz
INTEGER(I4B) :: test_flag = 0

if(test_flag > 0) write(u6, '(" computing phi grid...")')

tphi=0.0_dp
nsz = 2*size(phi)

allocate(wphi(nsz))

wphi=PI/real(nsz,dp)
do i=1,nsz/2
  phi(i)=PI*(i-0.5_dp)/real(nsz,dp)
  ! write(u6, '(2x,f20.16,2x,f20.16)') phi(i), wphi(i)
end do

!norm = 1._dp/2._dp/PI     !   this is the norm**2 for exp(i*k*phi)
norm = 1._dp/real(nsz,dp)  !  this is the same norm as above 
                           !  but INCLUDING weights which are uniform
do k= 0, ubound(tphi,1)
  tphi(k,:)= norm*cos(k*phi(:))
end do

deallocate(wphi)

if(test_flag > 0) write(u6, '(" done.")')

END SUBROUTINE


! this subroutine returns the same grid as grid_phi
! but values of sin(k x) instead of cos(k x)
! the arguments left the same as in grid_phi for compatibility
! it will over-write phi grid returned by grid_phi
! but that should not be a problem since the grid is the same

SUBROUTINE grid_phi_sin(phi,tphi,kmax,nphi)
USE types
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: kmax,nphi
REAL(DP), INTENT(OUT) :: phi(nphi),tphi(0:2*kmax,nphi)
REAL(DP), ALLOCATABLE :: wphi(:)
REAL(DP) :: norm
INTEGER(I4B) :: i,k, nsz
INTEGER(I4B) :: test_flag = 0

if(test_flag > 0) write(u6, '(" computing phi grid...")')

tphi=0.0_dp
nsz = 2*size(phi)

allocate(wphi(nsz))

wphi=PI/real(nsz,dp)
do i=1,nsz/2
  phi(i)=PI*(i-0.5_dp)/real(nsz,dp)
  ! write(u6, '(2x,f20.16,2x,f20.16)') phi(i), wphi(i)
end do

!norm = 1._dp/2._dp/PI     !   this is the norm**2 for exp(i*k*phi)
norm = 1._dp/real(nsz,dp)  !  this is the same norm as above 
                           !  but INCLUDING weights which are uniform
do k= 0, ubound(tphi,1)
  tphi(k,:)= norm*sin(k*phi(:))
end do

deallocate(wphi)

if(test_flag > 0) write(u6, '(" done.")')

END SUBROUTINE


SUBROUTINE grid_phi2(phi,tphi_cos,tphi_sin,kmax,nphi)
USE types
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: kmax,nphi
REAL(DP), INTENT(OUT) :: phi(nphi),tphi_cos(0:2*kmax,nphi),tphi_sin(0:2*kmax,nphi)
REAL(DP), ALLOCATABLE :: wphi(:)
REAL(DP) :: norm
INTEGER(I4B) :: i,k, nsz
INTEGER(I4B) :: test_flag = 0

if(test_flag > 0) write(u6, '(" computing phi grid...")')

tphi_cos=0.0_dp
tphi_sin=0.0_dp
nsz = 2*size(phi)

allocate(wphi(nsz))

wphi=PI/real(nsz,dp)
do i=1,nsz/2
  phi(i)=PI*(i-0.5_dp)/real(nsz,dp)
  ! write(u6, '(2x,f20.16,2x,f20.16)') phi(i), wphi(i)
end do

!norm = 1._dp/2._dp/PI      !  this is the norm**2 for exp(i*k*phi)
norm = 1._dp/real(nsz,dp)  !  this is the same norm as above 
                           !  but INCLUDING weights which are uniform

norm = sqrt(norm)          !  we need to take square root of norm because 
                           !  unlike in WAVR4 we use product of two functions

do k= 0, 2*kmax
  tphi_cos(k,:)= norm*cos(k*phi(:))
  tphi_sin(k,:)= norm*sin(k*phi(:))
end do

deallocate(wphi)

!debug! tests of grid orthogonality
!
! integrals over cos and sin give 1/4 due to halved integration region
! and due to the fact that the true norms are different for k=0 and k>0
! integral over cos(k phi) k=0 gives 1/2 which is correct owing to the
! fact that the integration region is halved
!
!do k= 0, kmax
!  norm = zero
!  do i=1,nsz/2
!    norm = norm + tphi_cos(k,i)**2
!  end do
!  write(6,*) 'cos',k, norm
!end do
!
!do k= 0, kmax
!  norm = zero
!  do i=1,nsz/2
!    norm = norm + tphi_sin(k,i)**2
!  end do
!  write(6,*) 'sin',k, norm
!end do
!
!do k= 0, kmax
!  norm = zero
!  do i=1,nsz/2
!    norm = norm + tphi_cos(k,i)*tphi_sin(k,i)
!  end do
!  write(6,*) 'c_s',k, norm
!end do


if(test_flag > 0) write(u6, '(" done.")')

END SUBROUTINE
