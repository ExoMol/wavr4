! dipole functions
!
!  USE param

MODULE dipole

INTERFACE mux
  FUNCTION mux(q1,q2,q3,theta1,theta,phi)
  USE types
  REAL(DP), INTENT(IN) :: q1,q2,q3,theta1,theta,phi
  REAL(DP) :: mux
  END FUNCTION mux
END INTERFACE

INTERFACE muy
  FUNCTION muy(q1,q2,q3,theta1,theta,phi)
  USE types
  REAL(DP), INTENT(IN) :: q1,q2,q3,theta1,theta,phi
  REAL(DP) :: muy
  END FUNCTION muy
END INTERFACE

INTERFACE muz
  FUNCTION muz(q1,q2,q3,theta1,theta,phi)
  USE types
  REAL(DP), INTENT(IN) :: q1,q2,q3,theta1,theta,phi
  REAL(DP) :: muz
  END FUNCTION muz
END INTERFACE

END MODULE

FUNCTION mux(q1,q2,q3,theta1,theta,phi)
  USE types
IMPLICIT NONE
  REAL(DP) :: mux
  REAL(DP), INTENT(IN) ::  q1,q2,q3,theta1,theta,phi
  !mux = one
  mux = q1*sin(theta1) + q2*sin(theta)*cos(phi)

END FUNCTION

FUNCTION muy(q1,q2,q3,theta1,theta,phi)
  USE types
IMPLICIT NONE
  REAL(DP) :: muy
  REAL(DP), INTENT(IN) ::  q1,q2,q3,theta1,theta,phi
  !muy = one
  muy = q2*sin(theta)*sin(phi)

END FUNCTION

FUNCTION muz(q1,q2,q3,theta1,theta,phi)
  USE types
IMPLICIT NONE
  REAL(DP) :: muz
  REAL(DP), INTENT(IN) ::  q1,q2,q3,theta1,theta,phi
  !muz = one
  muz = q1*cos(theta1) + q2*cos(theta)

END FUNCTION


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! muz = mu1* cos(theta1) + mu2* cos(theta)
! mux = mu1* sin(theta1) + mu2* sin(theta)*cos(phi)
! muy =                         sin(theta)*sin(phi)
