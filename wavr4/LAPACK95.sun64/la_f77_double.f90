      MODULE F77_LAPACK
!
!  cut down version of
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!

      INTERFACE

      FUNCTION DLAMCH( CMACH )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: DLAMCH
         CHARACTER(LEN=1), INTENT(IN) :: CMACH
      END FUNCTION DLAMCH

      END INTERFACE

      INTERFACE

      FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
         INTEGER :: ILAENV
         CHARACTER(LEN=*), INTENT(IN) :: NAME, OPTS
         INTEGER, INTENT(IN) :: ISPEC, N1, N2, N3, N4
      END FUNCTION ILAENV

      END INTERFACE


      INTERFACE

       SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,     &
     &                    ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL,     &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO, RANGE
         INTEGER, INTENT(IN) :: LDZ, N, IL, IU
         INTEGER, INTENT(OUT) :: INFO, IWORK(*), M, IFAIL(*)
         REAL(WP), INTENT(IN) :: VL, VU, ABSTOL
         REAL(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: W(*), Z(LDZ,*), WORK(*)
      END SUBROUTINE DSPEVX

       END INTERFACE


      INTERFACE

       SUBROUTINE DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK,    &
     &                    IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDZ, N, LWORK, LIWORK
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: W(*), Z(LDZ,*), WORK(*)
      END SUBROUTINE DSPEVD

       END INTERFACE


      INTERFACE

       SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDZ, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: W(*), Z(LDZ,*), WORK(*)
      END SUBROUTINE DSPEV

       END INTERFACE



      INTERFACE

       SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         ! USE SUNPERF
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSYEV

       END INTERFACE


      INTERFACE

       SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
     &                    LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDA, LIWORK, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSYEVD

       END INTERFACE


       INTERFACE

       SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
     &                    ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,    &
     &                    IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: N, IL, IU, LDZ, LDA, LWORK, LIWORK
         INTEGER, INTENT(OUT) :: M
         INTEGER, INTENT(OUT) :: ISUPPZ(*)
         REAL(WP), INTENT(IN) ::  ABSTOL, VL, VU
         INTEGER, INTENT(OUT) ::  IWORK(*)
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: W(*)
       END SUBROUTINE  DSYEVR

      END INTERFACE


      INTERFACE

       SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
     &                    ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,     &
     &                    IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: IL, IU, LDA, LDZ, LWORK, N
         INTEGER, INTENT(OUT) :: INFO, M
         INTEGER, INTENT(OUT) :: IFAIL(*), IWORK(*)
         REAL(WP), INTENT(IN) :: ABSTOL, VL, VU
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
      END SUBROUTINE DSYEVX

       END INTERFACE


      END MODULE F77_LAPACK
