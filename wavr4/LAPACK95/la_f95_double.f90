      MODULE F95_LAPACK
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!

      INTERFACE LA_SPEVX

       SUBROUTINE LA_SPEVX( AP, W, UPLO, Z, VL, VU, IL, IU, M, IFAIL, &
     &                        ABSTOL, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
           INTEGER, INTENT(IN), OPTIONAL :: IL, IU
           INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
           REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
           INTEGER, INTENT(OUT), OPTIONAL :: IFAIL(:)
           REAL(WP), INTENT(OUT), OPTIONAL :: Z(:,:)
           REAL(WP), INTENT(INOUT) :: AP(:)
           REAL(WP), INTENT(OUT) :: W(:)
        END SUBROUTINE

      END INTERFACE


      INTERFACE LA_SPEVD

       SUBROUTINE LA_SPEVD( AP, W, UPLO, Z, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
           INTEGER, INTENT(OUT), OPTIONAL :: INFO
           REAL(WP), INTENT(INOUT) :: AP(:)
           REAL(WP), INTENT(OUT) :: W(:)
           REAL(WP), INTENT(OUT), OPTIONAL :: Z(:,:)
        END SUBROUTINE

      END INTERFACE


      INTERFACE LA_SPEV

       SUBROUTINE LA_SPEV( AP, W, UPLO, Z, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: UPLO
           INTEGER, INTENT(OUT), OPTIONAL :: INFO
           REAL(WP), INTENT(INOUT) :: AP(:)
           REAL(WP), INTENT(OUT) :: W(:)
           REAL(WP), INTENT(OUT), OPTIONAL :: Z(:,:)
        END SUBROUTINE

      END INTERFACE


      INTERFACE LA_SYEV

       SUBROUTINE LA_SYEV( A, W, JOBZ, UPLO, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
           INTEGER, INTENT(OUT), OPTIONAL :: INFO
           REAL(WP), INTENT(INOUT) :: A(:,:)
           REAL(WP), INTENT(OUT) :: W(:)
        END SUBROUTINE

      END INTERFACE


      INTERFACE LA_SYEVD

       SUBROUTINE LA_SYEVD( A, W, JOBZ, UPLO, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
           INTEGER, INTENT(OUT), OPTIONAL :: INFO
           REAL(WP), INTENT(INOUT) :: A(:,:)
           REAL(WP), INTENT(OUT) :: W(:)
        END SUBROUTINE

      END INTERFACE


        INTERFACE LA_SYEVR

       SUBROUTINE LA_SYEVR( A, W, JOBZ, UPLO, VL, VU, IL, IU, M,      &
     &                        ISUPPZ, ABSTOL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP), INTENT(INOUT) :: A(:,:)
         REAL(WP), INTENT(OUT) :: W(:)
         CHARACTER(LEN=1), INTENT(IN), OPTIONAL ::  JOBZ, UPLO
         INTEGER, INTENT(OUT), OPTIONAL :: INFO
         REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
         INTEGER, INTENT(IN), OPTIONAL :: IL, IU
         INTEGER, INTENT(OUT), OPTIONAL :: M
         INTEGER, INTENT(OUT), OPTIONAL, TARGET :: ISUPPZ(:)
       END SUBROUTINE

       END INTERFACE


      INTERFACE LA_SYEVX

       SUBROUTINE LA_SYEVX( A, W, JOBZ, UPLO, VL, VU, IL, IU, M,      &
     &                        IFAIL, ABSTOL, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
           INTEGER, INTENT(IN), OPTIONAL :: IL, IU
           INTEGER, INTENT(OUT), OPTIONAL :: INFO, M
           REAL(WP), INTENT(IN), OPTIONAL :: ABSTOL, VL, VU
           INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IFAIL(:)
           REAL(WP), INTENT(INOUT) :: A(:,:)
           REAL(WP), INTENT(OUT) :: W(:)
        END SUBROUTINE

      END INTERFACE



      END MODULE F95_LAPACK
