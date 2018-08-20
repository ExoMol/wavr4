MODULE workarrays
  USE types
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: threej0(:,:,:)
  REAL(DP), ALLOCATABLE :: tp1(:),tp2(:),tp(:),v3(:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: t1(:,:),t2(:,:),t3(:,:)
  REAL(DP), ALLOCATABLE :: oner1(:,:),oner2(:,:),oner3(:,:)
  TYPE(triple), ALLOCATABLE :: krn(:)
!------------------------------------------------------------------
! removed from main. to be share ONLY with stage# and mainj subroutines
  INTEGER(I4B), ALLOCATABLE :: indx(:,:)
  TYPE(triple), ALLOCATABLE :: idx1(:,:,:),idx2(:,:),idx3(:)
!  TYPE(eigenvector), ALLOCATABLE :: eigvec1(:), eigvec2(:) ! , eigvec3(:)
!  INTEGER(I4B), ALLOCATABLE :: rec_num1(:,:,:,:),rec_num2(:,:,:)

  TYPE(matrixarray), ALLOCATABLE :: vectrs1(:,:,:), vectrs2(:,:)
  TYPE(vectorarray), ALLOCATABLE :: energs1(:,:,:), energs2(:,:)
!------------------------------------------------------------------
!  ONLY if the extra (HHSYM) permutation symmetry is present:
!  integer array for keeping the symmetry of the angular basis and
!  logical arrays for keeping the symmetry of stage1 & stage2
!  eigenfunctions with i2=i1:  symmetric .TRUE.;  asymmetric .FALSE.
!
  LOGICAL(LGT), ALLOCATABLE :: s1sym(:,:,:,:),s2sym(:,:,:)
  TYPE(triple), ALLOCATABLE :: idx2sym(:,:),idx2asy(:,:)
  INTEGER(I4B), ALLOCATABLE :: angsym(:,:)
!------------------------------------------------------------------
  !$OMP THREADPRIVATE(threej0,tp1,tp2,tp,v3)
  SAVE threej0,tp1,tp2,tp,v3

!$  REAL(DP) :: start_time, end_time

END MODULE
