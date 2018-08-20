MODULE workarrays
  USE types
  IMPLICIT NONE
  REAL(DP), ALLOCATABLE :: threej0(:,:,:)
  REAL(DP), ALLOCATABLE :: tp1(:),tp2(:),tp(:),v3(:,:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: t1(:,:),t2(:,:),t3(:,:)
  REAL(DP), ALLOCATABLE :: oner1(:,:),oner2(:,:),oner3(:,:)
  TYPE(triple), ALLOCATABLE :: krn(:)
! removed from main. to be share ONLY with mainj
  INTEGER(I4B), ALLOCATABLE :: indx(:,:)
  TYPE(triple), ALLOCATABLE :: idx1(:,:,:),idx2(:),idx3(:)
END MODULE
