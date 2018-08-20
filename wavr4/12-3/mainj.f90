SUBROUTINE mainj()
  USE types
  USE param, ONLY: jr, p, jp, jl_parity,l_parity,j_parity,enzero
  USE param, ONLY: nn, test_flag,oner_flag
  USE param, ONLY: icut4,encut4 ! ,imargin4,margin4
  USE workarrays, ONLY: oner3, krn
  USE workarrays, ONLY: indx, idx1,idx2,idx3
!
!  actually only indx and idx3 are really used; idx1 and idx2 are read afresh
!
  USE angles, ONLY: angleT3D
  USE base_lib, ONLY: time_stamp, diag_p, xn
  USE LA_PRECISION, ONLY: WP => DP
!  these are regular diagonalisers
!  USE F95_LAPACK, ONLY: LA_SYEV, LA_SYEVD, LA_SYEVX, LA_SYEVR
!  these are packed diagonalisers
  USE F95_LAPACK, ONLY: LA_SPEV, LA_SPEVD, LA_SPEVX
IMPLICIT NONE
!
!  common arrays
!
REAL(DP), ALLOCATABLE :: h3(:,:),e3(:),h5(:,:),e5(:),h6(:),e6(:),h7(:),e7(:)
REAL(DP), ALLOCATABLE :: vectors(:,:)
!
!  Primary stage specific arrays and variables: stages 1, 2, 3
!
TYPE(eigenvector), ALLOCATABLE :: eigvec1(:), eigvec2(:), eigvec3(:)
INTEGER(I4B), ALLOCATABLE :: rec_num1(:,:,:,:),rec_num2(:,:)
INTEGER(I4B) :: istage1,istage2,istage3
INTEGER(I4B) :: ist1_max,ist2_max ! ,ist3_max
!
!  Alternative stage specific arrays and variables: stages 1, 2, 3
!
TYPE(triple), ALLOCATABLE :: idx1A(:,:,:),idx2A(:)
TYPE(eigenvector), ALLOCATABLE :: eigvec1A(:), eigvec2A(:), eigvec3A(:)
INTEGER(I4B), ALLOCATABLE :: rec_num1A(:,:,:,:),rec_num2A(:,:)
INTEGER(I4B) :: istage1A,istage2A,istage3A
INTEGER(I4B) :: ist1_maxA,ist2_maxA ! ,ist3_maxA
!
!  vectors
!
REAL(DP), POINTER :: vec11(:),vec12(:), vec21(:),vec22(:), vec31(:),vec32(:)
!
!  aux arrays
!
TYPE(matrixarray), ALLOCATABLE :: vm(:,:,:)
!
!  auxiliary variables   NOTE:  kr is local variable
!
REAL(DP) :: temp, factor_jrkr  ! , memsize
INTEGER(I4B) :: n,n1, i1,i2,i3, nlev,nout,idiag, nf
INTEGER(I4B) :: ia,ia1,ib,ic,id,id1,im,im1
INTEGER(I4B) :: i, kr, krtop, icounter, nsize
! INTEGER(I4B) :: l,j, ib1, ic1

write(u6,'(" entering mainj: J=",i3)') jr

if(jr==0) then
  ! tell what to do if J=0
  ! presumably just print the levels
  write(6,*) ' mainj: J=', jr
  return
end if

! if J > 0  proceed further

allocate(idx1A(nn(1),nn(2),nn(3)),idx2A(nn(3)))

!  temporary external storage for off-diagonal in K matrix elements
!  it should be open by now
! open(unit=u16,file='x_h6.dat', form='UNFORMATTED')

!  for Ar2-HF the final matrix should not be very big
!  because vibration-rotation separation should be OK

nlev = idx3(jr)%end
allocate(h7(nlev*(nlev+1)/2))

h7 = zero

write(u6,'("  h7 size=",i6,/)') nlev

!  first we deal with the special case

kr = jr
write(u6,'("  ###  K=",i3,/)') kr

! ===================================================================
! stage 3 vectors are read first to place them at "the bottom" of free
! memory. So when all stage 1,1A,2,2A eigenvectors are deallocated
! this does not create any "gaps" and leave the memory contiguous.
! This should help avoid memory problems when reading h5 matrix

write(u6,'(" reading stage 3...")')
rewind(u31(kr))
read(u31(kr)) istage3,nsize
allocate(eigvec3(istage3))
if(test_flag >2) write(u6,'(" print out levels:")')

do icounter = 1,istage3
  allocate(eigvec3(icounter)%vec(nsize),eigvec3(icounter)%eng)
  read(u31(kr)) eigvec3(icounter)%eng,eigvec3(icounter)%vec(:)
end do
icounter= icounter-1

if(istage3== icounter) then
  write(u6,'(" done:  read in",i8," eigenvectors",/)') icounter
else
  write(u6,'(" error in reading stage 3:")')
  write(u6,'(" istage3=",i8," icounter=",i8)') istage3,icounter
  stop
end if

! ===================================================================

write(u6,'(" reading stage 1...")')
rewind(u11(kr))
read(u11(kr)) istage1, ist1_max
read(u11(kr)) idx1
allocate(eigvec1(istage1),rec_num1(ist1_max,nn(1),nn(2),nn(3)))
rec_num1 = 0
nsize  = krn(kr)%size
icounter = 0
if(test_flag >2) write(u6,'(" print out levels:")')
do i3=1,nn(3); do i2=1,nn(2); do i1=1,nn(1)
  do ia=1,idx1(i1,i2,i3)%size
    icounter = icounter+1
    rec_num1(ia,i1,i2,i3)=icounter
    allocate(eigvec1(icounter)%vec(nsize), eigvec1(icounter)%eng)
    read(u11(kr)) eigvec1(icounter)%eng,eigvec1(icounter)%vec(:)
    if(test_flag >2) write(u8,'(i6,4i4,f12.3)')  &
       icounter, i1,i2,i3, ia, eigvec1(icounter)%eng
  end do
end do; end do; end do

if(istage1== icounter) then
  write(u6,'(" done:  read in",i8," eigenvectors",/)') icounter
else
  write(u6,'(" error in reading stage 1:")')
  write(u6,'(" istage1=",i8," icounter=",i8)') istage1,icounter
  stop
end if

! ===================================================================

write(u6,'(" reading stage 2...")')
rewind(u21(kr))
read(u21(kr)) istage2,ist2_max
read(u21(kr)) idx2
allocate(eigvec2(istage2),rec_num2(ist2_max,nn(3)))
rec_num2 = 0
icounter = 0
if(test_flag >2) write(u6,'(" print out levels:")')
do id =1,nn(3)
  nsize = idx1(nn(1),nn(2),id)%end
  if(nsize < 1) cycle
  do n=1,idx2(id)%size
    icounter = icounter+1
    rec_num2(n,id) = icounter
    allocate(eigvec2(icounter)%vec(nsize),eigvec2(icounter)%eng)
    read(u21(kr)) eigvec2(icounter)%eng,eigvec2(icounter)%vec(:)
    if(test_flag >2)  &
      write(u8,'(i6,2i4,f12.3)') icounter, id, n, eigvec2(icounter)%eng
  end do
end do

if(istage2== icounter) then
  write(u6,'(" done:  read in",i8," eigenvectors",/)') icounter
else
  write(u6,'(" error in reading stage 2:")')
  write(u6,'(" istage2=",i8," icounter=",i8)') istage2,icounter
  stop
end if

! ===================================================================

!  we start formation of h7 from highest K

!  in principle
! krtop= min(krmax,jr)
!  but we assume krmax is unlimited so
krtop = jr

!  the sub-block kr1 = kr (=jr) is already diagonal

do i=1,idx3(jr)%size
  n=idx3(jr)%start+i-1
  ! write(6,*) ' test energy:',n, eigvec3(i)%eng
  h7(xn(n,n)) = eigvec3(i)%eng
end do

kr_loop: do kr=krtop-1,0,-1

  write(u6,'("  ###  K=",i3,/)') kr

  !  now we form upper off-diagonal block <K|  |K+1>
  !  the left-hand vectors refer to K
  !    indices are marked by 1 (for prime)
  !    vector arrays     by A
  !  the right-hand vectors refer to K+1
  !    indiced are NOT marked
  !    vector arrays are NOT marked
  ! ===================================================================

  write(u6,'(" reading stage 1A...")')
  rewind(u11(kr))
  read(u11(kr)) istage1A, ist1_maxA
  read(u11(kr)) idx1A
  allocate(eigvec1A(istage1A),rec_num1A(ist1_maxA,nn(1),nn(2),nn(3)))
  rec_num1A = 0
  nsize  = krn(kr)%size
  icounter = 0
  if(test_flag >2) write(u6,'(" print out levels:")')
  do i3=1,nn(3); do i2=1,nn(2); do i1=1,nn(1)
    do ia=1,idx1A(i1,i2,i3)%size
      icounter = icounter+1
      rec_num1A(ia,i1,i2,i3)=icounter
      allocate(eigvec1A(icounter)%vec(nsize), eigvec1A(icounter)%eng)
      read(u11(kr)) eigvec1A(icounter)%eng,eigvec1A(icounter)%vec(:)
      if(test_flag >2) write(u8,'(i6,4i4,f12.3)')  &
         icounter, i1,i2,i3, ia, eigvec1A(icounter)%eng
    end do
  end do; end do; end do

  if(istage1A== icounter) then
    write(u6,'(" done:  read in",i8," eigenvectors",/)') icounter
  else
    write(u6,'(" error in reading stage 1:")')
    write(u6,'(" istage1A=",i8," icounter=",i8)') istage1A,icounter
    stop
  end if

  ! ===================================================================

  write(u6,'(" reading stage 2A...")')
  rewind(u21(kr))
  read(u21(kr)) istage2A,ist2_maxA
  read(u21(kr)) idx2A
  allocate(eigvec2A(istage2A),rec_num2A(ist2_maxA,nn(3)))
  rec_num2A = 0
  icounter = 0
  if(test_flag >2) write(u6,'(" print out levels:")')
  do id =1,nn(3)
    nsize = idx1A(nn(1),nn(2),id)%end
    if(nsize < 1) cycle
    do n=1,idx2A(id)%size
      icounter = icounter+1
      rec_num2A(n,id) = icounter
      allocate(eigvec2A(icounter)%vec(nsize),eigvec2A(icounter)%eng)
      read(u21(kr)) eigvec2A(icounter)%eng,eigvec2A(icounter)%vec(:)
      if(test_flag >2)  &
        write(u8,'(i6,2i4,f12.3)') icounter, id, n, eigvec2A(icounter)%eng
    end do
  end do

  if(istage2A== icounter) then
    write(u6,'(" done:  read in",i8," eigenvectors",/)') icounter
  else
    write(u6,'(" error in reading stage 2:")')
    write(u6,'(" istage2A=",i8," icounter=",i8)') istage2A,icounter
    stop
  end if

  ! ===================================================================

  rewind(u16)
  ! allocate(h5(istage2,istage2A))
  ! to save memory use just a column and make h5 on disk
  allocate(h6(istage2A))

  !  temporary need h3 & e3 again
  allocate(h3(krn(kr)%size,krn(kr+1)%size),e3(krn(kr)%size))
  
  h3(:,:) = zero
  call angleT3D(h3,indx, &
    krn(kr)%start,krn(kr)%size,krn(kr+1)%start,krn(kr+1)%size)

  write(6,*) ' istage2 =', istage2 ,' idx2 =', idx2(nn(3))%end
  write(6,*) ' istage2A=', istage2A,' idx2A=', idx2A(nn(3))%end

  do id=1,nn(3)
    ! write(6,*) ' id=', id
    ! if idx2(id)%size < 1 there's nothing to do
    if(idx2(id)%size < 1) cycle

    !  for efficient computation we need to precompute off-diagonal
    !  in "id" vector-vector products

    !  if 1/R^2 is treated in DVR (oner_flag==0)
    !  then oner3(id1,id) = 0 if id1 /= id

    ! allocate auxiliary array of matrices to hold vector products
    allocate(vm(nn(1),nn(2),nn(3)))

      do id1=1,nn(3)

        if(oner_flag==0 .AND. id1 /= id) cycle

        if(idx2A(id1)%size < 1) cycle
        do ic=1,nn(2)
          do ib=1,nn(1)
            if(idx1(ib,ic,id)%size < 1 .OR. idx1A(ib,ic,id1)%size < 1) cycle

            allocate(vm(ib,ic,id1)%mat(idx1A(ib,ic,id1)%size,idx1(ib,ic,id)%size))
            do ia=1,idx1(ib,ic,id)%size
              vec12 => eigvec1(rec_num1(ia,ib,ic,id))%vec(:)

              e3 = matmul(h3,vec12)

              do ia1= 1,idx1A(ib,ic,id1)%size
                vec11 => eigvec1A(rec_num1A(ia1,ib,ic,id1))%vec(:)
                vm(ib,ic,id1)%mat(ia1,ia)= dot_product(vec11,e3)
              end do
            end do

          end do
        end do

      end do

      ! write(6,*) '   finished preparing aux matrices'

    do im=1,idx2(id)%size
      n = idx2(id)%start+im-1

      h6(:) = zero

        vec22 => eigvec2(rec_num2(im,id))%vec(:)

      do id1=1,nn(3)

        if(oner_flag==0 .AND. id1 /= id) cycle

        do im1=1,idx2A(id1)%size
          n1 = idx2A(id1)%start+im1-1

            vec21 => eigvec2A(rec_num2A(im1,id1))%vec(:)

          temp = 0.0D0
          do ic=1,nn(2)
            do ib=1,nn(1)
              if(idx1(ib,ic,id)%size < 1 .OR. idx1A(ib,ic,id1)%size < 1) cycle

              allocate(e6(idx1(ib,ic,id)%size))

              e6 = matmul(  &
                vec21(idx1A(ib,ic,id1)%start:idx1A(ib,ic,id1)%end), &
                vm(ib,ic,id1)%mat(:,:))
              temp = temp+ dot_product(e6, &
                vec22(idx1(ib,ic,id)%start:idx1(ib,ic,id)%end))

              deallocate(e6)
            end do
          end do

          h6(n1) = temp* oner3(id1,id)

        end do
      end do

      !  store the composed column of h6 to disk
      ! write(6,*) ' writing record=',n
      write(u16) h6(:)

    end do

    !  deallocate auxiliary array of matrices

      do id1=1,nn(3)

        if(oner_flag==0 .AND. id1 /= id) cycle

        if(idx2A(id1)%size < 1) cycle
        do ic=1,nn(2)
          do ib=1,nn(1)
            if(idx1(ib,ic,id)%size < 1 .OR. idx1A(ib,ic,id1)%size < 1) cycle
            deallocate(vm(ib,ic,id1)%mat)
          end do
        end do

      end do
      deallocate(vm)

    ! write(6,*) '  id end'
  end do

  deallocate(h3,e3,h6)

  nullify(vec22,vec21,vec11,vec12)

  !  deallocate all previous stage 1 & 2 eigenvectors and arrays

  do n=1,istage1;  deallocate(eigvec1(n)%eng,eigvec1(n)%vec);  end do
  do n=1,istage2;  deallocate(eigvec2(n)%eng,eigvec2(n)%vec);  end do

  deallocate(eigvec1,eigvec2,rec_num1,rec_num2)

  do n=1,istage1A;  deallocate(eigvec1A(n)%eng,eigvec1A(n)%vec);  end do
  do n=1,istage2A;  deallocate(eigvec2A(n)%eng,eigvec2A(n)%vec);  end do

  deallocate(eigvec1A,eigvec2A,rec_num1A,rec_num2A)

  !  hopefully this gave us some space to read in h5
  !  for final vector-matrix-vector operations

  write(u6,'(" allocating and reading in h5...")')

  allocate(h5(istage2A,istage2))
  rewind(u16)
  do n=1,istage2
    ! write(6,*) ' reading record=',n
    read(u16) h5(:,n)
  end do

  write(u6,'(" done.")')

  ! ===================================================================

  write(u6,'(" reading stage 3A...")')
  rewind(u31(kr))
  read(u31(kr)) istage3A,nsize
  allocate(eigvec3A(istage3A))
  if(test_flag >1) write(u6,'(" print out levels:")')

  do icounter = 1,istage3A
    allocate(eigvec3A(icounter)%vec(nsize),eigvec3A(icounter)%eng)
    read(u31(kr)) eigvec3A(icounter)%eng,eigvec3A(icounter)%vec(:)
  end do
  icounter= icounter-1

  if(istage3A== icounter) then
    write(u6,'(" done:  read in",i8," eigenvectors",/)') icounter
  else
    write(u6,'(" error in reading stage 3:")')
    write(u6,'(" istage3A=",i8," icounter=",i8)') istage3A,icounter
    stop
  end if

  ! ===================================================================

  ! remember istage3A = idx3(kr)%size;  istage3 = idx3(kr+1)%size
  ! write(6,*) ' istage3A=', istage3A,' idx3 K  =', idx3(kr)%size
  ! write(6,*) ' istage3 =', istage3 ,' idx3 K+1=', idx3(kr+1)%size
  write(6,*) ' idx3 K    size =', idx3(kr)%size, &
             ' start=', idx3(kr)%start,' end=', idx3(kr)%end
  write(6,*) ' idx3 K+1  size =', idx3(kr+1)%size, &
             ' start=', idx3(kr+1)%start,' end=', idx3(kr+1)%end

  ! final reduction h5 -> h7 sub-block
  ! aux array

  allocate(e5(istage2))

  do i1=1,idx3(kr)%size                 ! istage3A
    vec31 => eigvec3A(i1)%vec(:)
    n1=idx3(kr)%start+i1-1

    e5 = matmul(vec31(:),h5(:,:))

    do i=1,idx3(kr+1)%size                  ! istage3
      vec32 => eigvec3(i)%vec(:)
      n= idx3(kr+1)%start+i-1

      !  to test the influence of <K| |K+1> sub-block comment it
      h7(xn(n1,n))= dot_product(e5(:),vec32(:))
    end do
  end do

  deallocate(e5,h5)
  nullify(vec32,vec31)

  ! ===================================================================
  !  deallocate 3 and 3A eigenvectors

  do n = 1,istage3;  deallocate(eigvec3(n)%vec,eigvec3(n)%eng);  end do
  do n = 1,istage3A;  deallocate(eigvec3A(n)%vec,eigvec3A(n)%eng);  end do

  deallocate(eigvec3,eigvec3A)

  ! istage3= istage3A

  ! ===================================================================
  ! now we need to read in all the eigenvectors again

  write(u6,'(" reading stage 3...")')
  rewind(u31(kr))
  read(u31(kr)) istage3,nsize
  allocate(eigvec3(istage3))
  if(test_flag >2) write(u6,'(" print out levels:")')

  do icounter = 1,istage3
    allocate(eigvec3(icounter)%vec(nsize),eigvec3(icounter)%eng)
    read(u31(kr)) eigvec3(icounter)%eng,eigvec3(icounter)%vec(:)
  end do
  icounter= icounter-1

  if(istage3== icounter) then
    write(u6,'(" done:  read in",i8," eigenvectors",/)') icounter
  else
    write(u6,'(" error in reading stage 3:")')
    write(u6,'(" istage3=",i8," icounter=",i8)') istage3,icounter
    stop
  end if

  ! ===================================================================

  write(u6,'(" reading stage 1...")')
  rewind(u11(kr))
  read(u11(kr)) istage1, ist1_max
  read(u11(kr)) idx1
  allocate(eigvec1(istage1),rec_num1(ist1_max,nn(1),nn(2),nn(3)))
  rec_num1 = 0
  nsize  = krn(kr)%size
  icounter = 0
  if(test_flag >2) write(u6,'(" print out levels:")')
  do i3=1,nn(3); do i2=1,nn(2); do i1=1,nn(1)
    do ia=1,idx1(i1,i2,i3)%size
      icounter = icounter+1
      rec_num1(ia,i1,i2,i3)=icounter
      allocate(eigvec1(icounter)%vec(nsize), eigvec1(icounter)%eng)
      read(u11(kr)) eigvec1(icounter)%eng,eigvec1(icounter)%vec(:)
      if(test_flag >2) write(u8,'(i6,4i4,f12.3)')  &
         icounter, i1,i2,i3, ia, eigvec1(icounter)%eng
    end do
  end do; end do; end do

  if(istage1== icounter) then
    write(u6,'(" done:  read in",i8," eigenvectors",/)') icounter
  else
    write(u6,'(" error in reading stage 1:")')
    write(u6,'(" istage1=",i8," icounter=",i8)') istage1,icounter
    stop
  end if

  ! ===================================================================

  write(u6,'(" reading stage 2...")')
  rewind(u21(kr))
  read(u21(kr)) istage2,ist2_max
  read(u21(kr)) idx2
  allocate(eigvec2(istage2),rec_num2(ist2_max,nn(3)))
  rec_num2 = 0
  icounter = 0
  if(test_flag >2) write(u6,'(" print out levels:")')
  do id =1,nn(3)
    nsize = idx1(nn(1),nn(2),id)%end
    if(nsize < 1) cycle
    do n=1,idx2(id)%size
      icounter = icounter+1
      rec_num2(n,id) = icounter
      allocate(eigvec2(icounter)%vec(nsize),eigvec2(icounter)%eng)
      read(u21(kr)) eigvec2(icounter)%eng,eigvec2(icounter)%vec(:)
      if(test_flag >2)  &
        write(u8,'(i6,2i4,f12.3)') icounter, id, n, eigvec2(icounter)%eng
    end do
  end do

  if(istage2== icounter) then
    write(u6,'(" done:  read in",i8," eigenvectors",/)') icounter
  else
    write(u6,'(" error in reading stage 2:")')
    write(u6,'(" istage2=",i8," icounter=",i8)') istage2,icounter
    stop
  end if

  ! ===================================================================
  ! finally do K1=K diagonal sub-block
  ! (should be analogous to the one in main program

  factor_jrkr = (jr*(jr+one)-kr*(kr+one))

  if(oner_flag==0) then

    do i =1,idx3(kr)%size
      vec32 => eigvec3(i)%vec(:)
      n= idx3(kr)%start+i-1

      ! write(6,*) ' test energy:',n, eigvec3(i)%eng
      h7(xn(n,n)) = eigvec3(i)%eng

      do i1 =1,i
        vec31 => eigvec3(i1)%vec(:)
        n1 = idx3(kr)%start+i1-1

        temp = 0.0D0
        do id=1,nn(3)
          if(idx2(id)%size < 1) cycle

          temp = temp+ oner3(id,id)*dot_product( &
            vec31(idx2(id)%start:idx2(id)%end),  &
            vec32(idx2(id)%start:idx2(id)%end))

        end do

        !  to test the influence of <K| |K> sub-block comment it out
        h7(xn(n1,n)) = h7(xn(n1,n))+temp* factor_jrkr

      end do
    end do

  else

    rewind(u16)
    ! allocate(h6(istage2,istage2))
    ! that's a column
    allocate(h6(istage2))

    do id=1,nn(3)
      write(6,*) ' id=', id
      ! if idx2(id)%size < 1 there's nothing to do
      if(idx2(id)%size < 1) cycle

      ! for efficient computation we need to precompute off-diagonal
      ! in "id" vector-vector products

      ! allocate auxiliary array of matrices to hold vector products
      if(id > 1) then
        ! write(6,*) '   preparing aux matrices'
        allocate(vm(nn(1),nn(2),1:id-1))

        do id1=1,id-1

          ! if(oner_flag==0 .AND. id1 /= id) cycle

          if(idx2(id1)%size < 1) cycle
          do ic=1,nn(2)
            do ib=1,nn(1)
              if(idx1(ib,ic,id)%size < 1 .OR. idx1(ib,ic,id1)%size < 1) cycle

              allocate(vm(ib,ic,id1)%mat(idx1(ib,ic,id1)%size, &
                                       idx1(ib,ic,id)%size))
              do ia=1,idx1(ib,ic,id)%size
                vec12 => eigvec1(rec_num1(ia,ib,ic,id))%vec(:)

                do ia1=1,idx1(ib,ic,id1)%size
                  vec11 => eigvec1(rec_num1(ia1,ib,ic,id1))%vec(:)

                  vm(ib,ic,id1)%mat(ia1,ia)= dot_product(vec11,vec12)
                end do
              end do

            end do
          end do

        end do

        ! write(6,*) '   finished preparing aux matrices'
      end if

      do im=1,idx2(id)%size
        n = idx2(id)%start+im-1

        h6(:) = zero

          vec22 => eigvec2(rec_num2(im,id))%vec(:)

        ! diagonal contribution

        h6(n) = oner3(id,id)

        ! off-diagonal (in id) contributions

        do id1=1,id-1

          ! if(oner_flag==0 .AND. id1 /= id) cycle

          do im1=1,idx2(id1)%size
            n1 = idx2(id1)%start+im1-1

              vec21 => eigvec2(rec_num2(im1,id1))%vec(:)

            temp = 0.0D0
            do ic=1,nn(2)
              do ib=1,nn(1)
                if(idx1(ib,ic,id)%size < 1 .OR. idx1(ib,ic,id1)%size < 1) cycle

                allocate(e6(idx1(ib,ic,id)%size))

                e6 = matmul(  &
                  vec21(idx1(ib,ic,id1)%start:idx1(ib,ic,id1)%end), &
                  vm(ib,ic,id1)%mat(:,:))
                temp = temp+ dot_product(e6, &
                  vec22(idx1(ib,ic,id)%start:idx1(ib,ic,id)%end))

                deallocate(e6)
              end do
            end do

            h6(n1) = temp* oner3(id1,id)

          end do
        end do

        !  store the composed column of h6 to disk
        ! write(u16) h6(1:n)
        ! write(6,*) ' writing record=', n
        write(u16) h6(:)

      end do

      !  deallocate auxiliary array of matrices

      if(id > 1) then
        do id1=1,id-1

          ! if(oner_flag==0 .AND. id1 /= id) cycle

          if(idx2(id1)%size < 1) cycle
          do ic=1,nn(2)
            do ib=1,nn(1)
              if(idx1(ib,ic,id)%size < 1 .OR. idx1(ib,ic,id1)%size < 1) cycle

              deallocate(vm(ib,ic,id1)%mat)
            end do
          end do

        end do
        deallocate(vm)
      end if

      ! write(6,*) '  id end'
    end do

    deallocate(h6)

    nullify(vec22,vec21,vec11,vec12)

    !  stage2 -> stage3

    allocate(h5(istage2,istage2),e5(istage2))
    rewind(u16)
    do n=1,istage2
      read(u16) h5(:,n)
    end do

    !  however the stored h5 holds only half of the matrix
    !  fill in the missing part
    do n=1,istage2
      do n1=1,n-1
        h5(n,n1)=h5(n1,n)
      end do
    end do

    do i=1,idx3(kr)%size
      vec32 => eigvec3(i)%vec(:)
      n = idx3(kr)%start+i-1

      h7(xn(n,n)) = eigvec3(i)%eng

      e5 = matmul(h5(:,:),vec32(:))

      do i1=1,i
        vec31 => eigvec3(i1)%vec(:)
        n1 = idx3(kr)%start+i1-1

        !  to test the influence of <K| |K> sub-block comment it out
        h7(xn(n1,n))= h7(xn(n1,n))+dot_product(vec31(:),e5(:))* factor_jrkr
      end do
    end do

    deallocate(e5,h5)

  end if

  nullify(vec32,vec31)

end do kr_loop

!  deallocate eigenvectors and arrays

do n=1,istage1;  deallocate(eigvec1(n)%eng,eigvec1(n)%vec);  end do
do n=1,istage2;  deallocate(eigvec2(n)%eng,eigvec2(n)%vec);  end do
do n=1,istage3;  deallocate(eigvec3(n)%eng,eigvec3(n)%vec);  end do

deallocate(eigvec1,eigvec2,eigvec3,rec_num1,rec_num2, idx1A,idx2A)

write(u6,'(/," h7 done. diagonalizing... ")')

nlev = idx3(jr)%end

allocate(vectors(nlev,min(nlev,abs(icut4))))
allocate(e7(nlev))
call diag_p(h7,e7,vectors,nlev,nout,icut4,encut4,4,idiag)
! call la_spevx(h7,e7,UPLO='U',VU=encut4,M=nout,INFO=idiag)

if(idiag /= 0) then
  write(6,*) ' mainj: diag_p failed to diagonalise:',idiag
  stop
end if

write(u6,'(/," FINALLY: h7 energies (relative to zero energy):")')

do icounter=1,nout
  ! print out levels
  write(u6,'(2x,i4,f24.6)') icounter,e7(icounter)-enzero
  ! save vectors
  write(u41(jr)) e7(icounter),vectors(:,icounter)
end do

deallocate(vectors)
deallocate(h7,e7)

! CONTAINS

END SUBROUTINE
