! Transform wave functions to primitive basis in block fashion
!
! Details: the main idea is to increase efficiency by avoiding vector-matrix
! multiplication and using matrix-matrix multiplication instead.
! All eigenfunctions from all stages are read into memory but
! processing proceeds separately for every radial grid point.
! For each radial grid point the wavefunctions are converted
! to the primitive basis and written to disk at once for all eigenstates.
! If HHSYM is True then symmetric and asymmetric states are processed
! simultaneously.
!
! The format of created files is UNFORMATTED. The structure is the following
! nstates - number of states
! e3n(:)  - energies
! loop over K (not implemented yet)
!    loop over all radial grid points and write angular basis coefficients
!       do i=1,nstates
!          write(u09) pbasis(:,i)
!       enddo
!
!
! NOTE: this is for K=0 only!!
! ASSUMED jmax >= lmax


PROGRAM t2pb
  USE types
  USE param
  USE base_lib
  USE angles
  USE workarrays
IMPLICIT NONE

! stage 3 vectors and energies: non-sym, sym & asy
REAL(DP), ALLOCATABLE :: v3n(:,:), v3s(:,:), v3a(:,:)
REAL(DP), ALLOCATABLE :: e3n(:),   e3s(:),   e3a(:)
! intermediate temporary matrices
REAL(DP), ALLOCATABLE :: mtmp(:,:),h6(:,:)
!REAL(DP), ALLOCATABLE, TARGET :: h3(:,:)
REAL(DP), POINTER :: h3(:,:),h5(:,:)
! primitive basis in each i1,i2,3 grid point
REAL(DP), ALLOCATABLE :: pbasis(:,:)
! J wave-functions
REAL(DP), ALLOCATABLE :: jwfns(:,:)

!debug! 
REAL(DP), ALLOCATABLE :: fullv(:,:,:,:,:)


!TYPE(triple), ALLOCATABLE :: idx1(:,:,:),idx2(:)

INTEGER(I4B) :: istage1,istage2,istage3
INTEGER(I4B) :: ist1_max,ist2_max
INTEGER(I4B) :: idx1size_max,idx2size_max
!
INTEGER(I4B) :: n, ii,ni,i1,i2,i3
INTEGER(I4B) :: ia,ib,ic,id,im, m, ik
INTEGER(I4B) :: i, icounter, idiag, nsize, icut, istart
INTEGER(I4B) :: j, l, k, kr ! m, kp
INTEGER(I4B) :: j1, l1, k1  ! m1, kp1
INTEGER(I4B) :: nlev,nnmax,istage2s,istage2a,itmp
! INTEGER(I4B) :: j2min,j2max,l2min,l2max,j2sb,l2sb
! INTEGER(I4B) :: p2,k2,xk1,xk
!
!  parameters
!
! number of stage 3 states: non-sym, sym & asy
INTEGER(I4B) :: nstates, nss, nsa
!
REAL(DP) :: zpe, tmp,temp,tmp1,temp1, energy, x, step, gmf, gmf2
REAL(DP) :: one_over_sqrt2, ecut

INTEGER(I4B) :: u09=9,u11=11,u12=12,u13=13,u14=14
INTEGER(I4B) :: jwf_len, jwf_num, idx_k_start, idx_k_end, idx_k_size

! not used; it's only here for the sake of reading legacy input.txt
INTEGER(I4B) :: iauto

CHARACTER(LEN=80) title
CHARACTER(9) :: fname
!fname= 'x0000-K0-'

one_over_sqrt2 = one/sqrt(two)

! set the window for final eigenstates from state # istart to # icut 
! or Energy ecut (relative to zero energy) 
! - not implemented!
! istart = 0; icut =  10; ecut = 10000.D0

!open(6,file='out.txt')

! ====================================================================
! INPUT (part1)

open(5,file='input.txt',status='old')

write(u6,'(" Transform RV Wave-Functions to Primitive Basis",/, &
           & " INPUT: reading input...",/, &
           & 30("="),"< file:  input.txt >",30("="))')
read(u5,'(A80)') title; write(u6,'(A80)') title
!  masses:
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) mass1,mass2,mass3,mass4
write(u6,'(4(2x,f13.9))') mass1,mass2,mass3,mass4
!  coordinate system
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) opt
write(u6,'(i5)') opt
!  parameters of stretching basis functions  and the grid size
read(u5,'(A80)') title; write(u6,'(A80)') title
do i=1,3
  read(u5,*) igq(i),re(i), we(i), De(i), nn(i)
  write(u6,'(i4,f8.4,f12.3,f12.3,i6)') igq(i),re(i),we(i),De(i),nn(i)
end do
!  flags:
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) angular_problem_only_flag, optimize_flag, test_flag, expansion_flag
write(u6,'(11x,i1,23x,i1,13x,i1,13x,i1)') &
  angular_problem_only_flag, optimize_flag, test_flag, expansion_flag
!  angular basis functions
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) jmax,lmax,kmax, jrmax, krmax, &
  j_parity_max, l_parity_max, jl_parity_max
write(u6,'(1x,8i5)') jmax,lmax,kmax, jrmax, krmax, &
  j_parity_max, l_parity_max, jl_parity_max
!  expansion size
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) ne1, ne2, ne3
write(u6,'(1x,3i5)') ne1, ne2, ne3
!  angular grid points size
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) nt1, nt, nphi, iauto
write(u6,'(1x,4i5)') nt1, nt, nphi, iauto
!  zero energy:
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) enzero
write(u6,'(1x,f20.10)')  enzero
write(u6,'(80("="))')

close(5)

! ====================================================================
! get namax & mmax
! WARNING: angular_states_max  returns j_parity, l_parity, jl_parity
!          this is a by product to be aware of

call angular_states_max()

! ====================================================================
! INPUT (part2)

open(5,file='inp.transf.txt',status='old')

read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) jp, j_parity, l_parity, jl_parity
write(u6,'(1x,4i5)') jp, j_parity, l_parity, jl_parity

read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) jr
write(u6,'(1x,i5)') jr

! define missing quantum numbers
!jr=0; kr=0; 
p=jp-jr

close(5)

! ====================================================================

allocate(indx(namax,6),krn(0:krmax))

allocate(idx1(nn(1),nn(2),nn(3)),idx2(nn(3)))
allocate(vectrs1(nn(1),nn(2),nn(3)),vectrs2(nn(3)))

allocate(angsym(0:krmax,namax))

call angular_states_qnumbers(indx,krn)

! ====================================================================
!
! 27/04/2009
! two pass processing:
! 1 pass: transform w-f to primitive basis for each K
!         output - x0000-K?-all.dat (ie compatible with what was before)
! 2 pass: transform w-f to primitive basis for each J
!         output - x0000-J?-all.dat 
!         x0000-K?-all.dat can be re-read multiple times

! ====================================================================
! 1st pass

! loop over K
kr_loop: do kr=0,min(jr,krmax)

      ! if K>0 then eigenvectors for p=0 and p=1 are identical (by construction)
      ! K>0 jp=1 wave-functions are not computed and hence they are transformed

      n=jp
      ! the line below allows to re-use K>0 p=0 eigenvectors obtained for p=1
      if(kr>0) n=0
      ! then n needs to be used instead of jp when opening files

      ! ===================================================================
      ! define the prefix for the names of the files which store wave-functions (input)
      fname = 'x'//i2c(n)//i2c(j_parity)//i2c(l_parity)//i2c(jl_parity)//'-K'//i2c(kr)//'-'

      write(u6,'(80("="))')
      write(u6,*) ' transforming K=',kr,' wave-functions to primitive basis'
      write(u6,*) ' implied prefix for wave-function files: ', fname

      ! output
        open(u10,file=fname//'all.dat',status='OLD', &
          access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
      if(idiag==0) then
        write(u6,*) ' using existing file ',fname//'all.dat'
        cycle kr_loop
      else
        open(u10,file=fname//'all.dat', &
          access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
      end if
      write(u6,*) '  opened ',fname//'all.dat',' for output'

      ! input
      open(u11,file=fname//'1.dat', &
        status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
      if(idiag /=0) then 
        write(6,*) '  cant open ',fname//'1.dat';  stop
      end if
      write(u6,*) '  opened ',fname//'1.dat',' for input'

      open(u12,file=fname//'2.dat', &
        status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
      if(idiag /=0) then
        write(6,*) '  cant open ',fname//'2.dat';  stop
      end if
      write(u6,*) '  opened ',fname//'2.dat',' for input'

      open(u13,file=fname//'3.dat', &
        status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
      if(idiag /=0) then
        write(6,*) '  cant open ',fname//'3.dat'; stop
      end if
      write(u6,*) '  opened ',fname//'3.dat',' for input'

      ! ====================================================================
      ! stage 1

      write(u6,'("  reading in stage 1 eigenvectors...")')

      read(u11) istage1, ist1_max
      write(6,*) ' istage1=', istage1

      ! get indexing array
      read(u11) idx1

      !do i2=1,nn(2); do i1=1,nn(1); do i3=1,nn(3)
      !      write(u6,'(3i3,3i12)') i1,i2,i3, &
      !        idx1(i1,i2,i3)%size,idx1(i1,i2,i3)%start,idx1(i1,i2,i3)%end
      !end do; end do; end do

      ! the size of angular basis for each K
      na = krn(kr)%size

      allocate(e3n(na))

      do i3=1,nn(3); do i2=1,nn(2); do i1=1,nn(1)

          if(idx1(i1,i2,i3)%size < 1) cycle
          !!!if(HHSYM .and. i1 > i2) cycle 

          allocate(vectrs1(i1,i2,i3)%mat(1:na,1:idx1(i1,i2,i3)%size))

          do ia=1,idx1(i1,i2,i3)%size
            read(u11) e3n(ia), vectrs1(i1,i2,i3)%mat(:,ia)
          end do

          ! orthogonality check
          !do ia=1,idx1(i1,i2,i3)%size
          !  tmp = 0.0D0
          !  do i=1,na
          !    tmp = tmp + vectrs1(i1,i2,i3)%mat(i,ia)**2
          !  end do
          !  write(6,'(4i3,f14.10)') ia,i1,i2,i3, tmp
          !end do

      end do; end do; end do

      close(u11)

      deallocate(e3n)

      !  update stage 1 indices, count the total number of stage1 vectors and
      !  find max number of selected levels per grid

      call setup_idx1()

      ! ===================================================================
      ! stage 2

      write(u6,'("  reading in stage 2 eigenvectors...")')

      read(u12) istage2,ist2_max
      write(6,*) ' istage2=', istage2

      ! get indexing array
      read(u12) idx2

      do id =1,nn(3)

          if(idx2(id)%size < 1) cycle

          nnmax = idx1(nn(1),nn(2),id)%end
          if(nnmax < 1) cycle

          !  number of stage 2 eigenvectors at (id) grid 
          allocate(e3n(1:idx2(id)%size))
          allocate(vectrs2(id)%mat(1:nnmax,1:idx2(id)%size))

          do n=1,idx2(id)%size
            read(u12) e3n(n), vectrs2(id)%mat(:,n)
          end do

          deallocate(e3n)

          ! orthogonality check
          !do n=1,idx2(id)%size
          !  tmp = 0.0D0
          !  do i=1,idx1(nn(1),nn(2),id)%end
          !    tmp = tmp + vectrs2(id)%mat(i,n)**2
          !  end do
          !  write(6,'(2i3,f14.10)') n,id, tmp
          !end do

      end do

      close(u12)

      !  update stage 2 indices, count the total number of stage2 vectors and
      !  find max number of selected levels per grid

      call setup_idx2()

      ! ===================================================================
      ! stage 3

      write(6,'(" reading stage 3...")')

      istage2 = idx2(nn(3))%end


      !  if(ii < istart) cycle
      !  if(energy-zpe > ecut) cycle
      !  if(ii > icut) cycle

      !  c(i1,i2,q;t) c(i3,p;q) c(i,p)
      !  v3           vectrs2   vectrs1
      !  h6           h5        h3
      !  accumulator for (i1,i2,i3,K)  

      !!!if(HHSYM) then ! ====================================================
      !!!  <currently empty>
      !!!else ! HHSYM=.FALSE. ================================================

      write(u6,*) ' istage2=', istage2
      read(u13) nstates,itmp
      ! itmp above is effectively istage2 and must be equal to it

      ! allocate storage for eigen vectors and energies
      allocate(v3n(istage2,nstates),e3n(nstates))

      !debug! 
      !!allocate(fullv(na,nn(1),nn(2),nn(3),nstates))


      do ii=1,nstates
        read(u13) e3n(ii),v3n(:,ii)
      end do

      write(u10) nstates
      write(u6,*) ' nstates=', nstates
      write(u10) e3n(:)

      ! find max size of the matrices which are needed for transformation to truncated basis
      idx1size_max = 0
      do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
        if(idx1size_max < idx1(ib,ic,id)%size) idx1size_max = idx1(ib,ic,id)%size
      enddo; enddo; enddo
      allocate(h3(na,idx1size_max))

      idx2size_max = 0
      do id=1,nn(3)
        if(idx2size_max < idx2(id)%size) idx2size_max = idx2(id)%size
      enddo
      !!allocate(h5(idx1size_max,idx2size_max))
      allocate(h6(idx2size_max,nstates),mtmp(na,idx2size_max),pbasis(na,nstates))

      do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
      !  replace the nested loops above by one loop
      ! merged_loop do i=1,nn(3)*nn(2)*nn(1)
      
        !  unroll true indices if using merged loop
        !id = int(int((i-1)/nn(1))/nn(2))+1
        !ic = int((i-1)/nn(1))-nn(2)*(id-1)+1
        !ib = i-nn(1)*(ic-1)-nn(1)*nn(2)*(id-1)

        ! h3 => vectrs1(ib,ic,id)%mat(:,:)
        do ia =1,idx1(ib,ic,id)%size
          do ni=1,na
            h3(ni,ia) =  vectrs1(ib,ic,id)%mat(ni,ia)
          end do
        end do

        h5 => vectrs2(id)%mat(idx1(ib,ic,id)%start:idx1(ib,ic,id)%end,:)
        ! the pointer and the loop (below) are equivalent
        ! note: h5 needs to be allocated first
        !!do ik=1,idx2(id)%size
        !!  do im=1,idx1(ib,ic,id)%size
        !!    n = idx1(ib,ic,id)%start+im-1
        !!    h5(im,ik) = vectrs2(id)%mat(n,ik)
        !!  end do
        !!end do

      
        ! h6 => v3n(idx2(id)%start:idx2(id)%end,:)
        do ik=1,nstates
          do im=1,idx2(id)%size
            n = idx2(id)%start+im-1
            h6(im,ik) = v3n(n,ik)
          end do
        end do

        ! compute primitive basis (aggregation on the radial grid ib,ic,id)
        ! note that the order of the product is important and h3 x (h5 x h6) is wrong
        !  
        !  (h3 x h5) x h6 

        mtmp(1:na,1:idx2(id)%size) = &
      &   matmul(h3(1:na,1:idx1(ib,ic,id)%size),h5(1:idx1(ib,ic,id)%size,1:idx2(id)%size))

        pbasis(1:na,1:nstates) = &
      &   matmul(mtmp(1:na,1:idx2(id)%size),h6(1:idx2(id)%size,1:nstates))

        do ii=1,nstates
          write(u10) pbasis(:,ii)
          !debug! check the norm of the resulting w-f
          !!fullv(:,ib,ic,id,ii)=pbasis(:,ii)
        enddo
      
      !end do merged_loop
      end do; end do; end do

      !debug! check the norm of the resulting w-f
      !!do ii=1,nstates
      !!  tmp = 0.0D0
      !!  do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1); do ia=1,na
      !!    tmp = tmp + fullv(ia,ib,ic,id,ii)**2
      !!  end do; end do; end do; end do
      !!  write(6,*) ii, tmp
      !!end do

      !debug! check orthogonality of the resulting w-f
      !!do ii=1,nstates
      !!do ni=1,nstates
      !!  tmp = 0.0D0
      !!  do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1); do ia=1,na
      !!    tmp = tmp + fullv(ia,ib,ic,id,ii)*fullv(ia,ib,ic,id,ni)
      !!  end do; end do; end do; end do
      !!  if(abs(tmp) > 1.D-15) write(6,*) ii,ni, tmp
      !!end do
      !!end do

      !debug!
      !!deallocate(fullv)

      deallocate(h3,h6,mtmp,pbasis)
      deallocate(v3n,e3n)
      do i3=1,nn(3)
        if(associated(vectrs2(i3)%mat)) &
         deallocate(vectrs2(i3)%mat)
      end do
      
      do i3=1,nn(3); do i2=1,nn(2); do i1=1,nn(1)
        if(associated(vectrs1(i1,i2,i3)%mat)) &
         deallocate(vectrs1(i1,i2,i3)%mat)
      end do; end do; end do

      write(u6,'(" finished K=",i2)') kr

      !!!endif

end do kr_loop

      
close(u10)
!!!if(HHSYM) close(u09)

!!!if(HHSYM) deallocate(s1sym,s2sym)
deallocate(angsym)

! deallocate vectors

deallocate(idx1,idx2)
deallocate(vectrs1,vectrs2)

      write(u6,'(80("="))')

! ====================================================================
! 2nd pass

if(jr==0) then
  ! nothing more to do
  call time_stamp(6)
  stop
end if

! define the name of the file to store primitive wave-functions (output)
fname = 'x'//i2c(jp)//i2c(j_parity)//i2c(l_parity)//i2c(jl_parity)//'-J'//i2c(jr)//'-'

write(u6,*) ' transforming J=',jr,' wave-functions to primitive basis'
write(u6,*) ' prefix for output wave-function file: ', fname

open(u09,file=fname//'all.dat', &
        access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
write(u6,*) '  opened ',fname//'all.dat',' for output'


! define the name of the file to read wave-functions
fname = 'x'//i2c(jp)//i2c(j_parity)//i2c(l_parity)//i2c(jl_parity)//'-K'//i2c(jr)//'-'

! input: J wave-functions in contructed basis
open(u14,file=fname//'4.dat', &
  status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
if(idiag /=0) then 
  write(6,*) '  cant open ',fname//'4.dat';  stop
end if
write(u6,*) '  opened ',fname//'4.dat',' for input'

! J wave functions: the number of them and their vector length
read(u14) jwf_num, jwf_len

! read J wave-functions and energies

allocate(e3n(1:jwf_num))

allocate(jwfns(jwf_len,jwf_num))
do i=1,jwf_num
  read(u14) e3n(i), jwfns(:,i)
end do

close(u14)
write(u6,*) '  closed ',fname//'4.dat'

! write to new x????-J?-all.dat w-f file in a format compatible with x????-K?-all.dat
write(u09) jwf_num
write(u09) e3n(:)

deallocate(e3n)


! indexing construct similar to idx1 and idx2
idx_k_start = 0
idx_k_end   = 0
idx_k_size  = 0

write(u6,*) '          K  idx_k_size idx_k_start   idx_k_end'

do kr=0,jr

      ! the size of angular basis for each K
      na = krn(kr)%size

      ! the part of J w-functions relevant to current K is
      !jwfns(idx_k_start:idx_k_end,:)


      n=jp
      ! the line below allows to re-use K>0 p=0 eigenvectors obtained for p=1
      if(kr>0) n=0
      ! then n needs to be used instead of jp when opening files

      ! define the name of the file to read wave-functions
      fname = 'x'//i2c(n)//i2c(j_parity)//i2c(l_parity)//i2c(jl_parity)//'-K'//i2c(kr)//'-'

      ! input: K wave-functions in primitive basis
      open(u10,file=fname//'all.dat', &
        access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
      if(idiag /=0) then 
        write(u6,*) '  cant open ',fname//'all.dat';  stop
      end if
        write(u6,*) '  opened ',fname//'all.dat'

      ! we assume that nstates > 0 (otherwise the file would be not there)
      read(u10) nstates
      allocate(e3n(nstates))
      read(u10) e3n(:)
      deallocate(e3n)

      idx_k_size  = nstates
      idx_k_start = idx_k_end+1
      idx_k_end   = idx_k_start+idx_k_size-1

      write(u6,*) kr, idx_k_size, idx_k_start, idx_k_end

      ! the part of J w-functions relevant to current K is
      !jwfns(idx_k_start:idx_k_end,1:jwf_num)

      ! allocate temp buffer for reading primitive wave-functions for current K
      allocate(mtmp(na,nstates))

      ! allocate the buffer for resulting primitive wave-functions
      allocate(pbasis(na,jwf_num))

      do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
         do i=1,nstates
           read(u10) mtmp(:,i)
         end do

         pbasis(1:na,1:jwf_num) = matmul(mtmp(1:na,1:nstates),jwfns(idx_k_start:idx_k_end,1:jwf_num))
         do i=1,jwf_num
           write(u09) pbasis(:,i)
         end do
      end do; end do; end do

      deallocate(mtmp,pbasis)
      close(u10)
        write(u6,*) '  closed ',fname//'all.dat'

end do

write(u6,*) ' jwf_len  =', jwf_len, ' must be equal to the last idx_k_end'

deallocate(jwfns)
deallocate(indx,krn)
close(u09)

call time_stamp(6)

! end
! ====================================================================

CONTAINS

FUNCTION i2c(i)
  USE types
IMPLICIT NONE
INTEGER(I4B) :: i
CHARACTER*1 i2c

if(i<0 .OR. i>9) then
 write(6,*) i
 stop ' integer is out of range'
else
  i2c=char(48+i)
end if
END FUNCTION


SUBROUTINE setup_idx1()
  USE types
IMPLICIT NONE
INTEGER(I4B) :: n,i1,i2,i3

      ! input:  idx1(i1,i2,i3)%size - the number of selected angular functions
      ! output: idx1(i1,i2,i3)%start 
      !         idx1(i1,i2,i3)%end  - define bounds for a running index 
      !                               for all i1, i2, and given i3

      ! given i3, find the running index start and end positions
      ! as a function of i2 & i1 based on the length of selected
      ! angular basis functions
  
      do i3=1,nn(3)
        n=0
        do i2=1,nn(2)
          do i1=1,nn(1)
            if(idx1(i1,i2,i3)%size==0) then
              idx1(i1,i2,i3)%start=n
              idx1(i1,i2,i3)%end=n
              !write(u6,'(3i5,4i12)') i3,i2,i1, namax, &
              !  idx1(i1,i2,i3)%size,idx1(i1,i2,i3)%start,idx1(i1,i2,i3)%end
            else
              idx1(i1,i2,i3)%start=n+1
              n=n+idx1(i1,i2,i3)%size
              idx1(i1,i2,i3)%end=n
              !write(u6,'(3i5,4i12)') i3,i2,i1, namax, &
              !  idx1(i1,i2,i3)%size,idx1(i1,i2,i3)%start,idx1(i1,i2,i3)%end
            end if
          end do
        end do

      end do
END SUBROUTINE


SUBROUTINE setup_idx2()
  USE types
IMPLICIT NONE
INTEGER(I4B) :: n,id

      ! input:  idx2(id)%size - the number of selected functions from stage1 
      !                         for all ib, ic, and given id
      ! output: idx2(id)%start
      !         idx2(id)%end  - define bounds for a running index for given id


      !write(u6,'(" id  full basis  red. basis   start indx  end indx")')
      n=0
      do id=1,nn(3)
        if(idx2(id)%size==0) then
           idx2(id)%start=n
           idx2(id)%end=n
        else
           idx2(id)%start=n+1
           n=n+idx2(id)%size
           idx2(id)%end=n
        end if
        !write(u6,'(i3,4i12)') id, idx1(nn(1),nn(2),id)%end, &
        !  idx2(id)%size,idx2(id)%start,idx2(id)%end
      end do
END SUBROUTINE


END PROGRAM


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
