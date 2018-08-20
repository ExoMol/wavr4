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

!debug! REAL(DP), ALLOCATABLE :: fullv(:,:,:,:,:)


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

INTEGER(I4B) :: u09=9,u11=11,u12=12,u13=13

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

! define missing quantum numbers
jr=0; kr=0; p=jp-jr

fname = 'x'//i2c(jp)//i2c(j_parity)//i2c(l_parity)//i2c(jl_parity)//'-K'//i2c(kr)//'-'

write(u6,*) ' implied prefix for the file with wave-functions: ', fname
write(u6,'(80("="))')

close(5)

! ====================================================================

allocate(indx(namax,6),krn(0:krmax))

allocate(idx1(nn(1),nn(2),nn(3)),idx2(nn(1),nn(2)))
allocate(vectrs1(nn(1),nn(2),nn(3)),vectrs2(nn(1),nn(2)))

if(HHSYM) allocate(idx2sym(nn(1),nn(2)),idx2asy(nn(1),nn(2)))

allocate(angsym(0:krmax,namax))
if(HHSYM) then
  allocate( s1sym(nn(1),nn(3),0:krmax,namax), &
            s2sym(nn(1),0:krmax,nn(3)*namax))
end if

call angular_states_qnumbers(indx,krn)

! the size of angular basis for each K
na = krn(kr)%size

! ===================================================================

open(u11,file=fname//'1.dat', &
  status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
if(idiag /=0) then 
  write(6,*) '  cant open ',fname//'1.dat';  stop
end if

open(u12,file=fname//'2.dat', &
  status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
if(idiag /=0) then
  write(6,*) '  cant open ',fname//'2.dat';  stop
end if

open(u13,file=fname//'3.dat', &
  status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
if(idiag /=0) then
  write(6,*) '  cant open ',fname//'3.dat'; stop
end if

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

      allocate(e3n(na))

      do i3=1,nn(3); do i2=1,nn(2); do i1=1,nn(1)

          if(idx1(i1,i2,i3)%size < 1) cycle
          if(HHSYM .and. i1 > i2) cycle 

          read(u11) e3n(1:idx1(i1,i2,i3)%size)

          allocate(vectrs1(i1,i2,i3)%mat(1:na,1:idx1(i1,i2,i3)%size))

          do ia=1,idx1(i1,i2,i3)%size
            read(u11) vectrs1(i1,i2,i3)%mat(:,ia)
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
      !  this also takes care of H-H symmetry if present

      call setup_idx1(kr,istage1,ist1_max)

! ===================================================================
! stage 2

      write(u6,'("  reading in stage 2 eigenvectors...")')

      read(u12) istage2,ist2_max
      write(6,*) ' istage2=', istage2

      ! get indexing array
      read(u12) idx2

      do ic =1,nn(2); do ib =1,nn(1)

          if(idx2(ib,ic)%size < 1) cycle
          if(HHSYM .and. ib > ic) cycle 

          !  length of stage 2 eigenvectors at (ib,ic) grid
          !  1:nn(3) id grid points x reduced angular basis at each (ib,ic,id) grid point
          !!nnmax = 0
          !!do n=1,nn(3)
          !!  nnmax = nnmax+idx1(ib,ic,n)%size
          !!end do
          !  the length can be easily computed as above or just taken as
          nnmax = idx1(ib,ic,nn(3))%end
          if(nnmax < 1) cycle

          !  number of stage 2 eigenvectors at (ib,ic) grid
          allocate(e3n(1:idx2(ib,ic)%size))
          allocate(vectrs2(ib,ic)%mat(1:nnmax,1:idx2(ib,ic)%size))

          read(u12) e3n(1:idx2(ib,ic)%size)

          do n=1,idx2(ib,ic)%size
            read(u12) vectrs2(ib,ic)%mat(:,n)
          end do

          deallocate(e3n)

          ! orthogonality check
          !do n=1,idx2(ib,ic)%size
          !  tmp = 0.0D0
          !  do i=1,idx1(ib,ic,nn(3))%end
          !    tmp = tmp + vectrs2(ib,ic)%mat(i,n)**2
          !  end do
          !  write(6,'(3i3,f14.10)') n,ib,ic, tmp
          !end do

      end do; end do

      close(u12)

      !  update stage 2 indices, count the total number of stage2 vectors and
      !  find max number of selected levels per grid
      !  this also takes care of H-H symmetry if present

      call setup_idx2(kr,istage2,ist2_max)

! ===================================================================
! stage 3

      write(6,'(" reading stage 3...")')

      istage2 = idx2(nn(1),nn(2))%end


      !  if(ii < istart) cycle
      !  if(energy-zpe > ecut) cycle
      !  if(ii > icut) cycle

!  c(i1,i2,q;t) c(i3,p;q) c(i,p)
!  v3           vectrs2   vectrs1
!  h6           h5        h3
!  accumulator for (i1,i2,i3,K)  

if(HHSYM) then ! ====================================================

      istage2s = idx2sym(nn(1),nn(2))%end
      istage2a = idx2asy(nn(1),nn(2))%end

      ! action of P(HH) symmetry on angular basis:
      ! P(HH) |i> = (-1)^p |i'>  if K=0
      ! P(HH) |i> = (-1)^K |i'>  if K>0

      ! the sign factor
      gmf = one

      if (jp > 0) gmf = -one
      if (kr > 0 .AND. mod(kr,2) > 0 ) gmf = -one

      write(u6,*) ' istage2s=', istage2s
      ! read the number of symmetric states
      read(u13) nss,itmp

      ! allocate storage for eigen vectors and energies
      allocate(v3s(istage2s,nss),e3s(nss))

      do ii=1,nss
        read(u13) e3s(ii),v3s(:,ii)
      end do

      write(u6,*) ' istage2a=', istage2a
      ! read the number of asymmetric states
      read(u13) nsa,itmp

      ! allocate storage for eigen vectors and energies
      allocate(v3a(istage2a,nsa),e3a(nsa))

      do ii=1,nsa
        read(u13) e3a(ii),v3a(:,ii)
      end do


      ! if(nstates > 0)  
      open(u10,file=fname//'sym.dat', &
        access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
      open(u09,file=fname//'asy.dat', &
        access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)

      write(u10) nss
      write(u10) e3s(:)
      write(u09) nsa
      write(u09) e3a(:)

      ! find max size of the matrices which transform to truncated basis
      idx1size_max = 0
      do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
        if(idx1size_max < idx1(ib,ic,id)%size) idx1size_max = idx1(ib,ic,id)%size
      enddo; enddo; enddo
      allocate(h3(na,idx1size_max))

      idx2size_max = 0
      do ic=1,nn(2); do ib=1,nn(1)
        if(idx2size_max < idx2(ib,ic)%size) idx2size_max = idx2(ib,ic)%size
        ! if(idx2size_max < idx2sym(ib,ic)%size) idx2size_max = idx2sym(ib,ic)%size
        ! if(idx2size_max < idx2asy(ib,ic)%size) idx2size_max = idx2asy(ib,ic)%size
      enddo; enddo
      allocate(mtmp(na,idx2size_max),h6(idx2size_max,max(nss,nsa)),pbasis(na,max(nss,nsa)))

        !write(6,*) na, idx1size_max,":",nnmax, idx2size_max, max(nss,nsa)
        !write(6,*) "h6   :", idx2size_max, " x ", max(nss,nsa)
        !write(6,*) "pbasis", na, " x ", max(nss,nsa)

      ! do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
      !  replace the nested loops above by one loop
      do i=1,nn(3)*nn(2)*nn(1)
      
        !  unroll true indices
        id = int(int((i-1)/nn(1))/nn(2))+1
        ic = int((i-1)/nn(1))-nn(2)*(id-1)+1
        ib = i-nn(1)*(ic-1)-nn(1)*nn(2)*(id-1)

        !write(6,*) ib,ic,id

        ! symmetry factor
        gmf2 = one
        if(ib /= ic) gmf2 = gmf2* one_over_sqrt2
        if(ib > ic)  gmf2 = gmf2* gmf
      
        ! matrix of stage 2 vectors
        if(ib <= ic) then
          h5 => vectrs2(ib,ic)%mat(idx1(ib,ic,id)%start:idx1(ib,ic,id)%end,:)
        else ! ib > ic
          h5 => vectrs2(ic,ib)%mat(idx1(ic,ib,id)%start:idx1(ic,ib,id)%end,:)
        end if

        ! matrix of stage 1 vectors
        if(ib <= ic) then
          do ia =1,idx1(ib,ic,id)%size
            do ni=1,na
              h3(ni,ia) =  vectrs1(ib,ic,id)%mat(ni,ia)
            end do
          end do
        else ! ib > ic
          do ia =1,idx1(ic,ib,id)%size
            do ni=1,na
              h3(ni,ia) =  vectrs1(ic,ib,id)%mat(angsym(kr,ni),ia)
            end do
          end do
        end if

        ! matrix of stage 3 vectors - sym case
        h6(:,:) = zero
        if(ib /= ic) then
          do ik=1,nss
            do im=1,idx2sym(ib,ic)%size
              n = idx2sym(ib,ic)%start+im-1
              h6(im,ik) = v3s(n,ik)* gmf2
            end do
          end do
        else ! ib==ic
          do ik=1,nss
            im=0
            do ii=1,idx2(ib,ic)%size
              if(s2sym(ib,kr,ii)) then
                im=im+1
                n = idx2sym(ib,ic)%start+im-1
                h6(ii,ik) = v3s(n,ik)
              else
                h6(ii,ik) = zero
              endif
            end do
          end do
        endif

        mtmp(1:na,1:idx2(ib,ic)%size) = &
      &   matmul(h3(1:na,1:idx1(ib,ic,id)%size),h5(1:idx1(ib,ic,id)%size,1:idx2(ib,ic)%size))
        pbasis(1:na,1:nss) = matmul(mtmp(1:na,1:idx2(ib,ic)%size),h6(1:idx2(ib,ic)%size,1:nss))

        !debug! write(6,*) "=========================================================="

        do ii=1,nss
          write(u10) pbasis(:,ii)
          !debug! write(6,*) pbasis(:,ii)
        enddo

        !debug! write(6,*) "=========================================================="

        ! matrix of stage 3 vectors - asym case
        h6(:,:) = zero
        if(ib /= ic) then
          do ik=1,nsa
            do im=1,idx2asy(ib,ic)%size
              n = idx2asy(ib,ic)%start+im-1
              h6(im,ik) = -v3a(n,ik)* gmf2
            end do
          end do
        else ! ib==ic
          do ik=1,nsa
            im=0
            do ii=1,idx2(ib,ic)%size
              if(.NOT.s2sym(ib,kr,ii)) then
                im=im+1
                n = idx2asy(ib,ic)%start+im-1
                h6(ii,ik) = v3a(n,ik)
              else
                h6(ii,ik) = zero
              endif
            end do
          end do
        endif

        pbasis(1:na,1:nsa) = matmul(mtmp(1:na,1:idx2(ib,ic)%size),h6(1:idx2(ib,ic)%size,1:nsa))

        !debug! write(6,*) "=========================================================="

        do ii=1,nsa
          write(u09) pbasis(:,ii)
          !debug! write(6,*) pbasis(:,ii)
        enddo

        !debug! write(6,*) "=========================================================="

      end do
      ! end do; end do; end do

else ! HHSYM=.FALSE. ================================================

      write(u6,*) ' istage2=', istage2
      read(u13) nstates,itmp
      ! itmp above is effectively istage2 and must be equal to it

      ! allocate storage for eigen vectors and energies
      allocate(v3n(istage2,nstates),e3n(nstates))

      !debug! allocate(fullv(na,nn(1),nn(2),nn(3),nstates))


      do ii=1,nstates
        read(u13) e3n(ii),v3n(:,ii)
      end do

      ! if(nstates > 0)
      open(u10,file=fname//'all.dat', &
        access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)

      write(u10) nstates
      write(u10) e3n(:)

      ! find max size of the matrices which are needed for transformation to truncated basis
      idx1size_max = 0
      do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
        if(idx1size_max < idx1(ib,ic,id)%size) idx1size_max = idx1(ib,ic,id)%size
      enddo; enddo; enddo
      allocate(h3(na,idx1size_max))

      idx2size_max = 0
      do ic=1,nn(2); do ib=1,nn(1)
        if(idx2size_max < idx2(ib,ic)%size) idx2size_max = idx2(ib,ic)%size
      enddo; enddo
      !!allocate(h5(idx1size_max,idx2size_max))
      allocate(h6(idx2size_max,nstates),mtmp(na,idx2size_max),pbasis(na,nstates))

      ! do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
      !  replace the nested loops above by one loop
      do i=1,nn(3)*nn(2)*nn(1)
      
        !  unroll true indices
        id = int(int((i-1)/nn(1))/nn(2))+1
        ic = int((i-1)/nn(1))-nn(2)*(id-1)+1
        ib = i-nn(1)*(ic-1)-nn(1)*nn(2)*(id-1)

        ! h3 => vectrs1(ib,ic,id)%mat(:,:)
        do ia =1,idx1(ib,ic,id)%size
          do ni=1,na
            h3(ni,ia) =  vectrs1(ib,ic,id)%mat(ni,ia)
          end do
        end do

        !wrong! h5 => vectrs2(ib,ic)%mat(:,:)
        h5 => vectrs2(ib,ic)%mat(idx1(ib,ic,id)%start:idx1(ib,ic,id)%end,:)
        ! the pointer and the loop (below) are equivalent
        ! note: h5 needs to be allocated first
        !!do ik=1,idx2(ib,ic)%size
        !!  do im=1,idx1(ib,ic,id)%size
        !!    n = idx1(ib,ic,id)%start+im-1
        !!    h5(im,ik) = vectrs2(ib,ic)%mat(n,ik)
        !!  end do
        !!end do

      
        ! h6 => v3n(idx2(ib,ic)%start:idx2(ib,ic)%end,:)
        do ik=1,nstates
          do im=1,idx2(ib,ic)%size
            n = idx2(ib,ic)%start+im-1
            h6(im,ik) = v3n(n,ik)
          end do
        end do

        ! compute primitive basis (aggregation on the radial grid ib,ic,id)
        ! note that the order of the product is important and h3 x (h5 x h6) is wrong
        !  
        !  (h3 x h5) x h6 

        mtmp(1:na,1:idx2(ib,ic)%size) = &
      &   matmul(h3(1:na,1:idx1(ib,ic,id)%size),h5(1:idx1(ib,ic,id)%size,1:idx2(ib,ic)%size))

        pbasis(1:na,1:nstates) = &
      &   matmul(mtmp(1:na,1:idx2(ib,ic)%size),h6(1:idx2(ib,ic)%size,1:nstates))

        do ii=1,nstates
          write(u10) pbasis(:,ii)
          !debug! check the norm of the resulting w-f
          !fullv(:,ib,ic,id,ii)=pbasis(:,ii)
        enddo
      
      end do
      ! end do; end do; end do

      !debug! check the norm of the resulting w-f
      !do ii=1,nstates
      !  tmp = 0.0D0
      !  do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1); do ia=1,na
      !    tmp = tmp + fullv(ia,ib,ic,id,ii)**2
      !  end do; end do; end do; end do
      !  write(6,*) ii, tmp
      !end do


end if

      
close(u10)
if(HHSYM) close(u09)

if(HHSYM) deallocate(s1sym,s2sym)
deallocate(angsym)

! deallocate vectors
deallocate(indx,krn)

call time_stamp(6)

CONTAINS

FUNCTION i2c(i)
  USE types
IMPLICIT NONE
INTEGER(I4B) :: i
CHARACTER*1 i2c

if(i<0 .OR. i>9) then
 stop ' integer is out of range'
else
  i2c=char(48+i)
end if
END FUNCTION

END PROGRAM


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE setup_idx1(kr,istage1,ist1_max)
  USE types
  USE param
  USE workarrays, ONLY: idx1,s1sym,vectrs1,angsym
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: kr
INTEGER(I4B), INTENT(OUT) :: istage1,ist1_max
INTEGER(I4B) :: n,i1,i2,i3,icounter,ia
REAL(DP) :: temp

    if(HHSYM) then
      do i1=1,nn(1)
        i2=i1
          do i3=1,nn(3)
            !  find out the symmetry of i1=i2 states
            !  the algorithm works only if the eigenfuctions are
            !  either symmetric or asymmetric
            !
              do ia=1,idx1(i1,i2,i3)%size
                temp = zero; s1sym(i1,i3,kr,ia)= .TRUE.
                do n=1,na
                  temp = temp+ vectrs1(i1,i1,i3)%mat(n,ia)  &
                              *vectrs1(i1,i1,i3)%mat(angsym(kr,n),ia)
                end do
                if(temp < zero) s1sym(i1,i3,kr,ia)= .FALSE.
              end do
          end do
      end do
    end if

    ! -----------------------------------------------------------

    ! given i1 and i2, find the running index start and end positions
    ! as a function of i3 based on the length of selected
    ! angular basis functions

    if(test_flag >1) write(u6,'(" ib ic id    red. basis   start indx  end indx")')

      do i2=1,nn(2)
        do i1=1,nn(1)
          n=0

          do i3=1,nn(3)
            if(HHSYM .and. i1 > i2) idx1(i1,i2,i3)%size = idx1(i2,i1,i3)%size
            if(idx1(i1,i2,i3)%size==0) then
              idx1(i1,i2,i3)%start=n
              idx1(i1,i2,i3)%end=n
            else
              idx1(i1,i2,i3)%start=n+1
              n=n+idx1(i1,i2,i3)%size
              idx1(i1,i2,i3)%end=n
            end if

            if(test_flag >1) write(u6,'(3i3,3i12)') i1,i2,i3, &
              idx1(i1,i2,i3)%size,idx1(i1,i2,i3)%start,idx1(i1,i2,i3)%end
          end do
        end do
      end do

      ! count the total number of stage1 vectors and
      ! base vectors if additional symmetry is present
      icounter = 0; n=0
      do i2=1,nn(2); do i1=1,nn(1)
        icounter = icounter+idx1(i1,i2,nn(3))%end
        if(i1 <= i2)   n = n+idx1(i1,i2,nn(3))%end
      end do; end do
      istage1 = icounter
      write(u6,'(" stage 1 total:",i8," records")') istage1
        if(HHSYM) then
        istage1 = n
        write(u6,'(" of which basic",i8," records")') istage1
      end if

      !  find max number of selected levels per grid
      ist1_max=0
      do i3=1,nn(3); do i2=1,nn(2); do i1=1,nn(1)
          ist1_max= max(ist1_max,idx1(i1,i2,i3)%size)
      end do; end do; end do

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE setup_idx2(kr,istage2,ist2_max)
  USE types
  USE param
  USE workarrays, ONLY: idx1,idx2,idx2sym,idx2asy,s1sym,s2sym,vectrs2
IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: kr
INTEGER(I4B), INTENT(OUT) :: istage2,ist2_max
INTEGER(I4B) :: n,n1,ib,ic,id,icounter, nsym,nasy, ij,ii
REAL(DP) :: tmp0, tmp1

      !  find out the symmetry of ib=ic states
      !  the algorithm works only if the eigenfuctions are
      !  either symmetric or asymmetric
      !
      if(HHSYM) then
        do ib =1,nn(1)
          ic=ib
            do n1=1,idx2(ib,ic)%size
              tmp0 = zero; tmp1 = zero; s2sym(ib,kr,n1)= .TRUE.
              do id =1,nn(3)
                do ii =1,idx1(ib,ic,id)%size
                  n = idx1(ib,ic,id)%start+ii-1
                  if(s1sym(ib,id,kr,ii))then
                    tmp0 = tmp0+ vectrs2(ib,ib)%mat(n,n1)**2
                  else
                    tmp1 = tmp1+ vectrs2(ib,ib)%mat(n,n1)**2
                  endif
                end do
              end do
              if(tmp0 < tmp1) s2sym(ib,kr,n1)= .FALSE.
              !debug! write(6,*) 's2sym:',ib,kr,n1,  s2sym(ib,kr,n1)
            end do

        end do
      end if

    ! -----------------------------------------------------------

    if(test_flag >1) write(u6,'(" ib ic  full basis  red. basis   start indx  end indx")')
    n=0
    do ic=1,nn(2)
      do ib=1,nn(1)

        if(HHSYM .and. ib > ic) idx2(ib,ic)%size = idx2(ic,ib)%size
        if(idx2(ib,ic)%size==0) then
           idx2(ib,ic)%start=n
           idx2(ib,ic)%end=n
        else
           idx2(ib,ic)%start=n+1
           n=n+idx2(ib,ic)%size
           idx2(ib,ic)%end=n
        end if

        if(test_flag >1) write(u6,'(2i3,4i12)') ib,ic, idx1(ib,ic,nn(3))%end, &
          idx2(ib,ic)%size,idx2(ib,ic)%start,idx2(ib,ic)%end
      end do
    end do

    ! ---------------------------------------------------
    if(HHSYM) then
      nsym = 0  !  count symmetric functions
      nasy = 0  !  count asymmetric functions
      do ic =1,nn(2)
        do ib =1,ic
          idx2sym(ib,ic)%size=0
          idx2asy(ib,ic)%size=0
          if(ib == ic) then
            do ij=1,idx2(ic,ic)%size
              if(s2sym(ic,kr,ij)) then
                idx2sym(ic,ic)%size=idx2sym(ic,ic)%size+1
              else
                idx2asy(ic,ic)%size=idx2asy(ic,ic)%size+1
              end if
            end do
          else  ! ib /= ic
            idx2sym(ib,ic)%size= idx2(ib,ic)%size
            idx2asy(ib,ic)%size= idx2(ib,ic)%size
          end if

          if(idx2sym(ib,ic)%size==0) then
             idx2sym(ib,ic)%start=nsym
             idx2sym(ib,ic)%end=nsym
          else
             idx2sym(ib,ic)%start=nsym+1
             nsym=nsym+idx2sym(ib,ic)%size
             idx2sym(ib,ic)%end=nsym
          end if
              
          if(idx2asy(ib,ic)%size==0) then
             idx2asy(ib,ic)%start=nasy
             idx2asy(ib,ic)%end=nasy
          else
             idx2asy(ib,ic)%start=nasy+1
             nasy=nasy+idx2asy(ib,ic)%size
             idx2asy(ib,ic)%end=nasy
          end if

          if(ib /= ic) then
            idx2sym(ic,ib) = idx2sym(ib,ic)
            idx2asy(ic,ib) = idx2asy(ib,ic)
          end if

        end do
      end do
    end if
    ! ---------------------------------------------------

    ! the number of stage2 vectors
    ! istage2 = icounter
    istage2 = idx2(nn(1),nn(2))%end
    write(u6,'(" stage 2 total:",i8," records")') istage2
    if(HHSYM) then
      icounter = 0
      do ic=1,nn(2); do ib=1,nn(1)
        if(ib > ic) cycle
        icounter = icounter+idx2(ib,ic)%size
      end do; end do
      write(u6,'(" stage 2 basic:",i8," records")') icounter
      write(u6,'("           sym:",i8)') idx2sym(nn(1),nn(2))%end
      write(u6,'("           asy:",i8)') idx2asy(nn(1),nn(2))%end
    end if

    !  find max number of selected levels per grid
    ist2_max=0
    do ic=1,nn(2); do ib=1,nn(1)
        ist2_max= max(ist2_max,idx2(ib,ic)%size)
    end do; end do

END SUBROUTINE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
