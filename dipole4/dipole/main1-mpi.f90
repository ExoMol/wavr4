! 01/07/2008 INK
! dipole moment calculations
! NOTE: this program computes only integrals
! <f|mux|i>, <f|muy|i>, <f|muz|i>
!
! the code is meant to be used for band intecities calculations
! therefore J > 0 is not applicable

! convention: final state x1 <- initial state x2 or x
!
! beware of the possible confusion between Kr and Jr
! they should be the same since the code works so far for J=K=0 only
! but once we expand to J > 0, care must be taken to resolve everything properly

PROGRAM dipole
  USE types
  USE param
  USE workarrays
  USE threej
  USE mpi
IMPLICIT NONE

!  parameters
!
INTEGER(I4B), PARAMETER :: u11=11, u12=12
REAL(DP), PARAMETER :: TINY= 1.0D-14

REAL(DP), ALLOCATABLE :: dt(:,:),dm(:,:)
REAL(DP), ALLOCATABLE :: vm1(:,:),vm2(:,:)
REAL(DP), ALLOCATABLE :: e1(:), e2(:), mtmp(:,:)

REAL(DP), ALLOCATABLE :: theta1(:),ttheta1(:,:,:)
REAL(DP), ALLOCATABLE :: theta(:), ttheta(:,:,:)
!REAL(DP), ALLOCATABLE :: phi(:),   tphi(:,:)
REAL(DP), ALLOCATABLE :: phi(:),tphi_cos(:,:),tphi_sin(:,:)
REAL(DP), ALLOCATABLE :: q1(:),q2(:),q3(:)
REAL(DP), ALLOCATABLE :: tp(:),tp1(:),tp2(:),tpx(:),tpy(:)

REAL(DP), ALLOCATABLE :: mu_grid(:,:,:,:,:,:)
!!REAL(DP), ALLOCATABLE :: mux_grid(:,:,:,:,:,:)
!!REAL(DP), ALLOCATABLE :: muy_grid(:,:,:,:,:,:)
!!REAL(DP), ALLOCATABLE :: muz_grid(:,:,:,:,:,:)

REAL(DP) :: zpe, tmp
REAL(DP) :: one_over_sqrt2, alpha

INTEGER(I4B), ALLOCATABLE :: indx1(:,:),indx2(:,:)
TYPE(triple), ALLOCATABLE :: krn1(:),krn2(:)

INTEGER(I4B) :: istage1,istage2,istage3
INTEGER(I4B) :: ist1_max,ist2_max
!
INTEGER(I4B) :: n, ii,i1,i2,i3,  ii1, ii2
INTEGER(I4B) :: ia,ib,ic,id,im, m, ia1
INTEGER(I4B) :: i, icounter, idiag, nsize, icut
INTEGER(I4B) :: j1, l1, k1, m1, kp1, jr1, kr1
INTEGER(I4B) :: j2, l2, k2, m2, kp2, jr2, kr2
INTEGER(I4B) :: dkr, kr1_start, kr1_end, dkdif, dksum, dkr1
INTEGER(I4B) :: nlev,nnmax,istage2s,istage2a,itmp

INTEGER(I4B) :: jp1, j_parity1, l_parity1, jl_parity1, p1
INTEGER(I4B) :: jp2, j_parity2, l_parity2, jl_parity2, p2
INTEGER(I4B) :: nstates1, na1
INTEGER(I4B) :: nstates2, na2
!INTEGER(I4B) :: nstates     ! number of lowest states to be used
!
!TYPE(triple), ALLOCATABLE :: idx1(:,:,:),idx2(:)

CHARACTER(12) :: fname1, fname2
CHARACTER(3) :: smt
CHARACTER(1) :: dipole_projection
CHARACTER(LEN=80) title

! MPI relevant vars
INTEGER :: num_procs, proc_id, ierr, buffer_size
INTEGER, ALLOCATABLE :: proc2grid_map(:)
REAL(DP), ALLOCATABLE :: dt_buffer(:,:)

LOGICAL :: skip

! ====================================================================

call MPI_INIT(ierr)
call MPI_COMM_SIZE( MPI_COMM_WORLD, num_procs, ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, proc_id, ierr )

call MPI_BARRIER( MPI_COMM_WORLD, ierr )
!debug! write(u6,*) " im proc num", proc_id," of ", num_procs

! ====================================================================

one_over_sqrt2 = one/sqrt(two)

!call time_stamp(6)

! redirect output to a file
!open(u6,file='out.txt')

! ====================================================================
! INPUT - WAVR4 part
!

if(proc_id == 0) write(u6,'(" INPUT: reading WAVR4 control parameters",/)')

open(u5,file='input.txt',status='old')

call input_data()

close(u5)

! just set these for backward compatibility with WAVR4
qe(1)=re(1); qe(2)=re(2); qe(3)=re(3); qe(4)=zero; qe(5)=zero; qe(6)=zero

allocate(proc2grid_map(nn(3)*nn(2)*nn(1)))

! make sure there's no printout from all the nodes
if(proc_id /= 0) then
  test_flag=0
endif

! ====================================================================
! INPUT - dipole part
!

if(proc_id == 0) write(u6,'(" INPUT: reading DIPOLE4 control parameters",/)')

open(u5,file='inp.dipole.txt',status='old')

if(proc_id == 0) write(u6,'(80("="))')

! final state
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) jp1, j_parity1, l_parity1, jl_parity1, jr1, smt
if(proc_id == 0) write(u6,'(1x,5i5,A8)') jp1, j_parity1, l_parity1, jl_parity1, jr1, smt

fname1 = 'x'//i2c(jp1)//i2c(j_parity1)//i2c(l_parity1)//i2c(jl_parity1)//'-J'//i2c(jr1)//'-'//smt

if( HHSYM .AND. smt == "all") then
  if(proc_id == 0) write(6,*) ' INPUT error: incompatible file name, check HHSYM'
  stop
endif

! initial state
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) jp2, j_parity2, l_parity2, jl_parity2, jr2, smt
if(proc_id == 0) write(u6,'(1x,5i5,A8)') jp2, j_parity2, l_parity2, jl_parity2, jr2, smt

fname2 = 'x'//i2c(jp2)//i2c(j_parity2)//i2c(l_parity2)//i2c(jl_parity2)//'-J'//i2c(jr2)//'-'//smt

if( HHSYM .AND. smt == "all") then
  if(proc_id == 0) write(6,*) ' INPUT error: incompatible file name, check HHSYM'
  stop
endif

! dipole projection
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) dipole_projection
if(proc_id == 0) write(u6,'(1x,A8)') dipole_projection

if(proc_id == 0) write(u6,'(80("="))')

close(u5)

if(proc_id == 0) write(u6,'(/," INPUT: done.",/)')

! ====================================================================
! CHECK INPUT parameters

if(abs(jr1-jr2) > 1 ) then
  if(proc_id == 0) write(u6, '("  INPUT: delta J cannot be greater than 1")')
  stop
endif

if(j_parity1 > j_parity_max .OR. j_parity2 > j_parity_max) then
  if(proc_id == 0) write(u6, '("  INPUT: j_parity is higher than",i3)') j_parity_max
  stop
end if

if(l_parity1 > l_parity_max .OR. l_parity2 > l_parity_max) then
  if(proc_id == 0) write(u6, '("  INPUT: l_parity is higher than",i3)') l_parity_max
  stop
end if

if(jl_parity1 > jl_parity_max .OR. jl_parity2 > jl_parity_max) then
  if(proc_id == 0) write(u6, '("  INPUT: jl_parity is higher than",i3)') jl_parity_max
  stop
end if

! define missing (derived) quantum numbers
p1=mod((jp1+jr1),2)
p2=mod((jp2+jr2),2)

!debug! for norm testing comment out the if statement
if(p1==p2) then
  if(proc_id == 0) write(u6, '("  INPUT: transitions between the states with the same parity are forbidden")')
  stop
end if

if(proc_id == 0) write(u6,*) " external files to be opened..."
if(proc_id == 0) write(u6,*) "   w-f of   final states:  ", fname1
if(proc_id == 0) write(u6,*) "   w-f of initial states:  ", fname2

! ====================================================================

!  find max number of angular states (namax) and mmax
!
call angular_states_max()

! final state angular basis indexing
allocate(indx1(namax,6),krn1(0:krmax))
call angular_states_qnumbers(indx1,krn1, jr1,jp1,j_parity1,l_parity1,jl_parity1)

! initial state angular basis indexing
allocate(indx2(namax,6),krn2(0:krmax))
call angular_states_qnumbers(indx2,krn2, jr2,jp2,j_parity2,l_parity2,jl_parity2)

! radial grid points

allocate(t1(nn(1),nn(1)),t2(nn(2),nn(2)),t3(nn(3),nn(3)), &
         oner1(nn(1),nn(1)),oner2(nn(2),nn(2)),oner3(nn(3),nn(3)), &
         q1(nn(1)),q3(nn(3)),q2(nn(2)))

call radial_grids(q1,q2,q3,t1,t2,t3,oner1,oner2,oner3)

deallocate(t1,t2,t3,oner1,oner2,oner3)

allocate(theta1(nt1),ttheta1(0:jmax,0:mmax,nt1), &
         theta(nt),ttheta(0:lmax,0:kmax,nt),  &
         phi(nphi),tphi_cos(0:2*kmax,nphi),tphi_sin(0:2*kmax,nphi))

call grid_theta(theta1,ttheta1,jmax,mmax,nt1)
call grid_theta(theta,ttheta,lmax,kmax,nt)
! this is a modified version of grid_phi found in WAVR4
! in addition to cos grid it returns sin grid as well
call grid_phi(phi,tphi_cos,kmax,nphi)
call grid_phi_sin(phi,tphi_sin,kmax,nphi)
! allocate temp arrays to keep products of wf grids
allocate(tp(size(tphi_cos,2)),tp1(size(ttheta1,3)),tp2(size(ttheta,3)))

allocate(mu_grid(nt1,nt,nphi,0:1,0:1,0:1))

! ===================================================================

! ===================================================================

! open file with final eigen functions
open(u11,file=fname1//'.dat', &
  status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
if(idiag /=0) then 
  if(proc_id == 0) write(6,*) '  cant open ',fname1//'.dat';  stop
end if

! open file with initial eigen functions
open(u12,file=fname2//'.dat', &
  status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
if(idiag /=0) then 
  if(proc_id == 0) write(6,*) '  cant open ',fname2//'.dat';  stop
end if

! ===================================================================

if(proc_id == 0) call time_stamp(6)

read(u11) nstates1
read(u12) nstates2

! dipole transitions matrix
allocate(dt(nstates1,nstates2))
allocate(e1(nstates1),e2(nstates2))

dt(:,:)=ZERO

do i=1,nn(3)*nn(2)*nn(1)
  proc2grid_map(i) = mod(i,num_procs)
  !if(proc_id == 0) write(u6,*) " grid point", i, " mapped to ", proc2grid_map(i)
enddo

! Implement eqs. (32)-(50) of the Appendix.
! Note that X(J,J') of eq. (32) is computed here as a sum of 
!   X_x(J,J';K') which is a function of mu_x
!   X_y(J,J';K') which is a function of mu_y
!   X_z(J,J';K') which is a function of mu_z
! Therefore it is actually X_a(J,J') which is computed here, 
! "a" being specified as "dipole_projection" in the input. So
! X_a(J,J') = X_a(J,J';0)+X_a(J,J';1)+...
! The 3J symbol present in eqs. (33)-(50) is plugged into one of DGEMMs below
! as the Alpha coefficient  C = alpha*A*B+beta*C

dkr_loop: do dkr=-1,1,1

  ! special cases
  if( jr1==0 .AND. jr2==0 .AND. dkr /=0 ) cycle dkr_loop
  if( jr1==0 .AND. jr2==1 .AND. dkr ==-1 ) cycle dkr_loop
  if( jr1==1 .AND. jr2==0 .AND. dkr == 1 ) cycle dkr_loop

  if( dipole_projection == "x" .AND. abs(dkr) /=1 ) cycle dkr_loop
  if( dipole_projection == "y" .AND. abs(dkr) /=1 ) cycle dkr_loop
  if( dipole_projection == "z" .AND. abs(dkr) /=0 ) cycle dkr_loop

  !==============================================================
  ! rewind the data files to the start 
  ! (strictly speaking this is only necessary for mu_x and mu_y
  !  because for mu_z the data is read only once; therefore we have some
  !  redundant job for mu_z but this should not have high impact)
  rewind(u11)
  rewind(u12)
  !  read only to advance the record (we have read nstates already)
  read(u11) itmp
  read(u12) itmp
  !  read the energies (perhaps for the second time in case of mu_x and mu_y)
  read(u11) e1(:)
  read(u12) e2(:)

  !==============================================================

  kr1_start = 0
  if( dkr < 0 )   kr1_start = 1
  kr1_end = min(jr1, jr2-dkr)

  ! wind the data files forward if abs(dkr) = 1
  if( dkr == -1) then                   ! K'=1, K=0
    ! kr1 = 0
    if(proc_id == 0) write(6,*) "K1=0 is skipped"
    na1 = krn1(0)%size
    allocate(vm1(na1,nstates1))
    do i=1,nn(3)*nn(2)*nn(1)
        do ii=1,nstates1
          read(u11) vm1(:,ii)
          !debug!
          !!do ia=1,na1
          !!  write(6,*) "K1=0 skip:", ii, ia, vm1(ia,ii)
          !!enddo
        enddo
    enddo
    deallocate(vm1)
   elseif( dkr == 1) then                ! K'=0, K=1
    ! kr2 = 0
    if(proc_id == 0) write(6,*) "K2=0 is skipped"
    na2 = krn2(0)%size
    allocate(vm2(na2,nstates2))
    do i=1,nn(3)*nn(2)*nn(1)
        do ii=1,nstates2
          read(u12) vm2(:,ii)
          !debug!
          !!do ia=1,na2
          !!  write(6,*) "K2=0 skip:", ii, ia, vm2(ia,ii)
          !!enddo
        enddo
    enddo
    deallocate(vm2)
  endif

  kr1_loop: do kr1= kr1_start, kr1_end

    na1 = krn1(kr1)%size
    kr2 = kr1+dkr
    na2 = krn2(kr2)%size

    if(proc_id == 0) write(6,'( "mu_",A1,": K1=",i2," K2=",i2," na1=",i6," na2=",i6)') dipole_projection,kr1,kr2,na1,na2

    ! compute the factor which includes 3J symbol
    if    ( dkr== 0) then
      alpha = my3j(jr2,kr1,jr1,-kr1)*TWO
    elseif( dkr== 1) then
      alpha = my3j(jr2,(-kr1-1),jr1,kr1)*SQRT2
    elseif( dkr==-1) then
      alpha = my3j(jr2,(kr1-1),jr1,-kr1)*SQRT2
    endif

    !debug! for norm testing 3J factor  to 2
    !!alpha = 2.D0
    if(proc_id == 0) write(6,*) " alpha=", alpha

    skip = .FALSE.
    if(abs(alpha) < TINY)  skip = .TRUE.

    allocate(vm1(na1,nstates1))
    allocate(mtmp(na1,nstates2))
  
    allocate(vm2(na2,nstates2))
    allocate(dm(na1,na2))
  
    ! do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
    !  replace the nested loops above by one loop
    grid: do i=1,nn(3)*nn(2)*nn(1)
  
      !  unroll true indices
      id = int(int((i-1)/nn(1))/nn(2))+1
      ic = int((i-1)/nn(1))-nn(2)*(id-1)+1
      ib = i-nn(1)*(ic-1)-nn(1)*nn(2)*(id-1)
  
      ! read final states for i1,i2,i3,K
      !!write(6,*) "K1=",kr1
      do ii=1,nstates1
         read(u11) vm1(1:na1,ii)
         !debug!
         !!do ia=1,na1
         !!  write(6,*) ii, ia, vm1(ia,ii)
         !!enddo
      enddo
  
      ! read initial states
      !!write(6,*) "K2=",kr2
      do ii=1,nstates2
          read(u12) vm2(1:na2,ii)
         !debug!
         !!do ia=1,na2
         !!  write(6,*) ii, ia, vm2(ia,ii)
         !!enddo
      enddo

      if(skip) cycle grid
  
      ! all workers read but then skip if not mapped to grid point
      if(proc_id /= proc2grid_map(i)) cycle grid

      call dipolematrix(q1(ib),q2(ic),q3(id))
      ! passing arrays in context: theta1,ttheta1, theta,ttheta, phi,tphi_cos,tphi_sin 

      !debug!
      !!do ii=1,na1; do ia=1,na2
      !!  if(abs(dm(ia,ii)) > 1.D-10) write(6,*) ia,ii,dm(ia,ii)
      !!enddo; enddo
  
      ! vm1(nstates1,na1)^T x dm(na1,na2) x vm2(na2,nstates2)

      !write(6,*) dm(1,1)  

      !  gemm ([],[],m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
      !  performs matrix-matrix multiplication of the form 
      !  C = alpha [A] [B] + beta [C]
      !  
      !  [] is 'N', 'T' or 'C'
      !  [A] m.k ; A(lda,:)
      !  [B] k.n ; B(ldb,:)
      !  [C] m.n ; C(ldc,:)

      ! dm x vm2 = mtmp
      call dgemm('N','N',na1,nstates2,na2,1.D0,dm,na1,vm2,na2,0.D0,mtmp,na1)
      ! alpha* vm1 x mtmp + dt = dt  (accumulation over radial grid)
      call dgemm('T','N',nstates1,nstates2,na1,alpha,vm1,na1,mtmp,na1,1.D0,dt,nstates1)

      !debug! for norm testing uncomment the loops below and comment out the 2 dgemms above
      !call dgemm('T','N',nstates1,nstates2,na1,1.D0,vm1,nstates1,vm2,na1,1.D0,dt,nstates1)
      !do i1=1,nstates1
      !  do i2=1,nstates2
      !    do i3=1,na1
      !      dt(i1,i2)=dt(i1,i2)+vm1(i3,i1)*vm2(i3,i2)
      !    end do
      !  end do
      !end do
  
    end do grid
  
    deallocate(dm)
    deallocate(vm2)
    deallocate(vm1,mtmp)

    !debug!  dt(:,:)=ZERO

  end do kr1_loop

end do dkr_loop

! deallocate all unnecessary arrays

deallocate(indx1,krn1)
deallocate(indx2,krn2)

! gather dt matrices from all workers
allocate(dt_buffer(nstates1,nstates2))

!debug! write(u6,*) " proc ", proc_id, " done"
call MPI_BARRIER( MPI_COMM_WORLD, ierr )
if(proc_id == 0) write(u6,*) " job's done; about to reduce"
buffer_size = nstates1*nstates2
call MPI_REDUCE(dt, dt_buffer, buffer_size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

! TODO: create line frequency list and sort it

! S = (2J+1)(2J'+1) (X_x + X_y + X_z)^2  where dt = X

tmp = 0.0D0
do ii1=1,nstates1
  do ii2=1,nstates2

    !!write(6,*) ii1,ii2, dt_buffer(ii1,ii2)
    ! find the largest transition
    if (abs(dt_buffer(ii1,ii2)) > tmp) tmp = abs(dt_buffer(ii1,ii2))

    ! set threshold for non-zero transitions
    !if (abs(dt_buffer(ii1,ii2)) < 2.D-14) dt_buffer(ii1,ii2) = 0.0D0

    !debug!    if(proc_id == 0 .AND. abs(dt_buffer(ii1,ii2)) > 1.D-12) 
    if(proc_id == 0) &
      write(u6,'(i4,"<-",i4,f12.4,f12.2,es16.4)') ii1,ii2, e1(ii1)-e2(ii2),e2(ii2),dt_buffer(ii1,ii2)

  enddo
enddo

!!if(proc_id == 0) write(u6,'(" norm=",es16.4)') tmp
if(proc_id == 0) write(u6,'(" largest X=",es16.4)') tmp

if(proc_id == 0) call time_stamp(6)

close(u11)
close(u12)

call MPI_FINALIZE(ierr)


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


SUBROUTINE masses()
  USE types
  USE param
IMPLICIT NONE
REAL(DP) :: mu1, mu2, mu3, mu4
! INTEGER(I4B), INTENT(IN) :: opt

!  assuming:     q3 -> R,       q1 -> d1,       q2 -> d2
!  resp.      mu_q3 -> mu_R, mu_q1 -> mu_d1, mu_q2 -> mu_d2

msum = mass1+mass2+mass3+mass4
!  instead of multiplying back by msum we will work with dimensional masses
mu1 = mass1  ! /msum
mu2 = mass2  ! /msum
mu3 = mass3  ! /msum
mu4 = mass4  ! /msum

if (abs(opt)==1) then
  !  Jacobi vectors
  if(proc_id == 0) write(u6,'(" opt 1: Jacobi vectors",/)')
  muq(1) = mu1*mu2/(mu1+mu2)
  muq(2) = mu3*(mu1+mu2)/(mu1+mu2+mu3)
  muq(3) = mu4*(mu1+mu2+mu3) /msum
else if (abs(opt)==2) then
  !  Radau vectors
  if(proc_id == 0) write(u6,'(" opt 2: Radau vectors",/)')
  muq(1) = mu1
  muq(2) = mu2
  muq(3) = mu3
else if (abs(opt)==3) then
  !  Diatom-diatom vectors
!   these are original mu's
  if(proc_id == 0) write(u6,'(" opt 3: Diatom-diatom vectors",/)')
  muq(1) = mu1*mu2/(mu1+mu2)
  muq(2) = mu3*mu4/(mu3+mu4)
  muq(3) = (mu1+mu2)*(mu3+mu4) /msum
!   these mu have been modified to be compatible with option 4
!  write(u6,'(" opt 3: Diatom-diatom vectors, opt4 compatible",/)')
!  muq(1) = mu1*mu3/(mu1+mu3)
!  muq(2) = mu2*mu4/(mu2+mu4)
!  muq(3) = (mu1+mu3)*(mu2+mu4) /msum
else if (abs(opt)==4) then
  !  Orthogonalized satellite vectors
  if(proc_id == 0) write(u6,'(" opt 4: Orthogonal satellite vectors",/)')
  muq(1) = mu1
  muq(2) = mu2
  muq(3) = mu3*mu4/(mu3+mu4)
else
  if(proc_id == 0) write(u6,*) ' masses: no such option:', opt
  stop
end if

mud1 = muq(1)
mud2 = muq(2)
mur  = muq(3)

END SUBROUTINE


! INK 16/06/2007
! this is a trimmed down copy of input_data() in rv.sym (or WAVR4 proper)

SUBROUTINE input_data()
  USE types
  USE param
IMPLICIT NONE
REAL(DP) :: temp
INTEGER(I4B) :: i, iauto
CHARACTER(LEN=80) title

if(proc_id == 0) write(u6,'(" Four Atomic RoVibrational Program",/, &
           & " INPUT: reading input...",/, &
           & 30("="),"< file:  input.txt >",30("="))')
! ======================================================================
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
!  masses:
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) mass1,mass2,mass3,mass4
if(proc_id == 0) write(u6,'(4(2x,f13.9))') mass1,mass2,mass3,mass4
!  coordinate system
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) opt
if(proc_id == 0) write(u6,'(i5)') opt
!  parameters of stretching basis functions  and the grid size
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
do i=1,3
  read(u5,*) igq(i),re(i), we(i), De(i), nn(i)
  if(proc_id == 0) write(u6,'(i4,f8.4,f12.3,f12.3,i6)') igq(i),re(i),we(i),De(i),nn(i)
end do
!  flags:
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) angular_problem_only_flag, optimize_flag, test_flag, expansion_flag
if(proc_id == 0) write(u6,'(11x,i1,23x,i1,13x,i1,13x,i1)') &
  angular_problem_only_flag, optimize_flag, test_flag, expansion_flag
!  angular basis functions
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) jmax,lmax,kmax, jrmax, krmax, &
  j_parity_max, l_parity_max, jl_parity_max
if(proc_id == 0) write(u6,'(1x,8i5)') jmax,lmax,kmax, jrmax, krmax, &
  j_parity_max, l_parity_max, jl_parity_max
!  expansion size
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) ne1, ne2, ne3
if(proc_id == 0) write(u6,'(1x,3i5)') ne1, ne2, ne3
!  angular grid points size
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) nt1, nt, nphi, iauto
if(proc_id == 0) write(u6,'(1x,4i5)') nt1, nt, nphi, iauto
!  zero energy:
read(u5,'(A80)') title; if(proc_id == 0) write(u6,'(A80)') title
read(u5,*) enzero
if(proc_id == 0) write(u6,'(1x,f20.10)')  enzero

if(proc_id == 0) write(u6,'(80("="))')
! ======================================================================
if(proc_id == 0) write(u6,'(" INPUT: processing input...",/)')

call masses()
!   take Radau formulae for masses to test 3D angular + rotation block
!   because masses 1-3 will be the same as in the input and an analytical
!   solution is possible. Also set V=0
if(proc_id == 0) write(u6,'(" adapted masses:",/,3f14.8,/)')  muq

! basis size quantum numbers:
!
! this is just a quick fix:
! mmax will be defined in main -> angular_basis_size
mmax = -99

if(kmax>lmax) then
  stop ' input: kmax > lmax'
end if

if(krmax>jrmax) then
  stop ' input: krmax > jrmax'
end if

! check j, l, jl parity numbers

if(proc_id == 0) write(u6, '("  j-, l-, (j+l)-parity:")')
if(j_parity_max==0) then
  if(proc_id == 0) write(u6, '("  NO:  q1 -> -q1 symmetry (j-parity)")')
else if(j_parity_max==1) then
  if(proc_id == 0) write(u6, '("  YES: q1 -> -q1 symmetry (j-parity)")')
else
  if(proc_id == 0) write(u6, '("  INPUT: WRONG j_parity_max",i3)') j_parity_max
  stop
end if

if(l_parity_max==0) then
  if(proc_id == 0) write(u6, '("  NO:  q2 -> -q2 symmetry (l-parity)")')
else if(l_parity_max==1) then
  if(proc_id == 0) write(u6, '("  YES: q2 -> -q2 symmetry (l-parity)")')
else
  if(proc_id == 0) write(u6, '("  INPUT: WRONG l_parity_max=",i3)') l_parity_max
  stop
end if

if(jl_parity_max==0) then
  if(proc_id == 0) write(u6, '("  NO:  q3 -> -q3 symmetry (jl-parity)")')
else if(jl_parity_max==1) then
  if(proc_id == 0) write(u6, '("  YES: q3 -> -q3 symmetry (jl-parity)")')
else
  if(proc_id == 0) write(u6, '("  INPUT: WRONG jl_parity_max=",i3)') jl_parity_max
  stop
end if

!  can't have all j_parity_max=l_parity_max=jl_parity_max=1
if(j_parity_max==1 .AND. l_parity_max==1 .AND. jl_parity_max==1) then
  if(proc_id == 0) write(u6, '("  INPUT: when all parities = 1, jl_parity is redundant")')
  if(proc_id == 0) write(u6, '("  INPUT: jl_parity_max is reset to 0")')
  jl_parity_max=0
end if

if(proc_id == 0) write(u6, '(/," angular basis size:              ", &
  & " jmax=",i3," lmax=",i3," kmax=",i3,/)') jmax, lmax, kmax

!  number of grid points for angular integration
!  matrix integrals (expansion_flag=0,1)
!    see grid3d.f90 for details
!    full grid must be: nt1 >= 2*(jmax+1); nt >= 2*(lmax+1); nphi >= kmax+1
!   we use sym reduced: nt1 >= (jmax+1); nt >= (lmax+1); nphi >= (kmax+1+1)/2
!   the latter is because we want 2*nphi >= kmax+1
!  IF iauto =/= 1 then expansion grid is taken from the input

if (iauto == 1) then
  !  NOTE: this is a bare minimum
  if(proc_id == 0) write(u6,'(" (auto) signle ang grid, iauto=",i2,":")') iauto
  if(expansion_flag == 0) then
    nt1 = jmax+1
    nt  = lmax+1
    nphi = (kmax+1+1)/2
  else ! if(expansion_flag == 1) then
    nt1 = (ne1+1)
    nt  = (ne2+1)
    nphi = (ne3+1+1)/2
  end if
else if (iauto == 2) then
  if(proc_id == 0) write(u6,'(" (auto) double ang grid, iauto=",i2,":")') iauto
  if(expansion_flag == 0) then
    nt1 = 2*jmax+1
    nt  = 2*lmax+1
    nphi = kmax+1
  else ! if(expansion_flag == 1) then
    nt1 = 2*ne1+1
    nt  = 2*ne2+1
    nphi = ne3+1
  end if
else if (iauto == 3) then
  if(proc_id == 0) write(u6,'(" (auto) 4x     ang grid, iauto=",i2,":")') iauto
  if(expansion_flag == 0) then
    nt1 = 4*jmax+1
    nt  = 4*lmax+1
    nphi = 2*kmax+1
  else ! if(expansion_flag == 1) then
    nt1 = 4*ne1+1
    nt  = 4*ne2+1
    nphi = 2*ne3+1
  end if
else
  if(proc_id == 0) write(u6,'(" (input) angular   grid, iauto=",i2,":")') iauto
end if
if(proc_id == 0) write(u6,'(" nt1 =",i3," nt  =",i3," nphi=",i3)') nt1,nt,nphi

nagrid = nt1*nt*nphi
if(proc_id == 0) write(u6,'(/," full angular grid: nagrid = nt1*nt*nphi =",i6,/)') nagrid

! energy cutoffs and zero energy

if(proc_id == 0) write(u6,'(" input zero energy   (enzero)=",f10.1)') enzero

if(proc_id == 0) write(u6,'(/," INPUT: done.",/)')

END SUBROUTINE



SUBROUTINE dipolematrix(r1,r2,r3)
  USE types
  USE param
IMPLICIT NONE
!  REAL(DP), INTENT(OUT) :: dm(na1,na2)
  REAL(DP), INTENT(IN) :: r1,r2,r3
!  REAL(DP), INTENT(IN) :: theta1(nt1),ttheta1(0:jmax,0:mmax,nt1)
!  REAL(DP), INTENT(IN) :: theta(nt),   ttheta(0:lmax,0:kmax,nt)
!  REAL(DP), INTENT(IN) :: phi(nphi),tphi_cos(0:2*kmax,nphi),tphi_sin(0:2*kmax,nphi)
REAL(DP) ::  tmp1, integral, my3j
REAL(DP) ::  norm1, norm2, fct, mux, muy, muz
INTEGER(I4B) :: i1,i2,i3, dkp,dkm
! INTEGER(I4B) :: n, j, l, k, kp, kr, m
! INTEGER(I4B) :: kr1,kr2 ! inherited from the calling context
INTEGER(I4B) :: n1,j1,l1,k1,kp1,m1
INTEGER(I4B) :: n2,j2,l2,k2,kp2,m2
INTEGER(I4B) :: i, ii, ii1, ia1, ia2
INTEGER(I4B) :: is1, is2, is3
LOGICAL :: dkdif_negative, dksum_negative

! ia1,kr1, ia2,kr2, dkr  etc must be inherited from the calling routine

! generate symmetrized (sym & anti-sym grids for mux, muy, muz
!
! SELECT the right dipole projection

if(abs(dkr)==1) then
  if     (dipole_projection == 'x')  then
               call dipole_grid_x(r1,r2,r3)
  elseif (dipole_projection == 'y')  then
               call dipole_grid_y(r1,r2,r3)
  else
    ! this should not happen but if we get here make sure we return zero
    do ia1=1,krn1(kr1)%size
      do ia2=1,krn2(kr2)%size
        dm(ia1,ia2) = ZERO
      end do
    end do
    return
  endif
elseif(dkr==0) then
  if     (dipole_projection == 'z')  then
               call dipole_grid_z(r1,r2,r3)
  else
    ! this should not happen but if we get here make sure we return zero
    do ia1=1,krn1(kr1)%size
      do ia2=1,krn2(kr2)%size
        dm(ia1,ia2) = ZERO
      end do
    end do
    return
  endif
else
  if(proc_id == 0) write(6,*) " dipolematrix: delta K is greater than 1"
  stop
endif


! start the integrals

do ia1=1,krn1(kr1)%size

  n1 =ia1+krn1(kr1)%start-1

  l1 =indx1(n1,1)
  j1 =indx1(n1,2)
  kp1=indx1(n1,3)
  !kr1=indx1(n1,4)
  k1 =indx1(n1,5)
  m1 =indx1(n1,6)

  if(kr1==0 .AND. k1==0) then
    norm1 = 0.5_dp
  else
    norm1 = one_over_sqrt2
  end if

  if (kp1 == -1 .AND. mod(k1,2)==1) norm1 = -norm1

  !debug! write(6,'(" ia1 ia2  j1  m1  l1  k1  j2  m2  l2  k2")')

  do ia2=1,krn2(kr2)%size
    n2 =ia2+krn2(kr2)%start-1

    l2 =indx2(n2,1)
    j2 =indx2(n2,2)
    kp2=indx2(n2,3)
    !kr2=indx2(n2,4)
    k2 =indx2(n2,5)
    m2 =indx2(n2,6)

    if(kr2==0 .AND. k2==0) then
      norm2 = 0.5_dp
    else
      norm2 = one_over_sqrt2
    end if

    if (kp2 == -1 .AND. mod(k2,2)==1) norm2 = -norm2

    ! it seems like we dont need the extra factor of 4.D0
    fct = norm1*norm2 

    dksum = kp1*k1+kp2*k2
    dkdif = kp1*k1-kp2*k2

    if (mod(kr1+((kp1+1)/2)*k1+((kp2+1)/2)*k2, 2) ==1) fct = -fct

! 0)
! We are computing here dipole function integrals according to the formulae in the Appendix
!       1/(2*Pi) * mu * [...] * N(K',k') * N(K,k) * (-1)^(...)
! the 3J symbol and the factor of 2 or sqrt(2) are applied in the calling subroutine
! 
! cos_phi and sin_phi grids are normalised on the full range 0..2Pi
!
! the factor of 4 above comes from the fact that the cos and sin grids are computed
! assuming 1/(2*Pi) norm whereas the integral we are computing is
! i.e. 4 times bigger
!
! 1)
! since E^* mux =  mux;  mux(-phi) =  mux(phi)
!       E^* muy = -muy;  muy(-phi) = -muy(phi)
!       E^* muz =  muz;  muz(-phi) =  muz(phi)
!
! 2)
! integration over phi is in the range from 0..Pi
! the result is then doubled because the integrand must be symmetric
! i.r.t phi=Pi or it is zero
!
! 3)
! the symmetry i.r.t. Pi/2 axis
! mu_grid has always symmetrized values f(+)=f(x)+f(-x) for is1, is2, is3=0
!               anti-symmetrized values f(-)=f(x)-f(-x) for is1, is2, is3=1
! therefore there is no need for doubling the integral value after halving the
! integration regions.
! f(+) contains symmetric portion of mux, muy and muz 
! f(-) contains anti-symmetric portion of mux, muy and muz
!
! symmetry of cos(n phi) & sin(n phi) in the respect to phi= Pi/2 axis:
!        n= 1  2  3
! cos(nx)   a  s  a ...
! sin(nx)   s  a  s ...
!
! therefore when computing integrals [mu_z cos(n phi)] or [mu_x con(n phi)]
! for even n we need only symmetric part of mu_z or mu_x       => is3 = 0
! for odd  n we need only anti-symmetric part of mu_z or mu_x  => is3 = 1
! it is opposite for integrals of [mu_y sin(n phi)]
! for even n we need only anti-symmetric part of mu_y          => is3 = 1
! for odd  n we need only symmetric part of mu_y               => is3 = 0
! hence is3 needs to be modified for mu_y

    is1=mod(j1+m1+j2+m2,2)
    is2=mod(l1+k1+l2+k2,2)
    is3=mod(abs(dksum),2)
    if (dipole_projection == 'y') is3=mod(abs(dksum+1),2)

    ! We need to make sure dksum and dkdif are not negative because 
    ! cos(k x_i) and sin(k x_i) are computed only for k >= 1.
    ! We can safely take abs() for cos(dksum) and cos(dkdif) because cos(-x)=cos(x)
    ! However we need to know the sign for sin() which is used for mu_y
    ! because sin(-x) = -sin(x)
    dksum_negative = .FALSE.
    if(dksum < 0) dksum_negative = .TRUE.
    dkdif_negative = .FALSE.
    if(dkdif < 0) dkdif_negative = .TRUE.
    dksum = abs(dksum)
    dkdif = abs(dkdif)

    tp1(:) = ttheta1(j1,m1,:)* ttheta1(j2,m2,:)
    tp2(:) =  ttheta(l1,k1,:)*  ttheta(l2,k2,:)


    if (dipole_projection == 'z' .AND. dkr == 0) then
      ! main contribution
      do i=1,nphi; tp(i) = tphi_cos(dkdif,i); end do

      if     (kr1 == 0 .AND. jp1 == 1) then
        do i=1,nphi; tp(i) = tp(i) - tphi_cos(dksum,i); end do
      elseif (kr1 == 0 .AND. jp1 == 0) then
        do i=1,nphi; tp(i) = tp(i) + tphi_cos(dksum,i); end do
      endif

    elseif (dipole_projection == 'x' .AND. abs(dkr) == 1) then
      ! main contribution
      if     (dkr == 1 .AND. mod(jr1+jr2+1,2) == 0) then
        do i=1,nphi; tp(i) =  tphi_cos(dkdif,i); end do
      elseif (dkr == 1 .AND. mod(jr1+jr2+1,2) == 1) then
        do i=1,nphi; tp(i) = -tphi_cos(dkdif,i); end do
      elseif (dkr ==-1) then
        do i=1,nphi; tp(i) = -tphi_cos(dkdif,i); end do
      endif

      ! special cases
      if     (dkr == 1 .AND. kr1 == 0 .AND. mod(jp2,2) == 0) then
        do i=1,nphi; tp(i) = tp(i) + tphi_cos(dksum,i); end do
      elseif (dkr == 1 .AND. kr1 == 0 .AND. mod(jp2,2) == 1) then
        do i=1,nphi; tp(i) = tp(i) - tphi_cos(dksum,i); end do
      elseif (dkr ==-1 .AND. kr1 == 1 .AND. mod(jp2+1,2) == 0) then
        do i=1,nphi; tp(i) = tp(i) + tphi_cos(dksum,i); end do
      elseif (dkr ==-1 .AND. kr1 == 1 .AND. mod(jp2+1,2) == 1) then
        do i=1,nphi; tp(i) = tp(i) - tphi_cos(dksum,i); end do
      endif

    elseif (dipole_projection == 'y' .AND. abs(dkr) == 1) then
      ! main contribution
      if     (dkr == 1 .AND. mod(jr1+jr2,2) == 0) then
        if(dkdif_negative) then
          do i=1,nphi; tp(i) = -tphi_sin(dkdif,i); end do
        else
          do i=1,nphi; tp(i) =  tphi_sin(dkdif,i); end do
        endif
      elseif (dkr == 1 .AND. mod(jr1+jr2,2) == 1) then
        if(dkdif_negative) then
          do i=1,nphi; tp(i) =  tphi_sin(dkdif,i); end do
        else
          do i=1,nphi; tp(i) = -tphi_sin(dkdif,i); end do
        endif
      elseif (dkr ==-1) then
        if(dkdif_negative) then
          do i=1,nphi; tp(i) =  tphi_sin(dkdif,i); end do
        else
          do i=1,nphi; tp(i) = -tphi_sin(dkdif,i); end do
        endif
      endif

      ! special cases
      if     (dkr == 1 .AND. kr1 == 0 .AND. mod(jp2,2) == 0) then
        if(dksum_negative) then
          do i=1,nphi; tp(i) = tp(i) - tphi_sin(dksum,i); end do
        else
          do i=1,nphi; tp(i) = tp(i) + tphi_sin(dksum,i); end do
        endif
      elseif (dkr == 1 .AND. kr1 == 0 .AND. mod(jp2,2) == 1) then
        if(dksum_negative) then
          do i=1,nphi; tp(i) = tp(i) + tphi_sin(dksum,i); end do
        else
          do i=1,nphi; tp(i) = tp(i) - tphi_sin(dksum,i); end do
        endif
      elseif (dkr ==-1 .AND. kr1 == 1 .AND. mod(jp2+1,2) == 0) then
        if(dksum_negative) then
          do i=1,nphi; tp(i) = tp(i) - tphi_sin(dksum,i); end do
        else
          do i=1,nphi; tp(i) = tp(i) + tphi_sin(dksum,i); end do
        endif
      elseif (dkr ==-1 .AND. kr1 == 1 .AND. mod(jp2+1,2) == 1) then
        if(dksum_negative) then
          do i=1,nphi; tp(i) = tp(i) + tphi_sin(dksum,i); end do
        else
          do i=1,nphi; tp(i) = tp(i) - tphi_sin(dksum,i); end do
        endif
      endif

    else
      if(proc_id == 0) write(6,*) " dipolematrix:  unspecified case"
      stop
    endif

    integral = ZERO

      do i3=1,nphi
        do i2=1,nt
          do i1=1,nt1
            integral = integral+ tp1(i1)*tp2(i2)*tp(i3)*mu_grid(i1,i2,i3,is1,is2,is3)
          end do
        end do
      end do

    integral = integral * fct

    dm(ia1,ia2) = integral

  end do
end do

END SUBROUTINE


SUBROUTINE dipole_grid_x(r1,r2,r3)
  USE types
  USE param
IMPLICIT NONE
  REAL(DP), INTENT(IN) :: r1,r2,r3
!  REAL(DP), INTENT(IN) :: theta1(nt1)
!  REAL(DP), INTENT(IN) :: theta(nt)
!  REAL(DP), INTENT(IN) :: phi(nphi)
REAL(DP) :: ppp,mpp,pmp,mmp,ppm,mpm,pmm,mmm
REAL(DP) :: pps,mps,pms,mms,ppa,mpa,pma,mma
REAL(DP) :: pss,mss,pas,mas,psa,msa,paa,maa
REAL(DP) :: mux, muy, muz
INTEGER(I4B) :: i1,i2,i3

!  symmetry index: 0 means sym; 1 means asym
!  runs through theta1, theta and phi respectively

!write(6,*) 'in dipole_grid'

do i1=1, nt1; do i2=1, nt; do i3=1, nphi

    ! mux

    !  compute the full grid
    ppp = mux(r1,r2,r3,theta1(i1),    theta(i2),    phi(i3))
    mpp = mux(r1,r2,r3,PI-theta1(i1), theta(i2),    phi(i3))
    pmp = mux(r1,r2,r3,theta1(i1),    PI-theta(i2), phi(i3))
    mmp = mux(r1,r2,r3,PI-theta1(i1), PI-theta(i2), phi(i3))
    ppm = mux(r1,r2,r3,theta1(i1),    theta(i2),    PI-phi(i3))
    mpm = mux(r1,r2,r3,PI-theta1(i1), theta(i2),    PI-phi(i3))
    pmm = mux(r1,r2,r3,theta1(i1),    PI-theta(i2), PI-phi(i3))
    mmm = mux(r1,r2,r3,PI-theta1(i1), PI-theta(i2), PI-phi(i3))

    !  find sym/asym in the respect to phi
    pps = ppp+ppm;   ppa = ppp-ppm
    pms = pmp+pmm;   pma = pmp-pmm
    mps = mpp+mpm;   mpa = mpp-mpm
    mms = mmp+mmm;   mma = mmp-mmm

    !  find sym/asym in the respect to theta
    pss = pps+pms;  pas = pps-pms;  psa = ppa+pma;  paa = ppa-pma
    mss = mps+mms;  mas = mps-mms;  msa = mpa+mma;  maa = mpa-mma

    !  find sym/asym in the respect to theta1 and store it
    mu_grid(i1,i2,i3,0,0,0) = pss+mss;  mu_grid(i1,i2,i3,1,0,0) = pss-mss
    mu_grid(i1,i2,i3,0,1,0) = pas+mas;  mu_grid(i1,i2,i3,1,1,0) = pas-mas
    mu_grid(i1,i2,i3,0,0,1) = psa+msa;  mu_grid(i1,i2,i3,1,0,1) = psa-msa
    mu_grid(i1,i2,i3,0,1,1) = paa+maa;  mu_grid(i1,i2,i3,1,1,1) = paa-maa

end do; end do; end do

END SUBROUTINE dipole_grid_x


SUBROUTINE dipole_grid_y(r1,r2,r3)
  USE types
  USE param
IMPLICIT NONE
  REAL(DP), INTENT(IN) :: r1,r2,r3
!  REAL(DP), INTENT(IN) :: theta1(nt1)
!  REAL(DP), INTENT(IN) :: theta(nt)
!  REAL(DP), INTENT(IN) :: phi(nphi)
REAL(DP) :: ppp,mpp,pmp,mmp,ppm,mpm,pmm,mmm
REAL(DP) :: pps,mps,pms,mms,ppa,mpa,pma,mma
REAL(DP) :: pss,mss,pas,mas,psa,msa,paa,maa
REAL(DP) :: mux, muy, muz
INTEGER(I4B) :: i1,i2,i3

!  symmetry index: 0 means sym; 1 means asym
!  runs through theta1, theta and phi respectively

!write(6,*) 'in dipole_grid'

do i1=1, nt1; do i2=1, nt; do i3=1, nphi

    ! muy

    !  compute the full grid
    ppp = muy(r1,r2,r3,theta1(i1),    theta(i2),    phi(i3))
    mpp = muy(r1,r2,r3,PI-theta1(i1), theta(i2),    phi(i3))
    pmp = muy(r1,r2,r3,theta1(i1),    PI-theta(i2), phi(i3))
    mmp = muy(r1,r2,r3,PI-theta1(i1), PI-theta(i2), phi(i3))
    ppm = muy(r1,r2,r3,theta1(i1),    theta(i2),    PI-phi(i3))
    mpm = muy(r1,r2,r3,PI-theta1(i1), theta(i2),    PI-phi(i3))
    pmm = muy(r1,r2,r3,theta1(i1),    PI-theta(i2), PI-phi(i3))
    mmm = muy(r1,r2,r3,PI-theta1(i1), PI-theta(i2), PI-phi(i3))

    !  find sym/asym in the respect to phi
    pps = ppp+ppm;   ppa = ppp-ppm
    pms = pmp+pmm;   pma = pmp-pmm
    mps = mpp+mpm;   mpa = mpp-mpm
    mms = mmp+mmm;   mma = mmp-mmm

    !  find sym/asym in the respect to theta
    pss = pps+pms;  pas = pps-pms;  psa = ppa+pma;  paa = ppa-pma
    mss = mps+mms;  mas = mps-mms;  msa = mpa+mma;  maa = mpa-mma

    !  find sym/asym in the respect to theta1 and store it
    mu_grid(i1,i2,i3,0,0,0) = pss+mss;  mu_grid(i1,i2,i3,1,0,0) = pss-mss
    mu_grid(i1,i2,i3,0,1,0) = pas+mas;  mu_grid(i1,i2,i3,1,1,0) = pas-mas
    mu_grid(i1,i2,i3,0,0,1) = psa+msa;  mu_grid(i1,i2,i3,1,0,1) = psa-msa
    mu_grid(i1,i2,i3,0,1,1) = paa+maa;  mu_grid(i1,i2,i3,1,1,1) = paa-maa

end do; end do; end do

END SUBROUTINE dipole_grid_y


SUBROUTINE dipole_grid_z(r1,r2,r3)
  USE types
  USE param
IMPLICIT NONE
  REAL(DP), INTENT(IN) :: r1,r2,r3
!  REAL(DP), INTENT(IN) :: theta1(nt1)
!  REAL(DP), INTENT(IN) :: theta(nt)
!  REAL(DP), INTENT(IN) :: phi(nphi)
REAL(DP) :: ppp,mpp,pmp,mmp,ppm,mpm,pmm,mmm
REAL(DP) :: pps,mps,pms,mms,ppa,mpa,pma,mma
REAL(DP) :: pss,mss,pas,mas,psa,msa,paa,maa
REAL(DP) :: mux, muy, muz
INTEGER(I4B) :: i1,i2,i3

!  symmetry index: 0 means sym; 1 means asym
!  runs through theta1, theta and phi respectively

!write(6,*) 'in dipole_grid'

do i1=1, nt1; do i2=1, nt; do i3=1, nphi

    ! muz

    !  compute the full grid
    ppp = muz(r1,r2,r3,theta1(i1),    theta(i2),    phi(i3))
    mpp = muz(r1,r2,r3,PI-theta1(i1), theta(i2),    phi(i3))
    pmp = muz(r1,r2,r3,theta1(i1),    PI-theta(i2), phi(i3))
    mmp = muz(r1,r2,r3,PI-theta1(i1), PI-theta(i2), phi(i3))
    ppm = muz(r1,r2,r3,theta1(i1),    theta(i2),    PI-phi(i3))
    mpm = muz(r1,r2,r3,PI-theta1(i1), theta(i2),    PI-phi(i3))
    pmm = muz(r1,r2,r3,theta1(i1),    PI-theta(i2), PI-phi(i3))
    mmm = muz(r1,r2,r3,PI-theta1(i1), PI-theta(i2), PI-phi(i3))

    !  find sym/asym in the respect to phi
    pps = ppp+ppm;   ppa = ppp-ppm
    pms = pmp+pmm;   pma = pmp-pmm
    mps = mpp+mpm;   mpa = mpp-mpm
    mms = mmp+mmm;   mma = mmp-mmm

    !  find sym/asym in the respect to theta
    pss = pps+pms;  pas = pps-pms;  psa = ppa+pma;  paa = ppa-pma
    mss = mps+mms;  mas = mps-mms;  msa = mpa+mma;  maa = mpa-mma

    !  find sym/asym in the respect to theta1 and store it
    mu_grid(i1,i2,i3,0,0,0) = pss+mss;  mu_grid(i1,i2,i3,1,0,0) = pss-mss
    mu_grid(i1,i2,i3,0,1,0) = pas+mas;  mu_grid(i1,i2,i3,1,1,0) = pas-mas
    mu_grid(i1,i2,i3,0,0,1) = psa+msa;  mu_grid(i1,i2,i3,1,0,1) = psa-msa
    mu_grid(i1,i2,i3,0,1,1) = paa+maa;  mu_grid(i1,i2,i3,1,1,1) = paa-maa

end do; end do; end do

END SUBROUTINE dipole_grid_z


END PROGRAM
