! 11/05/2009 INK
! dipole moment calculations for any J
! requires: functions mux, muy, muz

! convention: final state x1 <- initial state x2 or x
!

PROGRAM dipole
  USE types
  USE param
  USE workarrays
IMPLICIT NONE

!  parameters
!
INTEGER(I4B), PARAMETER :: u11=11, u12=12

REAL(DP), ALLOCATABLE :: dt(:,:),dm(:,:)
REAL(DP), ALLOCATABLE :: vm1(:,:),vm2(:,:)
REAL(DP), ALLOCATABLE :: e1(:), e2(:), mtmp(:,:)

REAL(DP), ALLOCATABLE :: theta1(:),ttheta1(:,:,:)
REAL(DP), ALLOCATABLE :: theta(:), ttheta(:,:,:)
!REAL(DP), ALLOCATABLE :: phi(:),   tphi(:,:)
REAL(DP), ALLOCATABLE :: phi(:),tphi_cos(:,:),tphi_sin(:,:)
REAL(DP), ALLOCATABLE :: q1(:),q2(:),q3(:)
REAL(DP), ALLOCATABLE :: tp(:),tp1(:),tp2(:),tpx(:),tpy(:)

REAL(DP), ALLOCATABLE :: mux_grid(:,:,:,:,:,:)
REAL(DP), ALLOCATABLE :: muy_grid(:,:,:,:,:,:)
REAL(DP), ALLOCATABLE :: muz_grid(:,:,:,:,:,:)

REAL(DP) :: zpe, tmp
REAL(DP) :: one_over_sqrt2

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
CHARACTER(LEN=80) title
!fname1= 'x0000-J0'
!fname2= 'x0000-J0'

one_over_sqrt2 = one/sqrt(two)

!call time_stamp(6)

! ====================================================================
! INPUT - WAVR4 part
!
!open(u6,file='out.txt')

open(u5,file='input.txt',status='old')

call input_data()

close(u5)

! just set these for backward compatibility with WAVR4
qe(1)=re(1); qe(2)=re(2); qe(3)=re(3); qe(4)=zero; qe(5)=zero; qe(6)=zero

! ====================================================================
! INPUT - dipole part
!

write(u6,'(" INPUT DIPOLE: ",/)')

open(u5,file='inp.dipole.txt',status='old')

write(u6,'(80("="))')

! final state
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) jp1, j_parity1, l_parity1, jl_parity1, jr1, smt
write(u6,'(1x,5i5,A8)') jp1, j_parity1, l_parity1, jl_parity1, jr1, smt

fname1 = 'x'//i2c(jp1)//i2c(j_parity1)//i2c(l_parity1)//i2c(jl_parity1)//'-K'//i2c(jr1)//'-'//smt

if( HHSYM .AND. smt == "all") then
  write(6,*) ' INPUT error: incompatible file name, check HHSYM'
  stop
endif

! initial state
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) jp2, j_parity2, l_parity2, jl_parity2, jr2, smt
write(u6,'(1x,5i5,A8)') jp2, j_parity2, l_parity2, jl_parity2, jr2, smt

fname2 = 'x'//i2c(jp2)//i2c(j_parity2)//i2c(l_parity2)//i2c(jl_parity2)//'-K'//i2c(jr2)//'-'//smt

if( HHSYM .AND. smt == "all") then
  write(6,*) ' INPUT error: incompatible file name, check HHSYM'
  stop
endif

write(u6,'(80("="))')

close(u5)

write(u6,'(/," INPUT: done.",/)')

! ====================================================================
! CHECK INPUT parameters

write(u6,*) " external files to be opened..."
write(u6,*) "   w-f of   final states:  ", fname1
write(u6,*) "   w-f of initial states:  ", fname2

if(j_parity1 > j_parity_max .OR. j_parity2 > j_parity_max) then
  write(u6, '("  INPUT: j_parity is higher than",i3)') j_parity_max
  stop
end if

if(l_parity1 > l_parity_max .OR. l_parity2 > l_parity_max) then
  write(u6, '("  INPUT: l_parity is higher than",i3)') l_parity_max
  stop
end if

if(jl_parity1 > jl_parity_max .OR. jl_parity2 > jl_parity_max) then
  write(u6, '("  INPUT: jl_parity is higher than",i3)') jl_parity_max
  stop
end if


! ====================================================================

! define missing (derived) quantum numbers
p1=mod((jp1+jr1),2)
p2=mod((jp2+jr2),2)

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

!deallocate(t1,t2,t3,oner1,oner2,oner3)

! angular grid points
!    full grid must be: nt1 >= 2*(jmax+1); nt >= 2*(lmax+1); nphi >= kmax+1
!   we use sym reduced: nt1 >= (jmax+1); nt >= (lmax+1); nphi >= (kmax+1+1)/2
!   the latter is because we want 2*nphi >= kmax+1

! minimal grid
    nt1 = jmax+1
    nt  = lmax+1
    nphi = (kmax+1+1)/2
! double grid
!    nt1 = 2*jmax+1
!    nt  = 2*lmax+1
!    nphi = kmax+1
! quad grid
!    nt1 = 4*jmax+1
!    nt  = 4*lmax+1
!    nphi = 2*kmax+1
write(u6,'(" nt1 =",i3," nt  =",i3," nphi=",i3)') nt1,nt,nphi

allocate(theta1(nt1),ttheta1(0:jmax,0:mmax,nt1), &
         theta(nt),ttheta(0:lmax,0:kmax,nt),  &
         phi(nphi),tphi_cos(0:2*kmax,nphi),tphi_sin(0:2*kmax,nphi))

call grid_theta(theta1,ttheta1,jmax,mmax,nt1)
call grid_theta(theta,ttheta,lmax,kmax,nt)
! this is a modified version of grid_phi found in WAVR4
! in addition to cos grid it returns sin grid as well
call grid_phi2(phi,tphi_cos,tphi_sin,kmax,nphi)
! allocate temp arrays to keep products of wf grids
allocate(tp(size(tphi_cos,2)),tp1(size(ttheta1,3)),tp2(size(ttheta,3)))

allocate(mux_grid(nt1,nt,nphi,0:1,0:1,0:1), &
      &  muy_grid(nt1,nt,nphi,0:1,0:1,0:1), &
      &  muz_grid(nt1,nt,nphi,0:1,0:1,0:1))

! ===================================================================

! ===================================================================

! open file with final eigen functions
open(u11,file=fname1//'.dat', &
  status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
if(idiag /=0) then 
  write(6,*) '  cant open ',fname1//'.dat';  stop
end if

! open file with initial eigen functions
open(u12,file=fname2//'.dat', &
  status='OLD',access='SEQUENTIAL',form='UNFORMATTED',iostat=idiag)
if(idiag /=0) then 
  write(6,*) '  cant open ',fname2//'.dat';  stop
end if

! ===================================================================

call time_stamp(6)

read(u11) nstates1
read(u12) nstates2

!nstates1 = 1
!nstates2 = 1

! dipole transitions matrix
allocate(dt(nstates1,nstates2))
allocate(e1(nstates1),e2(nstates2))

dt(:,:)=ZERO

read(u11) e1(:)
read(u12) e2(:)

do kr1=0,jr1

  na1 = krn1(kr1)%size

  allocate(vm1(na1,nstates1))
  allocate(mtmp(na1,nstates2))

  !  THIS IS A TEMPORARY FIX AND MUST BE DONE PROPERLY FOR J>0
  !

  ! do kr2=0,jr2

  kr2 = kr1
  na2 = krn2(kr2)%size

  ! no contributions with delta K greater than 1
  if( abs(kr2-kr1) > 1 ) cycle

  allocate(vm2(na2,nstates2))
  allocate(dm(na1,na2))


  ! do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
  !  replace the nested loops above by one loop
  do i=1,nn(3)*nn(2)*nn(1)
      
    !  unroll true indices
    id = int(int((i-1)/nn(1))/nn(2))+1
    ic = int((i-1)/nn(1))-nn(2)*(id-1)+1
    ib = i-nn(1)*(ic-1)-nn(1)*nn(2)*(id-1)

    ! read final states for i1,i2,i3,K
      do ii=1,nstates1
        read(u11) vm1(:,ii)
        !debug! write(6,*) vm1(:,ii)
      enddo

    ! read initial states
      do ii=1,nstates2
        read(u12) vm2(:,ii)
        !debug! write(6,*) vm2(:,ii)
    enddo


    tmp = zero

    !call dipolematrix(dm, q1(ib),q2(ic),q3(id), theta1,ttheta1, theta,ttheta, phi,tphi_cos,tphi_sin)
    call dipolematrix(q1(ib),q2(ic),q3(id))

    ! compute all possible transition moments
    ! vm1(nstates1,na1)^T x dm(na1,na2) x vm2(na2,nstates2)

    !  C = alpha [A] [B] + beta [C]
    !  gemm ([],[],m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
    !  [] is 'N', 'T' or 'C'
    !  [A] m.k ; A(lda,:)
    !  [B] k.n ; B(ldb,:)
    !  [C] m.n ; C(ldc,:)

    call dgemm('N','N',na1,nstates2,na2,1.D0,dm,na1,vm2,na2,0.D0,mtmp,na1)
    call dgemm('T','N',nstates1,nstates2,na1,1.D0,vm1,na1,mtmp,na1,1.D0,dt,nstates1)

    ! norm testing
    !!call dgemm('T','N',nstates1,nstates2,na1,1.D0,vm1,nstates1,vm2,na1,1.D0,dt,nstates1)
    !do i1=1,nstates1
    !  do i2=1,nstates2
    !    do i3=1,na1
    !      dt(i1,i2)=dt(i1,i2)+vm1(i3,i1)*vm2(i3,i2)
    !    end do
    !  end do
    !end do

    !    do ia2 =1, krn2(kr2)%size
    !    do ia1 =1, krn1(kr2)%size

  end do

  deallocate(dm)
  deallocate(vm2)

  deallocate(vm1,mtmp)

end do

! deallocate all unnecessary arrays

deallocate(indx1,krn1)
deallocate(indx2,krn2)

! TODO: create line frequency list and sort it

! S = (2J+1)(2J'+1) X^2  where dt = X
!     since J=J'=0, S = dt**2

do ii1=1,nstates1
  do ii2=1,nstates2
    write(u6,'(i4,"<-",i4,f12.4,es16.4)') ii1,ii2,abs(e1(ii1)-e2(ii2)),dt(ii1,ii2)**2
  enddo
enddo

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
  write(u6,'(" opt 1: Jacobi vectors",/)')
  muq(1) = mu1*mu2/(mu1+mu2)
  muq(2) = mu3*(mu1+mu2)/(mu1+mu2+mu3)
  muq(3) = mu4*(mu1+mu2+mu3) /msum
else if (abs(opt)==2) then
  !  Radau vectors
  write(u6,'(" opt 2: Radau vectors",/)')
  muq(1) = mu1
  muq(2) = mu2
  muq(3) = mu3
else if (abs(opt)==3) then
  !  Diatom-diatom vectors
!   these are original mu's
  write(u6,'(" opt 3: Diatom-diatom vectors",/)')
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
  write(u6,'(" opt 4: Orthogonal satellite vectors",/)')
  muq(1) = mu1
  muq(2) = mu2
  muq(3) = mu3*mu4/(mu3+mu4)
else
  write(u6,*) ' masses: no such option:', opt
  stop
end if

mud1 = muq(1)
mud2 = muq(2)
mur  = muq(3)

END SUBROUTINE


! INK 16/06/2007
! this is a trimmed down copy of inpu_data() in rv.sym (or WAVR4 proper)

SUBROUTINE input_data()
  USE types
  USE param
IMPLICIT NONE
REAL(DP) :: temp
INTEGER(I4B) :: i, iauto
CHARACTER(LEN=80) title

write(u6,'(" Four Atomic RoVibrational Program",/, &
           & " INPUT: reading input...",/, &
           & 30("="),"< file:  input.txt >",30("="))')
! ======================================================================
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
! ======================================================================
write(u6,'(" INPUT: processing input...",/)')

call masses()
!   take Radau formulae for masses to test 3D angular + rotation block
!   because masses 1-3 will be the same as in the input and an analytical
!   solution is possible. Also set V=0
write(u6,'(" adapted masses:",/,3f14.8,/)')  muq

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

write(u6, '("  j-, l-, (j+l)-parity:")')
if(j_parity_max==0) then
  write(u6, '("  NO:  q1 -> -q1 symmetry (j-parity)")')
else if(j_parity_max==1) then
  write(u6, '("  YES: q1 -> -q1 symmetry (j-parity)")')
else
  write(u6, '("  INPUT: WRONG j_parity_max",i3)') j_parity_max
  stop
end if

if(l_parity_max==0) then
  write(u6, '("  NO:  q2 -> -q2 symmetry (l-parity)")')
else if(l_parity_max==1) then
  write(u6, '("  YES: q2 -> -q2 symmetry (l-parity)")')
else
  write(u6, '("  INPUT: WRONG l_parity_max=",i3)') l_parity_max
  stop
end if

if(jl_parity_max==0) then
  write(u6, '("  NO:  q3 -> -q3 symmetry (jl-parity)")')
else if(jl_parity_max==1) then
  write(u6, '("  YES: q3 -> -q3 symmetry (jl-parity)")')
else
  write(u6, '("  INPUT: WRONG jl_parity_max=",i3)') jl_parity_max
  stop
end if

!  can't have all j_parity_max=l_parity_max=jl_parity_max=1
if(j_parity_max==1 .AND. l_parity_max==1 .AND. jl_parity_max==1) then
  write(u6, '("  INPUT: when all parities = 1, jl_parity is redundant")')
  write(u6, '("  INPUT: jl_parity_max is reset to 0")')
  jl_parity_max=0
end if

write(u6, '(/," angular basis size:              ", &
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
  write(u6,'(" (auto) signle ang grid, iauto=",i2,":")') iauto
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
  write(u6,'(" (auto) double ang grid, iauto=",i2,":")') iauto
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
  write(u6,'(" (auto) 4x     ang grid, iauto=",i2,":")') iauto
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
  write(u6,'(" (input) angular   grid, iauto=",i2,":")') iauto
end if
write(u6,'(" nt1 =",i3," nt  =",i3," nphi=",i3)') nt1,nt,nphi

nagrid = nt1*nt*nphi
write(u6,'(/," full angular grid: nagrid = nt1*nt*nphi =",i6,/)') nagrid

! energy cutoffs and zero energy

write(u6,'(" input zero energy   (enzero)=",f10.1)') enzero

write(u6,'(/," INPUT: done.",/)')

END SUBROUTINE



!SUBROUTINE dipolematrix(dm, r1,r2,r3, theta1,ttheta1, theta,ttheta, phi,tphi_cos,tphi_sin)
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

! ia1,kr1, ia2,kr2  etc must be inherited from the calling routine

! generate symmetrized (sym & anti-sym grids for mux, muy, muz
call dipole_grid(r1,r2,r3)

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
    norm1 = 1._dp/SQRT2
  end if

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
      norm2 = 1._dp/SQRT2
    end if

    integral = ZERO

    fct = norm1*norm2
    if (mod(k1+k2, 2) /= 0) fct = -fct

! since E^* mux =  mux;  mux(-phi) =  mux(phi)
!       E^* muy = -muy;  muy(-phi) = -muy(phi)
!       E^* muz =  muz;  muz(-phi) =  muz(phi)
! i suspect the only non-zero integrals are going to be between states 
! with different parity for muy
! and with same parity for mux & muz

    is1=mod(j1+m1+j2+m2,2)
    is2=mod(l1+k1+l2+k2,2)
    is3=mod(k1+k2,2)

    tp1(:) = ttheta1(j1,m1,:)* ttheta1(j2,m2,:)
    tp2(:) =  ttheta(l1,k1,:)*  ttheta(l2,k2,:)

! DEBUG: if set to .TRUE. this will do just orthogonality test 
!                         final and initial w-f need to be the same
!
if(.FALSE.) then
!if(.TRUE.) then  ! DEBUG

    ! delta K = 0
    !if(kr1 == kr2) then

      if (kr1 == 0) then                       ! special case K=0
        dkp = k1+k2
        dkm = abs(k1-k2)
        if (p1 /= 0) then     !  (-1)**(J+p) factor
          do i=1,nphi; tp(i) = tphi_cos(dkm,i)-tphi_cos(dkp,i); end do
        else
          do i=1,nphi; tp(i) = tphi_cos(dkm,i)+tphi_cos(dkp,i); end do
          ! just a test to see if SIN(phi) is orthogonal do i=1,nphi; tp(i) = tphi_sin(dkp,i); end do
        end if
      else if (kr1 /= 0 .AND. kp1 == kp2) then  ! case I
        dkm = abs(k1-k2)
          do i=1,nphi; tp(i) = tphi_cos(dkm,i); end do
      else                                    ! case II
        dkp = k1+k2
        ! if (mod(jp+kr, 2) /= 0) then     !  (-1)**(J+K+p) factor removed
        !   do i=1,nphi; tp(i) = -tphi(dkp,i); end do
        ! else
          do i=1,nphi; tp(i) =  tphi_cos(dkp,i); end do
        ! end if
      end if

      integral = 0.0D0

      ! due to symmetry the integral is non-zero only if
      !if(is1 == 0 .AND. is2 == 0 .AND. is3 == 0) then
        do i3=1,nphi
          do i2=1,nt
            do i1=1,nt1
              integral = integral+ tp1(i1)*tp2(i2)*tp(i3)*muz_grid(i1,i2,i3,is1,is2,is3)
      !        integral = integral+ tp1(i1)*tp2(i2)*tp(i3)*8.0_dp
            end do
          end do
        end do
      !end if

    dm(ia1,ia2) = 2._dp*fct* integral

    !! write(6,'(" ia1=",i4," ia2=",i4,es18.8)') ia1,ia2,dm(ia1,ia2)
    !debug! write(6,'(10i4,f18.14)') ia1,ia2,j1,m1,l1,k1,j2,m2,l2,k2, dm(ia1,ia2)

! DEBUG END
else

    ! delta K = 0
    if(kr1 == kr2) then

      ! symmetry of cos & sin in the respect to phi= Pi/2 axis:
      !        n= 1  2  3
      ! cos(nx)   a  s  a
      ! sin(nx)   s  a  s
      !
      ! since (-kp1*k1+kp2*k2)  and (kp1*k1+kp2*k2)  must have different parity
      ! the sum of cos(-kp1*k1+kp2*k2) + sin(kp1*k1+kp2*k2) retains symmetry i.r.t. phi=Pi/2

      tmp1= my3j(jr2,kr2,jr1,-kr1)*2.D0
      if(p1 == p2) then
        do i=1,nphi; tp(i) = tphi_cos(-kp1*k1+kp2*k2,i); end do
      else
        do i=1,nphi; tp(i) = tphi_sin(-kp1*k1+kp2*k2,i); end do
      end if

      if(kr1 == 0) then
        if(mod(jr1+p1,2) /= 0) then
          if(p1 == p2) then
            do i=1,nphi; tp(i) = tp(i)- tphi_sin(kp1*k1+kp2*k2,i); end do
          else
            do i=1,nphi; tp(i) = tp(i)- tphi_cos(kp1*k1+kp2*k2,i); end do
          end if
        else
          if(p1 == p2) then
            do i=1,nphi; tp(i) = tp(i)+ tphi_sin(kp1*k1+kp2*k2,i); end do
          else
            do i=1,nphi; tp(i) = tp(i)+ tphi_cos(kp1*k1+kp2*k2,i); end do
          end if
        end if
      end if

      ! straightforward implementation
      !      integral = integral+ tp1(i1)*tp2(i2)*tp(i3)*muz(r1,r2,r3,theta1(i1),    theta(i2),    phi(i3))

      ! alternative implementation using symmetrised grid values

      ! due to muz invertion symmetry
      if(p1 == p2) then
        do i3=1,nphi
          do i2=1,nt
            do i1=1,nt1
              integral = integral+ tp1(i1)*tp2(i2)*tp(i3)*muz_grid(i1,i2,i3,is1,is2,is3)
            end do
          end do
        end do
      end if

    ! delta K = 1
    else 
      if(kr1 == kr2-1) then

        tmp1= my3j(jr2,-kr2,jr1,kr2-1)*SQRT2
        if(mod(jr2+jr1,2) /= 0) tmp1 = -tmp1
        if(p1 == p2) then
          do i=1,nphi; tpy(i) =  tphi_cos(kp1*k1-kp2*k2,i); end do
          do i=1,nphi; tpx(i) =  tphi_sin(kp1*k1-kp2*k2,i); end do
        else
          do i=1,nphi; tpx(i) = -tphi_cos(kp1*k1-kp2*k2,i); end do
          do i=1,nphi; tpy(i) =  tphi_sin(kp1*k1-kp2*k2,i); end do
        end if

      else if(kr1 == 1-kr2) then

        tmp1= my3j(jr2,-kr2,jr1,kr2-1)*SQRT2
        if(mod(jr2+p2+kr2,2) /= 0) tmp1 = -tmp1
        if(p1 == p2) then
          do i=1,nphi; tpy(i) = -tphi_cos(kp1*k1+kp2*k2,i); end do
          do i=1,nphi; tpx(i) =  tphi_sin(kp1*k1+kp2*k2,i); end do
        else
          do i=1,nphi; tpx(i) = -tphi_cos(kp1*k1+kp2*k2,i); end do
          do i=1,nphi; tpy(i) = -tphi_sin(kp1*k1+kp2*k2,i); end do
        end if

      else if(kr1 == kr2+1) then

        tmp1= my3j(jr2,-kr2,jr1,kr2-1)*SQRT2

        if(p1 == p2) then
          do i=1,nphi; tpy(i) = -tphi_cos(-kp1*k1+kp2*k2,i); end do
          do i=1,nphi; tpx(i) = -tphi_sin(-kp1*k1+kp2*k2,i); end do
        else
          do i=1,nphi; tpx(i) = -tphi_cos(-kp1*k1+kp2*k2,i); end do
          do i=1,nphi; tpy(i) =  tphi_sin(-kp1*k1+kp2*k2,i); end do
        end if

      end if

      ! straightforward implementation
      !      integral = integral+ tp1(i1)*tp2(i2)*(  &
      !        &   tpx(i3)* mux(r1,r2,r3,theta1(i1),    theta(i2),    phi(i3))  &
      !        &  +tpy(i3)* muy(r1,r2,r3,theta1(i1),    theta(i2),    phi(i3)))

      ! alternative implementation using symmetrised grid values

      ! due to mux and muy invertion symmetry
      if(p1 == p2) then
        do i3=1,nphi
          do i2=1,nt
            do i1=1,nt1
              integral = integral+ tp1(i1)*tp2(i2)*tpx(i3)* mux_grid(i1,i2,i3,is1,is2,is3)  
            end do
          end do
        end do
      else
        do i3=1,nphi
          do i2=1,nt
            do i1=1,nt1
              integral = integral+ tp1(i1)*tp2(i2)*tpy(i3)* muy_grid(i1,i2,i3,is1,is2,is3)
            end do
          end do
        end do
      end if

    end if

    integral = integral*tmp1

    dm(ia1,ia2) = integral

end if

  end do
end do

END SUBROUTINE


SUBROUTINE dipole_grid(r1,r2,r3)
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
    mux_grid(i1,i2,i3,0,0,0) = pss+mss;  mux_grid(i1,i2,i3,1,0,0) = pss-mss
    mux_grid(i1,i2,i3,0,1,0) = pas+mas;  mux_grid(i1,i2,i3,1,1,0) = pas-mas
    mux_grid(i1,i2,i3,0,0,1) = psa+msa;  mux_grid(i1,i2,i3,1,0,1) = psa-msa
    mux_grid(i1,i2,i3,0,1,1) = paa+maa;  mux_grid(i1,i2,i3,1,1,1) = paa-maa


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
    muy_grid(i1,i2,i3,0,0,0) = pss+mss;  muy_grid(i1,i2,i3,1,0,0) = pss-mss
    muy_grid(i1,i2,i3,0,1,0) = pas+mas;  muy_grid(i1,i2,i3,1,1,0) = pas-mas
    muy_grid(i1,i2,i3,0,0,1) = psa+msa;  muy_grid(i1,i2,i3,1,0,1) = psa-msa
    muy_grid(i1,i2,i3,0,1,1) = paa+maa;  muy_grid(i1,i2,i3,1,1,1) = paa-maa


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
    muz_grid(i1,i2,i3,0,0,0) = pss+mss;  muz_grid(i1,i2,i3,1,0,0) = pss-mss
    muz_grid(i1,i2,i3,0,1,0) = pas+mas;  muz_grid(i1,i2,i3,1,1,0) = pas-mas
    muz_grid(i1,i2,i3,0,0,1) = psa+msa;  muz_grid(i1,i2,i3,1,0,1) = psa-msa
    muz_grid(i1,i2,i3,0,1,1) = paa+maa;  muz_grid(i1,i2,i3,1,1,1) = paa-maa

end do; end do; end do

END SUBROUTINE dipole_grid


END PROGRAM

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Three J symbols. Only
! [ J'    1    J]
! [ M' (-M'-M) M]
! i.e. for dipole calculations
!
! need to make sure that the formulae are valid for M < 0 too.
!
! [j1 j2 j3]       (j1+j2+j3) [j3 j2 j1]
! [m1 m2 m3] = (-1)           [m3 m2 m1]
!
! the formulae are valid for J' >= J.
! therefore if J'=J or J'=J+1 everything is fine but 
! if J'=J-1 we need to swap the columns (j1,m1) and (j3,m3).
! however the factor (-1)^(j1+j2+j3) = (-1)^(J-1+1+J) = 1
! so we actually just swap the number and that is all.

FUNCTION my3j(jj2,mm2,jj,mm)
  USE types
IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: jj2,mm2,jj,mm
  INTEGER(I4B) :: j2,m2,dk,j,m
  REAL(DP) :: my3j, factor

factor = 1.D0

! selection rules
!
! |(-M'-M)| =< 1
!
if (abs(mm2+mm) > 1) then
  write(6,*) " 3J: delta M > 1"
  stop
endif

! |J'-J| <= 1 <= J'+J
!
if (abs(jj2-jj) > 1) then
  write(6,*) " 3J: delta J > 1"
  stop
endif

if (jj2+jj < 1) then
!  write(6,*) " 3J: J'+J > 1"
!  temporary fix
  my3j = 1.D0
  return  
!  stop
endif

! the firmulae works only for J' >= J
! swap (j2,m2) and (j,m) if j > j2
if (jj > jj2) then
  j2=jj
  m2=mm
  j=jj2
  m=mm2
else
  j2=jj2
  m2=mm2
  j=jj
  m=mm
endif

dk = -m2-m

if (j2==j) then
  if (dk==0) then
    my3j = m/(sqrt(j*(j+1.D0)*(2*j+1.D0)))
    if(mod(j-m,2) /=0) factor= -factor
  elseif (dk==  1) then
    my3j = sqrt((j-m)*(j+m+1.D0)/(2.D0*j*(j+1.D0)*(2*j+1.D0)))
    if(mod(j-m,2) /=0) factor= -factor
  elseif (dk== -1) then
    my3j = sqrt((j+m)*(j-m+1.D0)/(2.D0*j*(j+1.D0)*(2*j+1.D0)))
    if(mod(j-m+1,2) /=0) factor= -factor
  endif
else ! j2 > d
  if (dk==0) then
    my3j = sqrt(2.D0*(j+m+1.D0)*(j-m+1.D0)/((2.D0*j+1.D0)*(2.D0*j+2.D0)*(2.D0*j+3.D0)))
    if(mod(j-m+1,2) /=0) factor= -factor
  elseif (dk==  1) then
    my3j = sqrt((j+m+1.D0)*(j+m+2.D0)/((2.D0*j+1.D0)*(2.D0*j+2.D0)*(2.D0*j+3.D0)))
    if(mod(j-m,2) /=0) factor= -factor
  elseif (dk== -1) then
    my3j = sqrt((j-m+1.D0)*(j-m+2.D0)/((2.D0*j+1.D0)*(2.D0*j+2.D0)*(2.D0*j+3.D0)))
    if(mod(j-m,2) /=0) factor= -factor
  endif
endif

my3j = factor* my3j

END FUNCTION
