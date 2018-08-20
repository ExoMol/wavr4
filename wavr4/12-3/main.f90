! WAVR4 Program: Wide Amplitude Vibration-Rotation 4-atomic code
!
! 24/11/02  final release.  truncation order: (1,2),3
!
! 3D stretch (DVR) + 3D angles (FBR) + rotation
!
! optimization_flag =   : redundant
!
! expansion_flag = 0 : compute angular matrix elements using quadrature
! expansion_flag = 1 : compute angular matrix elements thru the
!                          expansion of potential energy function at every
!                          stretching grid point (stored as array)
!
! common quantum numbers:
! j, m; l, k; jr=J, kr=K, p=parity
! jp is parity of J+p:
!    if jp=0 => J+p is even
!    if jp=1 => J+p is odd
!
! i1/j1/l1 is row/first/left index/quantum number
! i/j/l   is column/second/right index/quantum number
!
! stage specific arrays and variables
!
! stage 1:
!   istage1   number of selected eigenvectors
!   eigvec1   array of eigenvectors and eigenvalues
!   rec_num1  number of the (eigenvector) record as referenced by
!             the radial grid i1,i2,i3 and order number i in angular basis
!   ist_max1  maximum possible i (i.e. max selected number of angular 
!             functions). If known in advance, it helps to allocate smaller 
!             rec_num1 size
!   idx1      helps to reference stage1 eigenstates from within stage2.
!             For each radial triple point it indicates:
!   idx1()%size   number of selected eigenstates
!   idx1()%start  comulative order number of the first eigenstate
!   idx1()%end    comulative order number of the last eigenstate
!   Note that idx1(nn(1),nn(2),nn(3))%end = istage1
!
! stage 2:
!   istage2   number of selected eigenvectors
!   eigvec2   array of eigenvectors and eigenvalues
!   rec_num2  number of the (eigenvector) record as referenced by
!             the radial grid i3 and order number i in angular basis
!   ist_max2  maximum possible i (i.e. max selected number of angular 
!             functions). If known in advance, it helps to allocate smaller 
!             rec_num2 size
!   idx2      helps to reference stage2 eigenstates from within stage3.
!             For each radial point it indicates:
!   idx2()%size   number of selected eigenstates
!   idx2()%start  comulative order number of the first eigenstate
!   idx2()%end    comulative order number of the last eigenstate
!   Note that idx2(nn(3))%end = istage2
!
! stage 3:
!   istage3   number of selected eigenvectors
!   eigvec3   array of eigenvectors and eigenvalues
!   idx3      helps to reference stage3 eigenstates from within stage4.
!             For each K it indicates:
!   idx3()%size   number of selected eigenstates
!   idx3()%start  comulative order number of the first eigenstate
!   idx3()%end    comulative order number of the last eigenstate
!
! some arrays:
! grids:
!   q1, q2, q3, theta1, theta, phi
! precomputed at grid basis functions:
!   t1, t2, t3, ttheta1, ttheta, tphi
!
!   indx  angular quantum numbers for every basis function
!   idx   3 numbers to identify stretching grid and 1 angular basis
!   krn   number of angular basis functions for every K
!   kri   number of angular basis functions for every K in reduced basis
!         since reduced_basis option has been removed it is now redundant
!         and set equal to krn. TO BE REMOVED
!   vex   angular potential V expansion
!
!   h3, e3  matrix and energy used for the 1st stage problem (angular)
!   eigvec1  eigenvectors + enrgies of the 1st stage
!   h5, e5  matrix and energy used for the 2nd stage problem
!   eigvec2  eigenvectors + enrgies of the 2nd stage
!   h6, e6  matrix and energy used for the final problem
!
!   nsmax  full streching basis (grid)
!   na     full angular basis for given J,p

PROGRAM rv4
  USE types
  USE potential
  USE param
  USE workarrays
  USE grid3d
  USE laguer_dvr
  USE angles
  USE expansion
  USE base_lib, ONLY: time_stamp, diag, diag_p, xn
IMPLICIT NONE
!
!  common arrays
!
REAL(DP), ALLOCATABLE :: theta1(:),ttheta1(:,:,:)
REAL(DP), ALLOCATABLE :: theta(:), ttheta(:,:,:)
REAL(DP), ALLOCATABLE :: phi(:),   tphi(:,:)
REAL(DP), ALLOCATABLE :: q1(:),q2(:),q3(:)
REAL(DP), ALLOCATABLE :: h3(:,:),e3(:),h5(:,:),e5(:)
REAL(DP), ALLOCATABLE :: h6(:),e6(:),vectors(:,:)
REAL(DP), ALLOCATABLE :: vex(:,:,:,:,:,:),h3ta(:)
!
!  stage specific arrays and variables: stages 1, 2, 3
!
TYPE(matrixarray), ALLOCATABLE :: vm(:,:,:)
REAL(DP), POINTER :: vec11(:),vec12(:), vec21(:),vec22(:)
TYPE(eigenvector), ALLOCATABLE :: eigvec1(:), eigvec2(:) ! , eigvec3(:)
INTEGER(I4B), ALLOCATABLE :: rec_num1(:,:,:,:),rec_num2(:,:)
INTEGER(I4B) :: istage1,istage2   ! ,istage3
INTEGER(I4B) :: ist1_max,ist2_max ! ,ist3_max
!
!  auxiliary variables
!
REAL(DP) :: temp,eshift, eprev,emin,en0  ! memsize
INTEGER(I4B) :: n,n1, i1,i2,i3, nlev,nout,idiag, nf,nnmax
INTEGER(I4B) :: ia,ia1,ib,ib1,ic,ic1,id,id1,im,im1
INTEGER(I4B) :: i,l,j, kr, icounter
! CHARACTER*2 :: i2c
CHARACTER*15 :: fname
CHARACTER*5 :: fpar
!
!  logical array storing the results of test_potential()
!  stage 1 calculation at a grid point will NOT be performed
!  if Vmin > Encut1 for that grid point => testpot = .TRUE.
!
LOGICAL(LGT), ALLOCATABLE :: testpot(:,:,:)

open(5,file='input.txt',status='old')
!open(5,file='input1.txt',status='old')
open(6,file='output.txt')
! open(7,file='output7.txt')
! open(8,file='output8.txt')

write(u6,'(" WAVR4 version 1.0: (1,2),3 h6swap",/)')

!  temporary external storage for h6
open(unit=u16,file='x_h6.dat', form='UNFORMATTED')
! open(18,file='data.txt')

!------------------------------------------------------
!  Ar2HF data
!open(15,file='ar2hf1.txt',status='old')
!temp=v_init()
!------------------------------------------------------

call input_data()
!  find max number of angular states (namax) and mmax
call angular_states_max()

allocate(testpot(nn(1),nn(2),nn(3)))
allocate(t1(nn(1),nn(1)),t2(nn(2),nn(2)),t3(nn(3),nn(3)), &
         oner1(nn(1),nn(1)),oner2(nn(2),nn(2)),oner3(nn(3),nn(3)), &
         q1(nn(1)),q3(nn(3)),q2(nn(2)), &
         idx1(nn(1),nn(2),nn(3)),idx2(nn(3)),idx3(0:krmax))

call radial_grids(q1,q2,q3,t1,t2,t3,oner1,oner2,oner3)

call time_stamp(u6)
if(expansion_flag == 0) then
  !  allocate and compute angular grid
  allocate(theta1(nt1),ttheta1(0:jmax,0:mmax,nt1), &
         theta(nt),ttheta(0:lmax,0:kmax,nt),  &
         phi(nphi),tphi(0:2*kmax,nphi))

  call angular_grids(theta1,theta,phi,ttheta1,ttheta,tphi)

  testpot = .FALSE.
  ! test if a calculation on a radial grid point needs to be done
  do i3=1,nn(3); do i2=1,nn(2); do i1=1,nn(1)
    if(test_potential(q1(i1),q2(i2),q3(i3), theta1,theta,phi)) then
      write(6,'(" Vmin > Encut1 for grid point",3i4)') i1,i2,i3
      testpot(i1,i2,i3) = .TRUE.
    end if
  end do; end do; end do

else ! if(expansion_flag == 1) then
  !  expansion of the angular potential in spherical harmonics
  !  for every radial triple point
  !  since 3J symbols will be needed later, allocate & precompute binom
  !
  !  quick fix to write vex to disk and read it later
  !
  !  NEED to check whether this is really MAX nbin possible
    nbin1 = max(jmax,lmax)
    nbin2 = max(ne1,ne2)
  write(6,*) '  nbin1=',nbin1,'  nbin2=',nbin2
  allocate(threej0(0:nbin2,0:nbin1,0:nbin1),tp1(0:nbin2),tp2(0:nbin2))
  allocate(vex(0:ne1,0:ne2,0:ne3,nn(1),nn(2),nn(3)))

  call setfac()
  call vexpand(vex,q1,q2,q3,testpot)

  ! storage of expansion terms
  n=(ne1+1)*(ne2+1)*(ne3+1)
  if(mod(n,2)==1) n=n+1
  n= rec_factor*n
  open(unit=u10,file='x_expansion.dat',access='DIRECT',recl=n,form='UNFORMATTED')
  ! open(unit=u10, form='UNFORMATTED')
  do i3=1,nn(3); do i2=1,nn(2); do i1=1,nn(1)
    write(u10,rec=(i1+(i2-1)*nn(1)+(i3-1)*nn(1)*nn(2))) vex(:,:,:,i1,i2,i3)
  end do; end do; end do
  ! endfile u10

  deallocate(vex)
  !  this will be much smaller temporary storage when reading in vex
  allocate(vex(0:ne1,0:ne2,0:ne3,1,1,1))

end if

write(u6,'(" preparatory stage")')
call time_stamp(u6)

jp_loop: do jp=0,1

j_parity_loop: do j_parity =0,j_parity_max
l_parity_loop: do l_parity =0,l_parity_max
jl_parity_loop: do jl_parity =0,jl_parity_max

  ! skip A1 symmetry for Ar2-HF
  ! if(jp==0 .AND. j_parity==0) cycle

  ! skip B2 symmetry for Ar2-HF
  ! if(jp==0 .AND. j_parity==1) cycle

  ! skip A2 symmetry for Ar2-HF
  ! if(jp==1 .AND. j_parity==0) cycle

  ! skip B1 symmetry for Ar2-HF
  ! if(jp==1 .AND. j_parity==1) cycle

  write(u6,'(//,"   %%%%     jp=",i2,"    %%%%")') jp
  if(j_parity_max > 0)  write(u6,'("   %%%%  j_parity=",i2," %%%%")') j_parity
  if(l_parity_max > 0)  write(u6,'("   %%%%  l_parity=",i2," %%%%")') l_parity
  if(jl_parity_max > 0) write(u6,'("   %%%%  jl_parity=",i2," %%%%")') jl_parity
  write(u6,'(/)')

  !  allocate and clean auxiliary arrays:
  allocate(indx(namax,6),krn(0:krmax))

  !  find how many and what possible angular states are: fill in indx(:,6)
  call angular_states_qnumbers(indx,krn)

  kr_loop: do kr=0,krmax

    ! to keep consistency with old code
    ! J and K are effectively the same in this loop
    jr = kr
    p = abs(mod(jr,2)-jp)

    write(u6,'(/,"   %%%%  K=",i2,"  p=",i2,"  %%%%",/)') kr,p
    call time_stamp(u6)

    !  angular sub-block
    na  = krn(kr)%size

    write(u6,*) '  angular K sub-block: na =', na

    ! setup the files to store eigenvalues and eigenvectors
    ! the file names change according to their symmetry and K value

    n=jp
    ! the line below allows to re-use K>0 p=0 eigenvectors obtained for p=1
    if(kr>0) n=0

    ! fpar = 'x'//i2c(jp)//i2c(j_parity)//i2c(l_parity)//i2c(jl_parity)
    fpar = 'x'//i2c(n)//i2c(j_parity)//i2c(l_parity)//i2c(jl_parity)

    ! if(mod(n,2)==1) n=n+1
    ! n= rec_factor*n
    ! open(unit=u11(jr),file=fname,access='DIRECT',recl=n,form='UNFORMATTED')
    ! write(u6,'(" opened direct file with RECL=",i10)') n

    ! stage_flag number identifies what stage to begin with computing
    ! (i.e. all previous stages must have been already computed).
    ! stage_flag = 0 means all stages must be computed regardless of the
    ! file existence.

    ! stage 1 eigenvectors file
    fname = fpar//'-K'//i2c(jr)//'-1.dat'
    if(stage_flag < 1) then
        write(u6,*) ' stage1 vectors will be computed'
        open(unit=u11(jr),file=fname,access='SEQUENTIAL',form='UNFORMATTED')
    else !  attempt to find if stage1 vectors exist
      write(u6,*) ' attempting to open stage1 file'
      open(unit=u11(jr),file=fname,status='OLD',access='SEQUENTIAL', &
        form='UNFORMATTED',iostat=idiag)
      if(idiag==0) then
        write(u6,*) ' stage1 vectors seem to exist'
      else
        write(u6,*) ' unable to open file'
        stop
      end if
    end if

    ! stage 2 eigenvectors file
    fname = fpar//'-K'//i2c(jr)//'-2.dat'
    if(stage_flag < 2) then
        write(u6,*) ' stage2 vectors will be computed'
        open(unit=u21(jr),file=fname,access='SEQUENTIAL',form='UNFORMATTED')
    else !  attempt to find if stage2 vectors exist
      write(u6,*) ' attempting to open stage2 file'
      open(unit=u21(jr),file=fname,status='OLD',access='SEQUENTIAL', &
        form='UNFORMATTED',iostat=idiag)
      if(idiag==0) then
        write(u6,*) ' stage2 vectors seem to exist'
      else
        write(u6,*) ' unable to open file'
        stop
      end if
    end if

    ! stage 3 eigenvectors file
    fname = fpar//'-K'//i2c(jr)//'-3.dat'
    if(stage_flag < 3) then
        write(u6,*) ' stage3 vectors will be computed'
        open(unit=u31(jr),file=fname,access='SEQUENTIAL',form='UNFORMATTED')
    else !  attempt to find if stage3 vectors exist
      write(u6,*) ' attempting to open stage3 file'
      open(unit=u31(jr),file=fname,status='OLD',access='SEQUENTIAL', &
        form='UNFORMATTED',iostat=idiag)
      if(idiag==0) then
        write(u6,*) ' stage3 vectors seem to exist'
      else
        write(u6,*) ' unable to open file'
        stop
      end if
    end if

    if(stage_flag > 2) then
      read(u31(kr)) nout,nlev
      rewind(u31(kr))
      idx3(kr)%size = nout
      if(kr==0) then
        idx3(kr)%start= 1
        idx3(kr)%end = idx3(kr)%start+idx3(kr)%size-1
        cycle  kr_loop
      else
        idx3(kr)%start= idx3(kr-1)%end+1
        idx3(kr)%end = idx3(kr)%start+idx3(kr)%size-1
      end if
    end if

    if( kr>0 ) then
      fpar = 'x'//i2c(jp)//i2c(j_parity)//i2c(l_parity)//i2c(jl_parity)

      ! stage 4 eigenvectors file
      fname = fpar//'-K'//i2c(jr)//'-4.dat'
      if(stage_flag < 4) then
          write(u6,*) ' stage4 vectors will be computed'
          open(unit=u41(jr),file=fname,access='SEQUENTIAL',form='UNFORMATTED')
      else !  attempt to find if stage4 vectors exist
        write(u6,*) ' attempting to open stage4 file'
        open(unit=u41(jr),file=fname,status='OLD',access='SEQUENTIAL', &
          form='UNFORMATTED',iostat=idiag)
        if(idiag==0) then
          write(u6,*) ' stage4 vectors seem to exist'
        else
          write(u6,*) ' unable to open file'
          stop
        end if
      end if

    end if

    if(stage_flag > 3) then
      !  we have nothing to do
      cycle  kr_loop
    else if(stage_flag==3) then  ! if K=0 we should have already cycled
      call mainj()
      cycle kr_loop
    end if

! ===================================================================

    write(u6,'(" starting stage 1...",/)')

    emin = 0.0_dp
    if(stage_flag >= 1) then

      read(u11(jr)) istage1, ist1_max
      if(istage1 < 1) then
        write(u6,'("  skipping: istage1 < 1")')
        ! go to the next J
        cycle kr_loop
      end if

      write(6,*) ' istage1=', istage1
      read(u11(jr)) idx1
      allocate(eigvec1(istage1),rec_num1(ist1_max,nn(1),nn(2),nn(3)))
      rec_num1 = 0

      write(u6,'("  reading in stage 1 eigenvectors...")')

      ! retain only the eigenvectors conforming  Ecut1
      allocate(e3(na))
      icounter = 0
      do i3=1,nn(3)
        do i2=1,nn(2)
          do i1=1,nn(1)

            if(idx1(i1,i2,i3)%size < 1) cycle
            nlev=idx1(i1,i2,i3)%size
            idx1(i1,i2,i3)%size=0
            do ia=1,nlev
              read(u11(jr)) en0,e3

            !  selection algorithm: takes only the sates below encut1

              ! if(en0 > margin1) cycle
              if(select(ia,en0,encut1,margin1,imargin1)) cycle
              idx1(i1,i2,i3)%size=idx1(i1,i2,i3)%size+1
              icounter = icounter+1
              rec_num1(ia,i1,i2,i3)=icounter
              allocate(eigvec1(icounter)%vec(na), &
                       eigvec1(icounter)%eng)
              eigvec1(icounter)%eng = en0
              eigvec1(icounter)%vec(:) = e3
              if(test_flag >1)  &
                write(u8,'(i6,4i4,f12.3)') icounter, i1,i2,i3, ia, en0
            end do

          end do
        end do
      end do
      deallocate(e3)
      call setup_idx1()
      write(u6,'("  read in",i8," eigenvectors...")') icounter
      istage1= icounter

      ! overwrite stored eigenvectors to bring them in accord
      rewind(u11(jr))
      write(u11(jr)) istage1,ist1_max
      write(u11(jr)) idx1
      do n=1,istage1
        write(u11(jr)) eigvec1(n)%eng,eigvec1(n)%vec
      end do

    else

      idx1(:,:,:)%size=0
      !  allocate and clean auxiliary arrays and
      allocate(h3(na,na),e3(na), &
        eigvec1(nsmax*na),rec_num1(na,nn(1),nn(2),nn(3)))

      !  call to angleT3D is moved from angle3d and angle3de to here
      !  so that angular kinetic matrix can be saved and re-used
      !  to revert things back (e.g. if memory is low) bypass h3ta
      !  and move call to angleT3D inside of the triple loop below

      allocate(h3ta(1:(na*(na+1))/2))

      !  h3(:,:) = zero  ! it will be cleaned up inside angleT3D anyway
      call angleT3D(h3,indx, &
        krn(kr)%start,krn(kr)%size,krn(kr)%start,krn(kr)%size)
      do n=1,na; do n1=1,n
        h3ta(xn(n1,n))=h3(n1,n)
      end do; end do

      rec_num1 = 0

      if(expansion_flag == 0) then
        allocate(tp1(nt1),tp2(nt),tp(nphi),v3(nt1,nt,nphi,0:1,0:1,0:1))
      else ! if(expansion_flag == 1) then
        rewind(u10)
      end if

      icounter = 0; ist1_max=0; emin = huge(one); eprev = -huge(one)

      do i3=1,nn(3)
        do i2=1,nn(2)
          do i1=1,nn(1)

            if(expansion_flag == 1) then
              read(u10,rec=(i1+(i2-1)*nn(1)+(i3-1)*nn(1)*nn(2))) vex(:,:,:,1,1,1)
            end if

            ! test if the calculation on this radial grid point is required

            if(testpot(i1,i2,i3)) cycle

            !  initiate angular kinetic matrix
            do n=1,na; do n1=1,n
              h3(n1,n)=h3ta(xn(n1,n))
            end do; end do
            !  alternatively call angleT3D here
            ! call angleT3D(h3,indx, &
            !   krn(kr)%start,krn(kr)%size,krn(kr)%start,krn(kr)%size)

            ! idx1(i1,i2,i3)%size=0
            ! call fdr(q1(i1),q2(i2),q3(i3),frr,fd1r,fd2r)
            eshift = oner3(i3,i3)*jr*(jr+1.0_dp)

            if(expansion_flag == 0) then
              call angle3d(q1(i1),q2(i2),q3(i3), &
                oner1(i1,i1),oner2(i2,i2),oner3(i3,i3), &
                h3,e3, eshift,nlev,  &
                theta1,ttheta1, theta,ttheta, phi,tphi, indx)
            else ! if(expansion_flag == 1) then
              call angle3de(q1(i1),q2(i2),q3(i3), &
                oner1(i1,i1),oner2(i2,i2),oner3(i3,i3), &
                h3,e3, eshift,nlev, &
                vex(:,:,:,1,1,1), indx)
              ! uncomment this if _full_ vex is stored in memory
              ! vex(:,:,:,i1,i2,i3), indx)
            end if

            if(angular_problem_only_flag ==1) cycle

            do ia=1,nlev
              ! if(e3(ia) > encut1) exit
              if(select(ia,e3(ia),encut1,margin1,imargin1)) exit
              idx1(i1,i2,i3)%size=idx1(i1,i2,i3)%size+1
              icounter = icounter+1
              rec_num1(ia,i1,i2,i3)=icounter
              allocate(eigvec1(icounter)%vec(na), &
                       eigvec1(icounter)%eng)
              eigvec1(icounter)%eng = e3(ia)
              eigvec1(icounter)%vec(:) = h3(:,ia)
              if(test_flag >1)  &
                write(u8,'(i6,4i4,f12.3)') icounter, i1,i2,i3, ia, e3(ia)
            end do

            if(emin>e3(1)) emin=e3(1)

            ist1_max= max(ist1_max,idx1(i1,i2,i3)%size)
            ! if(idx1(i1,i2,i3)%size == icut1 .AND. icut1 /= 0) write(u6, &
            !   & '(" warning: max number of levels has been reached:", &
            !   & 3i4,f12.3)') i1,i2,i3,e3(nlev)

          end do
        end do
      end do

      write(u6,'(" stage 1: made",i8," records")') icounter
      call time_stamp(u6)
      istage1 = icounter

      call setup_idx1()

      if(test_flag >0) then
        write(u7,'("  direct product basis=",i10)') na*nn(1)*nn(2)*nn(3)
        write(u7,'("  optimized      basis=",i10)') istage1
        write(u7,'("  effective energy window for optimised basis was",/, &
                 & "    from",f16.3," to",f16.3)') emin,encut1
        ! that must be the minimum number of requested eigenvalues
        write(u7,'("  max number of selected levels=",i10)') ist1_max
      end if
  
      if(expansion_flag == 0) then
        deallocate(tp1,tp2,tp,v3)
      end if

      deallocate(h3,e3,h3ta)

      !  angular problem only  => cycle
      if(angular_problem_only_flag == 1) then
        !  no eigvec1's allocated, so we don't dealloc them
        deallocate(eigvec1,rec_num1)
        call time_stamp(u6)
        cycle
      end if

      write(u11(jr)) istage1,ist1_max

      ! if encut1 is too low => no levels => cycle
      if(istage1 < 1) then
        !  no eigvec1's allocated, so we don't dealloc them
        deallocate(eigvec1,rec_num1)
        call time_stamp(u6)
        cycle
      end if

      write(u11(jr)) idx1
      do n=1,istage1
        write(u11(jr)) eigvec1(n)%eng,eigvec1(n)%vec
      end do

      write(6,*) ' istage1=', istage1

      call time_stamp(u6)
      write(u6,'(" end of stage 1",/)')
  
    end if

! ===================================================================

!  at this point we can still choose the order of further
!  diagonalization-contraction. we choose:
!  i3 -> id
!  i2 -> ic \ these will be treated together in the next step
!  i1 -> ib /

    write(u6,'(" starting stage 2...",/)')

    if(stage_flag >= 2) then
      read(u21(jr)) istage2,ist2_max
      if(istage2 < 1) then
        write(u6,'("  skipping: istage2 < 1")')
        ! go to the next J
        cycle kr_loop
      end if

      write(6,*) ' istage2=', istage2
      read(u21(jr)) idx2
      allocate(eigvec2(istage2),rec_num2(ist2_max,nn(3)))
      rec_num2 = 0

      write(u6,'("  reading in stage 2 eigenvectors...")')

      ! retain only the eigenvectors conforming  Ecut2
      icounter = 0
      do id =1,nn(3)

        if(idx2(id)%size < 1) cycle
        nnmax = idx1(nn(1),nn(2),id)%end
        if(nnmax < 1) cycle
        allocate(e3(nnmax))
        nlev=idx2(id)%size
        idx2(id)%size=0
        do n=1,nlev
          read(u21(jr)) en0,e3

          ! if(en0 > encut2) cycle
          if(select(n,en0,encut2,margin2,imargin2)) cycle
          idx2(id)%size = idx2(id)%size+1
          icounter = icounter+1
          rec_num2(n,id) = icounter
          allocate(eigvec2(icounter)%vec(nnmax), &
                   eigvec2(icounter)%eng)
          eigvec2(icounter)%eng = en0
          eigvec2(icounter)%vec = e3
          if(test_flag >1)  &
            write(u8,'(i6,2i4,f12.3)') icounter, id, n, en0
        end do

        deallocate(e3)
      end do
      call setup_idx2()
      write(u6,'("  read in",i8," eigenvectors...")') icounter
      istage2=icounter

      ! overwrite stored eigenvectors to bring them in accord
      rewind(u21(jr))
      write(u21(jr)) istage2,ist2_max
      write(u21(jr)) idx2
      do icounter=1,istage2
        write(u21(jr)) eigvec2(icounter)%eng,eigvec2(icounter)%vec
      end do

    else

      !  find out upper maximum size of records per i3 coordinate
      nf=0
      do i3=1,nn(3)
        nf=max(nf,idx1(nn(1),nn(2),i3)%end)
      end do
      write(6,*) ' max number of records per i3, nf =', nf

      idx2(:)%size = 0
      !  since we do not know how many recodrs we might need
      !  our safe upper guess will be istage1 >= istage2
      allocate(eigvec2(istage1),rec_num2(nf,nn(3)))
      rec_num2 = 0
      icounter = 0; ist2_max = 0; emin = huge(one); eprev = -huge(one)

      do id =1,nn(3)
        ! write(u6,'(" id=",i4)') id
        ! idx2(id)%size = 0
        nnmax = idx1(nn(1),nn(2),id)%end
        if(nnmax < 1) cycle
        allocate(h5(nnmax,nnmax),e5(nnmax))

        ! write(u6,'(" initiate")')
        h5(:,:) = zero
  
        ! loop over columns
        do ic =1,nn(2)
          do ib =1,nn(1)
            do ia =1,idx1(ib,ic,id)%size

              ! read(u11,rec=rec_num1(ia,ib,ic,id)) en2,vec12(:)
              en0  =  eigvec1(rec_num1(ia,ib,ic,id))%eng
              vec12 => eigvec1(rec_num1(ia,ib,ic,id))%vec(:)
              n = idx1(ib,ic,id)%start+ia-1

              ! diagonal contribution
              h5(n,n) = en0 + t1(ib,ib)+t2(ic,ic)

              ! loop over rows:
              ! off-diag contribution from t1, therefore
              ic1 = ic

              do ib1 =1,ib-1
                do ia1 =1,idx1(ib1,ic1,id)%size
                  !   read(u11,rec=rec_num1(ia1,ib1,ic1,id)) en1,vec11(:)
                    vec11 => eigvec1(rec_num1(ia1,ib1,ic1,id))%vec(:)
                  n1 = idx1(ib1,ic1,id)%start+ia1-1
                  h5(n1,n)=h5(n1,n)+ dot_product(vec11,vec12)*t1(ib1,ib)
                end do
              end do

              ! off-diag contribution from t2, therefore
              ib1 = ib

              do ic1 =1,ic-1
                do ia1 =1,idx1(ib1,ic1,id)%size
                  !   read(u11,rec=rec_num1(ia1,ib1,ic1,id)) en1,vec11(:)
                    vec11 => eigvec1(rec_num1(ia1,ib1,ic1,id))%vec(:)
                  n1 = idx1(ib1,ic1,id)%start+ia1-1
                  h5(n1,n)=h5(n1,n)+ dot_product(vec11,vec12)*t2(ic1,ic)
                end do
              end do

            end do
          end do

        end do

        ! write(u6,'(" diagonalizing")')
        nlev= nnmax
        call diag(h5,e5,nlev,nout,icut2,encut2,2,idiag)
        if(idiag /= 0) then
          write(6,*) ' diag failed to diagonalise on stage 2'
          write(6,*) ' trace: kr=',kr,' id=',id,' id1=',id1
          stop
        end if

        ! write(u6,'(/," h5 energies:")')

        do n=1,nout
          ! if(e5(n) > encut2) exit
          if(select(n,e5(n),encut2,margin2,imargin2)) exit
          idx2(id)%size = idx2(id)%size+1
          icounter = icounter+1
          rec_num2(n,id) = icounter
          allocate(eigvec2(icounter)%vec(nnmax), &
                   eigvec2(icounter)%eng)
          eigvec2(rec_num2(n,id))%eng = e5(n)
          eigvec2(rec_num2(n,id))%vec = h5(:,n)
          if(test_flag >1)  &
           write(u6,'(i6,2i4,f12.3)') icounter,id,n,e5(n)
        end do

        if(emin>e5(1)) emin=e5(1)

        ist2_max= max(ist2_max,idx2(id)%size)
        ! if(idx2(id)%size == icut2 .AND. icut2 /= 0) write(u6, &
        !   & '(" warning: max number of levels has been reached:", &
        !   & i6,f12.3)') id,e5(nout)

        deallocate(h5,e5)

      end do   !  end of "id"  loop

      call setup_idx2()

      if(test_flag >0) then
        write(u7,'("  effective energy window for optimised basis was",/, &
                 & "    from",f16.3," to",f16.3)') emin,encut2
        ! that must be the minimum number of requested eigenvalues
        write(u7,'("  max number of selected levels=",i10)') ist2_max
      end if

      istage2 = idx2(nn(3))%end
      write(u6,'(" stage 2: made",i8," records")') istage2

      write(u21(jr)) istage2,ist2_max

      ! if encut2 is too low => no levels => cycle
      if(istage2 < 1) then
        !  no eigvec2's allocated, so we don't dealloc them
        deallocate(eigvec2,rec_num2)
        cycle
      end if

      write(u21(jr)) idx2
      do icounter=1,istage2
        write(u21(jr)) eigvec2(icounter)%eng,eigvec2(icounter)%vec
      end do

      call time_stamp(u6)

    end if

    write(u6,'(" end of stage 2",/)')

! ===================================================================

    write(u6,'(" starting stage 3...",/)')
    write(u6,'("   matrix size=",i5)') istage2

    ! h6 setup in packed storage: only Upper triangular part columnwise
    !   only one column of h6 is used to create full h6 on the disk
    !   and then read it in
    rewind(u16)
    ! allocate(h6(istage2,istage2))
    allocate(h6(istage2))

    !  temporary need h3 & e3 again
    allocate(h3(na,na),e3(na))

    h3(:,:) = zero
    ! call angleT3D(h3,indx,1,na,1,na)  <- that must be wrong
    !  remember that na=krn(kr)%size
    call angleT3D(h3,indx,krn(kr)%start,na,krn(kr)%start,na)

    ! rewind(u22)
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

          if(idx2(id1)%size < 1) cycle
          do ic=1,nn(2)
            do ib=1,nn(1)
              if(idx1(ib,ic,id)%size < 1 .OR. idx1(ib,ic,id1)%size < 1) cycle

              ! write(6,'(3i4,i8,3i8))') ib,ic,id1, id, &
              !         idx1(ib,ic,id1)%size,idx1(ib,ic,id)%size, &
              !         idx1(ib,ic,id1)%size*idx1(ib,ic,id)%size

              allocate(vm(ib,ic,id1)%mat(idx1(ib,ic,id1)%size, &
                                         idx1(ib,ic,id)%size))
              do ia=1,idx1(ib,ic,id)%size
                !    read(u11,rec=rec_num1(ia,ib,ic,id)) en2,vec12(:)
                vec12 => eigvec1(rec_num1(ia,ib,ic,id))%vec(:)
              !==============================================================
              ! this selection makes it possible either to include or not
              ! the off-diagonal in R contributions from 1/R^2

              if(oner_flag==0) then

                ! USE this to ignore the off-diag contributions
                ! this flag should also reset oner(i,i) to temp/r(i)/r(i)
                ! in laguerre_dvr.90

                e3 = t3(id1,id)*vec12

              else

                ! USE this to include the off-diag contributions
                ! this flag should retain oner(i,i) in laguerre_dvr.90

                e3 = oner3(id1,id) * matmul(h3,vec12)
                ! same as above using BLAS
                ! call dgemv('N',size(h3,1),size(h3,2),oner3(id1,id),h3, &
                !               size(h3,1),vec12,1,zero,e3,1)
                !  action of angular kinetic operator on the right function
                do i=1,na
                  n= krn(kr)%start+i-1
                  l =indx(n,1)
                  j =indx(n,2)
                  temp= oner3(id1,id)*(j*(j+one)+l*(l+one)+jr*(jr+one))+t3(id1,id)
                  e3(i) = temp*vec12(i) + e3(i)
                end do

              end if
              !==============================================================

                do ia1=1,idx1(ib,ic,id1)%size
                  !    read(u11,rec=rec_num1(ia1,ib,ic,id1)) en1,vec11(:)
                  vec11 => eigvec1(rec_num1(ia1,ib,ic,id1))%vec(:)
                  vm(ib,ic,id1)%mat(ia1,ia)= dot_product(vec11,e3)
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

        ! read(u22) en22,vec22(:)
          en0   =  eigvec2(rec_num2(im,id))%eng
          vec22 => eigvec2(rec_num2(im,id))%vec(:)

        ! diagonal contribution

        h6(n) = en0 + t3(id,id)

        ! off-diagonal contributions

        ! rewind(u21)
        do id1=1,id-1
          if(idx2(id1)%size < 1) cycle
          ! allocate(vec21(idx1(nn(1),nn(2),id1)%end))
          do im1=1,idx2(id1)%size
            n1 = idx2(id1)%start+im1-1

            ! write(6,'(2(2x,3i4))') id,im,n, id1,im1,n1

            ! read(u21) en21,vec21(:)
              vec21 => eigvec2(rec_num2(im1,id1))%vec(:)

            temp = 0.0D0
            do ic=1,nn(2)
              do ib=1,nn(1)
                if(idx1(ib,ic,id)%size < 1 .OR. idx1(ib,ic,id1)%size < 1) cycle

                allocate(e6(idx1(ib,ic,id)%size))

                e6 = matmul(  &
                   vec21(idx1(ib,ic,id1)%start:idx1(ib,ic,id1)%end), &
                   vm(ib,ic,id1)%mat(:,:))
                ! same as above using BLAS
                ! call dgemv('T',idx1(ib,ic,id1)%size,idx1(ib,ic,id)%size,one, &
                !           vm(ib,ic,id1)%mat(:,:),idx1(ib,ic,id1)%size, &
                !           vec21(idx1(ib,ic,id1)%start:idx1(ib,ic,id1)%end),1, &
                !           zero,e6,1)
                temp = temp+ dot_product(e6, &
                  vec22(idx1(ib,ic,id)%start:idx1(ib,ic,id)%end))

                deallocate(e6)
              end do
            end do

            h6(n1) = temp

          end do
        end do

        !  store the composed column of h6 to disk
        write(u16) h6(1:n)

      end do

      !  deallocate auxiliary array of matrices

      if(id > 1) then
        do id1=1,id-1

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

    deallocate(h3,e3,h6)

    nullify(vec22,vec21,vec11,vec12)

    do n=1,istage1
      deallocate(eigvec1(n)%eng,eigvec1(n)%vec)
    end do

    do n=1,istage2
      deallocate(eigvec2(n)%eng,eigvec2(n)%vec)
    end do

    deallocate(eigvec1,eigvec2,rec_num1,rec_num2)

    nlev= idx2(nn(3))%end
    ! nlev and istage2 appear to be the same

    write(6,*) ' allocating and reading in h6...'

    allocate(h6((istage2*(istage2+1))/2),e6(istage2),vectors(nlev,abs(icut3)))

    rewind(u16)
    ! h6(:) = zero
    nf=0
    do n=1,istage2
      read(u16) h6(nf+1:nf+n)
      nf=nf+n
    end do
    write(6,*) ' read', nf,' out of',(istage2*(istage2+1))/2
    call time_stamp(u6)
    write(6,*) ' diagonalising h6...'

    nlev=istage2
    call diag_p(h6,e6,vectors,nlev,nout,icut3,encut3,3,idiag)
    if(idiag /= 0) then
      write(6,*) ' diag_p failed to diagonalise on stage 3'
      write(6,*) ' trace: kr=',kr
      stop
    end if

    if(nout < 1) then
      write(u6,*) ' stage 3: no levels'
      write(6,*) ' trace: kr=',kr
      stop
    end if

    ! the zero energy level must be totally symmetric
     if(jr==0 .AND. p==0 .AND. j_parity==0 .AND. &
         l_parity==0 .AND. jl_parity==0) then
      enzero=e6(1)
      write(u6,'(/," zero energy=",f26.8,/)') e6(1)
     end if

    write(u6,'(/," h6 energies (relative to zero energy):")')

    ! write number of eigenvectors and the matrix size
    write(u31(jr)) nout,nlev

    do icounter=1,nout
      write(u6,'(2x,i4,f24.6)') icounter,e6(icounter)-enzero
      write(u31(jr)) e6(icounter),vectors(:,icounter)
    end do

    deallocate(h6,e6,vectors)

    idx3(kr)%size = nout
    if(kr==0) then
      idx3(kr)%start= 1
      idx3(kr)%end = idx3(kr)%start+idx3(kr)%size-1
      ! job's done, do nothing
    else
      idx3(kr)%start= idx3(kr-1)%end+1
      idx3(kr)%end = idx3(kr)%start+idx3(kr)%size-1
      call mainj()
    end if

    call time_stamp(u6)
  end do kr_loop

  ! deallocate vectors
  deallocate(indx,krn)

  !  end of the preparation of K-eigenvectors

  do kr=0,krmax
    close(u41(kr))
    close(u31(kr))
    close(u21(kr))
    close(u11(kr))
  end do

end do jl_parity_loop
end do l_parity_loop
end do j_parity_loop

end do jp_loop

if(expansion_flag == 0) then
  deallocate(phi,tphi,theta,ttheta,theta1,ttheta1)
else ! if(expansion_flag == 1) then
  deallocate(vex, threej0, tp1, tp2)
  close(u10)
end if

deallocate(q1,q2,q3,t1,t2,t3,oner1,oner2,oner3, idx1,idx2,idx3)
deallocate(testpot)

call time_stamp(u6)

CONTAINS


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


SUBROUTINE input_data()
  USE types
  USE param
  USE potential
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
!  icuts, energy cutoffs, margins
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) icut0, icut1, icut2, icut3, icut4
write(u6,'(1x,5i10)') icut0, icut1, icut2, icut3, icut4
read(u5,*)  encut0, encut1, encut2, encut3, encut4
write(u6,'(2x,5f10.0)')  encut0, encut1, encut2, encut3, encut4
read(u5,*) margin0,margin1,margin2,margin3,margin4
write(u6,'(2x,5f10.0)') margin0,margin1,margin2,margin3,margin4
read(u5,*) imargin0,imargin1,imargin2,imargin3,imargin4
write(u6,'(1x,5i10)') imargin0,imargin1,imargin2,imargin3,imargin4
!   equilibrium/reference configuration
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) qe(1),qe(2),qe(3),qe(4),qe(5),qe(6)
write(u6,'(1x,6f9.5)') qe(1),qe(2),qe(3),qe(4),qe(5),qe(6)
!   flags
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) stage_flag
write(u6,'(1x,i10)') stage_flag
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) oner_flag
write(u6,'(1x,i10)') oner_flag
write(u6,'(80("="))')
! ======================================================================
write(u6,'(" INPUT: processing input...",/)')

call masses()
!   take Radau formulae for masses to test 3D angular + rotation block
!   because masses 1-3 will be the same as in the input and an analytical
!   solution is possible. Also set V=0
write(u6,'(" adapted masses:",/,3f14.8,/)')  muq

if (angular_problem_only_flag==1) then
  write(u6,'(" WARNING: stretching coordinates are frozen to given lengths",/)')
  nn=1
  ! open(u8,file='3Dlevels.txt')
end if

if (optimize_flag >0) then
  write(u6,'(" Basis Optimization is ENABLED",/)')
else
  write(u6,'(" Basis Optimization is DISABLED",/)')
end if

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

!  angular potential expansion in Pjk*Plk*cos(k*phi)
!  valid range:  ne1 <= 2*jmax; ne2 <= 2*lmax;
!  ne3 <= min(2*kmax,ne1,ne2)

if(ne1 > 2*jmax) then
  write(u6, '(" ne1 is unnecessary large. CORRECTED")')
  ne1 = 2*jmax
end if
if(ne2 > 2*lmax) then
  write(u6, '(" ne2 is unnecessary large. CORRECTED")')
  ne2 = 2*lmax
end if
if(ne3 > min(ne1,ne2)) then
  write(u6, '(" ne3 is unnecessary large. CORRECTED")')
  ne3 = min(2*kmax,ne1,ne2)
end if

write(u6, '(/," angular basis size:              ", &
  & " jmax=",i3," lmax=",i3," kmax=",i3,/)') jmax, lmax, kmax
write(u6, '(" angular potential expansion size:", &
  & " ne1 =",i3," ne2 =",i3," ne3 =",i3,/)') ne1, ne2, ne3

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

write(u6,'(" stage #;      icut;  E cutoff  E margin")')

write(u6,'(" stage 0:",i10,2f10.1,i10)') icut0,encut0,margin0,imargin0
write(u6,'(" stage 1:",i10,2f10.1,i10)') icut1,encut1,margin1,imargin1
write(u6,'(" stage 2:",i10,2f10.1,i10)') icut2,encut2,margin2,imargin2
write(u6,'(" stage 3:",i10,2f10.1,i10)') icut3,encut3,margin3,imargin3
write(u6,'(" stage 4:",i10,2f10.1,i10)') icut4,encut4,margin4,imargin4

write(u6,'(/,"   Equilibrium/reference Valence configuration:")')
  write(u6,'("   Re1=",f20.8," AA")')   qe(1)
  write(u6,'("   Re2=",f20.8," AA")')   qe(2)
  write(u6,'("   Re3=",f20.8," AA")')   qe(3)
  write(u6,'("   Th1=",f20.8," rad")')  qe(4)
  write(u6,'("   Th2=",f20.8," rad")')  qe(5)
  write(u6,'("   Phi=",f20.8," rad")')  qe(6)

!call coords(tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6), &
!   qe(1),qe(2),qe(3),qe(4),qe(5),qe(6), mass1, mass2, mass3, mass4, opt)

temp = v(qe(1),qe(2),qe(3),qe(4),qe(5),qe(6))

write(u6,'(/,"   Equilibrium/reference energy:",f16.6)') temp

write(u6,'(/," INPUT: done.",/)')

END SUBROUTINE

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

SUBROUTINE setup_idx1()
  USE types
IMPLICIT NONE
INTEGER(I4B) :: n,i1,i2,i3
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
      write(u6,'(" id  full basis  red. basis   start indx  end indx")')
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
        write(u6,'(i3,4i12)') id, idx1(nn(1),nn(2),id)%end, &
          idx2(id)%size,idx2(id)%start,idx2(id)%end
      end do
END SUBROUTINE


! function select returns 
! .FALSE. if the level needs to be selected and 
! .TRUE. otherwise
! it makes the choice based on 
! i - level number, en - level energy, encut -  enrgy cut off,
! margin and imargin - input parameters
!
FUNCTION select(i,en,encut,margin,imargin)
REAL(DP) :: en,margin,encut
INTEGER(I4B) :: imargin,i
LOGICAL(LGT) :: select

! use the default logic without margins
!
if(en > encut) then
  select= .TRUE.
else
  select = .FALSE.
end if

! more complicated logic with margins
!if(i == 1 .AND. en > margin) then
!  select= .TRUE.
!else if(i > imargin .AND. en > encut) then
!  select= .TRUE.
!else
!  select = .FALSE.
!end if

END FUNCTION

END PROGRAM
