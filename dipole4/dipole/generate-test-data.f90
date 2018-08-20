! 15/05/2009 INK
!
! generate fake eigen-functions

PROGRAM generate
  USE types
  USE param
IMPLICIT NONE

INTEGER(I4B), PARAMETER :: u11=11

REAL(DP) :: zpe

INTEGER(I4B), ALLOCATABLE :: indx(:,:)
TYPE(triple), ALLOCATABLE :: krn(:)

REAL(DP), ALLOCATABLE :: vec(:)

!  cut off related
  INTEGER(I4B) :: icut0, icut1, icut2, icut3, icut4

INTEGER(I4B) :: ia,ib,ic,id,i,n
INTEGER(I4B) :: jr, kr

INTEGER(I4B) :: jp, j_parity, l_parity, jl_parity
INTEGER(I4B) :: nstates, na, icounter, iacc

CHARACTER(12) :: fname
CHARACTER(3) :: smt
CHARACTER(LEN=80) title


! ====================================================================
! INPUT - WAVR4 part
!
!open(u6,file='out.txt')

open(u5,file='input.txt',status='old')

call input_data()

close(u5)

! just set these for backward compatibility with WAVR4
qe(1)=re(1); qe(2)=re(2); qe(3)=re(3); qe(4)=zero; qe(5)=zero; qe(6)=zero

!  find max number of angular states (namax) and mmax
!
call angular_states_max()

write(6,*) " namax=", namax

! final state angular basis indexing
allocate(indx(namax,6),krn(0:krmax))

smt='all'

! ====================================================================

jp_loop:         do jp=0,1
j_parity_loop:   do j_parity =0,j_parity_max
l_parity_loop:   do l_parity =0,l_parity_max
jl_parity_loop:  do jl_parity =0,jl_parity_max

jr_loop: do jr=0,jrmax

  fname = 'x'//i2c(jp)//i2c(j_parity)//i2c(l_parity)//i2c(jl_parity)//'-J'//i2c(jr)//'-'//smt
  open(unit=u11,file=fname//".dat",access='SEQUENTIAL',form='UNFORMATTED')
  write(6,*) fname

  !call angular_states_qnumbers(indx,krn)
  call angular_states_qnumbers(indx,krn, jr,jp,j_parity,l_parity,jl_parity)

  ! compute the size of the primitive basis
  nstates=0
  do kr=0,jr
      nstates = nstates+ krn(kr)%size*nn(1)*nn(2)*nn(3)
  enddo
  write(6,*) " full primitive basis size=", nstates
  ! obviously we cant have more states than the basis
  nstates = min(nstates,icut4)

  write(u11) nstates

  allocate(vec(1:nstates))
  
  do i=1,nstates
    vec(i) = dble(i)
  enddo
  write(u11) vec(1:nstates)

  deallocate(vec)

  !!write(6,*) " #  l  j kp kr  k  m ib ic id"

  iacc = 0

  kr_loop: do kr=0,jr
    !  angular sub-block
    na  = krn(kr)%size
    allocate(vec(1:na))
    do id=1,nn(3); do ic=1,nn(2); do ib=1,nn(1)
      do icounter=1,nstates

        vec(1:na) = zero
        do ia=1,na
          ! write 1 just once for each wavefunction
          i=iacc+ia
          if(i==icounter) then
            !! l1 =indx1(n1,1)
            !! j1 =indx1(n1,2)
            !! kp1=indx1(n1,3)
            !! kr1=indx1(n1,4)
            !! k1 =indx1(n1,5)
            !! m1 =indx1(n1,6)
            !!n =ia+krn(kr)%start-1
            !!write(6,'(10i3)') icounter,indx(n,1),indx(n,2),indx(n,3),indx(n,4),indx(n,5),indx(n,6),ib,ic,id
            vec(ia)= one
          endif
        end do
        write(u11) vec(1:na)

      enddo
      iacc = iacc + na
    enddo; enddo; enddo
    deallocate(vec)
  enddo kr_loop

  close(u11)

end do jr_loop

end do jl_parity_loop
end do l_parity_loop
end do j_parity_loop
end do jp_loop


deallocate(indx,krn)


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
!  icuts, energy cutoffs, margins
read(u5,'(A80)') title; write(u6,'(A80)') title
read(u5,*) icut0, icut1, icut2, icut3, icut4
write(u6,'(1x,5i10)') icut0, icut1, icut2, icut3, icut4

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

END PROGRAM
