program merge
  
implicit none

INTEGER(4) :: narg,iargc, i
CHARACTER*15 :: arg
CHARACTER*80 :: linebuffer
! REAL(8), ALLOCATABLE :: mux(:), muy(:), muz(:)
INTEGER(4), PARAMETER :: n=10000
REAL(8) :: mux(n), muy(n), muz(n)
CHARACTER*25 :: buffer(n)

!===============================================
! input processing
narg=iargc()

! check command line input
if (narg == 0) then
  write(6,*) "for example:"
  write(6,*) "./merge 0000-0000-J1-J0"
  stop
end if

! get diagonilizer algorithm
if (narg > 0) then
  call getarg(1,arg)
  write(6,*) " arg1=", arg
end if

open(10,file=arg//'.mux',status='old')

do i=1,64
  read(10,'(A80)') linebuffer
  !write(6,'(A80)') linebuffer
enddo

do i=1,n
  read(10,'(A22,es16.4)') buffer(i), mux(i)
enddo

write(6,*) ' done'

close(10)


open(10,file=arg//'.muy',status='old')

do i=1,64
  read(10,'(A80)') linebuffer
  !write(6,'(A80)') linebuffer
enddo

do i=1,n
  read(10,'(A22,es16.4)') buffer(i), muy(i)
enddo

close(10)


open(10,file=arg//'.muz',status='old')

do i=1,64
  read(10,'(A80)') linebuffer
  !write(6,'(A80)') linebuffer
enddo

do i=1,n
  read(10,'(A22,es16.4)') buffer(i), muz(i)
enddo

close(10)

do i=1,n
  write(6,'(A25,4es16.4)') buffer(i), mux(i), muy(i), muz(i), (mux(i)+muy(i)+muz(i))**2
enddo


END program
