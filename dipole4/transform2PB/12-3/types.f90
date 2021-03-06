MODULE types
!  Symbolic names or kind types of 4-,2-,and 1-byte integers:
INTEGER,PARAMETER ::I4B =SELECTED_INT_KIND(9)
INTEGER,PARAMETER ::I2B =SELECTED_INT_KIND(4)
INTEGER,PARAMETER ::I1B =SELECTED_INT_KIND(2)
!  Symbolic names or kind types of single-and double-precision reals:
INTEGER,PARAMETER ::SP =KIND(1.0)
INTEGER,PARAMETER ::DP =KIND(1.0D0)
! Symbolic name or kind type of default logical:
INTEGER,PARAMETER ::LGT =KIND(.true.)
! numbers
REAL(DP),PARAMETER :: zero=0.0_dp
REAL(DP),PARAMETER :: half=0.5_dp
REAL(DP),PARAMETER :: one=1.0_dp
REAL(DP),PARAMETER :: two=2.0_dp
REAL(DP),PARAMETER :: three=3.0_dp
! Frequently used mathematical constants (with precision to spare):
REAL(DP),PARAMETER ::PI=3.141592653589793238462643383279502884197_dp
REAL(DP),PARAMETER ::PIO2=1.57079632679489661923132169163975144209858_dp
! REAL(DP),PARAMETER ::PIO4=PIO2/2.0_dp
REAL(DP),PARAMETER ::TWOPI=6.283185307179586476925286766559005768394_dp
REAL(DP),PARAMETER ::SQRT2=1.41421356237309504880168872420969807856967_dp
REAL(DP),PARAMETER ::SQRT3=1.732050807568877293527446341505872366943_dp
REAL(DP),PARAMETER ::SQRT5=2.236067977499789696409173668731276235441_dp
! REAL(DP),PARAMETER :: one_over_sqrtpi= one/sqrt(PI)
REAL(DP),PARAMETER :: one_over_sqrtpi= 0.56418958354775628695_dp
! REAL(DP),PARAMETER :: one_over_sqrt2pi= one/sqrt(TWOPI)
REAL(DP),PARAMETER :: one_over_sqrt2pi= 0.39894228040143267794_dp
! REAL(DP),PARAMETER ::EULER=0.5772156649015328606065120900824024310422_dp
! mword is an amount of memory (in MB) taken up by REAL(DP)
REAL(DP),PARAMETER ::mword = 8./1024./1024.

!  conversion factor for hbar^2/M/R^2
!  where hbar in CI, M in atomic mass, R in Angstrom
!  it needs to be devided by hbar*c
REAL(DP), PARAMETER :: hbar=1.05457266E-34_dp           ! Plank's constant in J*sec
REAL(DP), PARAMETER :: atomic_mass=1.6605402E-27_dp     ! kg
REAL(DP), PARAMETER :: speed_of_light=29979245800._dp   ! cm/sec
REAL(DP), PARAMETER :: angstrom=1.0E-10_dp              ! metres
! REAL(DP), PARAMETER :: convf=hbar/(2*PI*atomic_mass*speed_of_light*angstrom**2)
REAL(DP), PARAMETER :: convf=33.71526230147670_dp

!  new data types
!
TYPE vectorarray
  REAL(DP), POINTER :: vec(:)
END TYPE
TYPE matrixarray
  REAL(DP), POINTER :: mat(:,:)
END TYPE
TYPE eigenvector
  REAL(DP), POINTER :: vec(:)
  REAL(DP), POINTER :: eng
END TYPE

TYPE triple
  INTEGER(I4B) :: size, start, end
END TYPE

!  file units
INTEGER(I4B), PARAMETER :: u5=5, u6=6, u7=6, u8=6, u10=10, u16=16

END MODULE
