!-------------------------------------------------------------------------------------------------------------------------------------------
! Contains historically global variables shared by some other modules
!-------------------------------------------------------------------------------------------------------------------------------------------
module general_vars
  implicit none
  real*8, allocatable :: grho2(:), sintet2(:)
  complex*16, allocatable :: freq1z(:), der1z(:, :)

  ! Masses
  integer :: mol        ! Molecule code, for example, 686
  real*8 :: m0          ! Mass of central atom
  real*8 :: m1          ! Mass of 1st terminal atom
  real*8 :: m2          ! Mass of 2nd terminal atom
  real*8 :: mtot        ! Total mass
  real*8 :: mu          ! Three-particle reduced mass
  real*8 :: mu0         ! Reduced mass of m1 and m2
  real*8 :: mu1         ! Reduced mass of m0 and m2
  real*8 :: mu2         ! Reduced mass of m0 and m1

  ! Grids
  integer :: n1, n2, n3 ! Problem dimensions
  real*8, allocatable :: g1(:), g2(:), g3(:)
  real*8, allocatable :: jac1(:), jac2(:), jac3(:)
  real*8 :: alpha1, alpha2, alpha3

  ! Directories
  character(:), allocatable :: outdir

  ! Network
  integer :: myid        ! Process id
  integer :: nprocs      ! Number of processes
end module
