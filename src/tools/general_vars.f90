!-------------------------------------------------------------------------------------------------------------------------------------------
! Contains historically global variables shared by some other modules
!-------------------------------------------------------------------------------------------------------------------------------------------
module general_vars
  use iso_fortran_env, only: real64
  implicit none
  real(real64), allocatable :: grho2(:), sintet2(:)
  complex(real64), allocatable :: freq1z(:), der1z(:, :)

  ! Masses
  integer :: mol        ! Molecule code, for example, 686
  real(real64) :: m0          ! Mass of central atom
  real(real64) :: m1          ! Mass of 1st terminal atom
  real(real64) :: m2          ! Mass of 2nd terminal atom
  real(real64) :: mtot        ! Total mass
  real(real64) :: mu          ! Three-particle reduced mass
  real(real64) :: mu0         ! Reduced mass of m1 and m2
  real(real64) :: mu1         ! Reduced mass of m0 and m2
  real(real64) :: mu2         ! Reduced mass of m0 and m1

  ! Grids
  integer :: n1, n2, n3 ! Problem dimensions
  real(real64), allocatable :: g1(:), g2(:), g3(:)
  real(real64), allocatable :: jac1(:), jac2(:), jac3(:)
  real(real64) :: alpha1, alpha2, alpha3

  ! Directories
  character(:), allocatable :: outdir
end module
