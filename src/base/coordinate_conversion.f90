module coordinate_coversion_mod
  use general_utils
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts APH coordinates to mass-scaled Jacobi coordinates.
! *coord_list*: 3 x n array. Each column is a point in APH coordinates in order: rho (bohr), theta (rad), phi (rad).
! The APH coordinates are rather abstract and difficult to visualize. The the code below for exact equations, but qualitatively 
! rho represents the overall size of the "bounding box" of the system,
! theta (ranges from 0 to pi/2) corresponds to bending motion, where 0 - equilateral triangle, and pi/2 - linear configuration,
! phi (ranges from 0 to 2pi) corresponds to "atom permutation" motion.
! See Fig. 6 in V. Kokoouline and C.H. Greene, Phys. Rev. A 68, 012703 (2003).
! *new_coords_list*: 3 x n array. i-th column is a point in mass-scaled Jacobi coordinates, corresponding
! to the i-th column of *coord_list* in order: s (bohr), S (bohr), Theta (rad).
! s is a mass-scaled length of vector from the 2nd terminal atom to 1st terminal atom of symmetric equilibrium configuration.
! S is a mass-scaled length of vector from the center of mass between the terminal atoms to the central atom.
! Theta is an angle between vectors corresponding to s and S.
! See Eqs. 52a-52c (phi = 2*chi) in R.T. Pack and G.A. Parker, J. Chem. Phys. 87, 3888 (1987).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_aph_to_mass_jacobi(coord_list) result(new_coord_list)
    real(real64), intent(in) :: coord_list(:, :)
    real(real64) :: new_coord_list(size(coord_list, 1), size(coord_list, 2))
    integer :: i

    call assert(size(coord_list, 1) == 3, 'Wrong size of coord_list')
    do i = 1, size(coord_list, 2)
      new_coord_list(1, i) = coord_list(1, i) / sqrt(2d0) * sqrt(1 - sin(coord_list(2, i)) * cos(coord_list(3, i)))
      new_coord_list(2, i) = coord_list(1, i) / sqrt(2d0) * sqrt(1 + sin(coord_list(2, i)) * cos(coord_list(3, i)))
      new_coord_list(3, i) = acos(sin(coord_list(2, i)) * sin(coord_list(3, i)) / sqrt(1 - sin(coord_list(2, i))**2 * cos(coord_list(3, i))**2))
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts mass-scaled Jacobi coordinates to regular Jacobi coordinates.
! *coord_list*: 3 x n array. Each column is a point in mass-scaled Jacobi coordinates in order: s (bohr), S (Bohr), Theta (rad).
! See `convert_aph_to_mass_jacobi` for more details.
! *mass*: the masses of the three atoms in order: terminal1, central, terminal2 (of symmetric equilibrium configuration).
! *new_coords_list*: 3 x n array. i-th column is a point in Jacobi coordinates, corresponding
! to the i-th column of *coord_list* in order: r (bohr), R (bohr), Theta (rad).
! r is a length of vector from the 2nd terminal atom to 1st terminal atom of symmetric equilibrium configuration.
! R is a length of vector from the center of mass between the terminal atoms to the central atom.
! Theta is an angle between vectors corresponding to r and R.
! See Eqs. 4-7 in R.T. Pack and G.A. Parker, J. Chem. Phys. 87, 3888 (1987).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_mass_jacobi_to_jacobi(coord_list, mass) result(new_coord_list)
    real(real64), intent(in) :: coord_list(:, :)
    real(real64), intent(in) :: mass(:)
    real(real64) :: new_coord_list(size(coord_list, 1), size(coord_list, 2))
    integer :: i
    real(real64) :: M, mu, d

    call assert(size(coord_list, 1) == 3, 'Wrong size of coord_list')
    call assert(size(mass) == 3, 'Wrong size of mass')
    M = sum(mass)
    mu = sqrt(product(mass) / M)
    d = sqrt(mass(2) / mu * (1 - mass(2) / M))
    do i = 1, size(coord_list, 2)
      new_coord_list(1, i) = coord_list(1, i) * d
      new_coord_list(2, i) = coord_list(2, i) / d
      new_coord_list(3, i) = coord_list(3, i)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts Jacobi coordinates to Caresian coordinates.
! *coord_list*: 3 x n array. Each column is a point in Jacobi coordinates in order: r (Bohr), R (Bohr), Theta (rad).
! See `convert_mass_jacobi_to_jacobi` for more details.
! *mass*: the masses of the three atoms in order: terminal1, central, terminal2 (of symmetric equilibrium configuration).
! *new_coords_list*: 3 x n array. i-th column is a point in Cartesian coordinates, corresponding to the i-th column of *coord_list*.
! Z-coordinate is redundant for 3 atoms, so (z1, z2, z3) = (0, 0, 0).
! The first atom (terminal1) is fixed in (x1, y1) = (0, 0).
! The third atom (terminal2) has fixed y3 = 0.
! The remaining 3 non-zero coordinates are given in order: x2 (Bohr), y2 (Bohr), x3 (Bohr).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_jacobi_to_cartesian(coord_list, mass) result(new_coord_list)
    real(real64), intent(in) :: coord_list(:, :)
    real(real64), intent(in) :: mass(:)
    real(real64) :: new_coord_list(size(coord_list, 1), size(coord_list, 2))
    integer :: i
    real(real64) :: mass_center_13_x ! x coordinate of center of mass between atoms 1 and 3 (terminal atoms)

    call assert(size(coord_list, 1) == 3, 'Wrong size of coord_list')
    call assert(size(mass) == 3, 'Wrong size of mass')
    do i = 1, size(coord_list, 2)
      mass_center_13_x = mass(3) / (mass(1) + mass(3)) * coord_list(1, i)
      new_coord_list(1, i) = mass_center_13_x - coord_list(2, i) * cos(coord_list(3, i))
      new_coord_list(2, i) = coord_list(2, i) * sin(coord_list(3, i))
      new_coord_list(3, i) = coord_list(1, i)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts Cartesian coordinates to all bonds coordinates.
! *coord_list*: 3 x n array. Each column is a point in Cartesian coordinates in order: x2 (Bohr), y2 (Bohr), x3 (Bohr).
! See `convert_jacobi_to_cartesian` for more details.
! *new_coords_list*: 3 x n array. i-th column is a point in all bonds coordinates, corresponding to the i-th column of *coord_list*
! in order: bond12 (Bohr), bond13 (Bohr), bond23 (Bohr).
! The order of atoms is: terminal1, center, terminal2 in symmetric equilibrium configuration.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_cartesian_to_all_bonds(coord_list) result(new_coord_list)
    real(real64), intent(in) :: coord_list(:, :)
    real(real64) :: new_coord_list(size(coord_list, 1), size(coord_list, 2))
    integer :: i

    call assert(size(coord_list, 1) == 3, 'Wrong size of coord_list')
    do i = 1, size(coord_list, 2)
      new_coord_list(1, i) = sqrt(coord_list(1, i)**2 + coord_list(2, i)**2)
      new_coord_list(2, i) = coord_list(3, i)
      new_coord_list(3, i) = sqrt((coord_list(1, i) - coord_list(3, i))**2 + coord_list(2, i)**2)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts all bonds coordinates to internal coordinates.
! *coord_list*: 3 x n array. Each column is a point in all bonds coordinates in order: bond12 (Bohr), bond13 (Bohr), bond23 (Bohr).
! See `convert_cartesian_to_all_bonds` for more details.
! *new_coords_list*: 3 x n array. i-th column is a point in internal coordinates, corresponding to the i-th column of *coord_list*
! in order: bond12 (Bohr), bond23 (Bohr), angle123 (rad).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_all_bonds_to_internal(coord_list) result(new_coord_list)
    real(real64), intent(in) :: coord_list(:, :)
    real(real64) :: new_coord_list(size(coord_list, 1), size(coord_list, 2))
    integer :: i

    call assert(size(coord_list, 1) == 3, 'Wrong size of coord_list')
    do i = 1, size(coord_list, 2)
      new_coord_list(1, i) = coord_list(1, i)
      new_coord_list(2, i) = coord_list(3, i)
      new_coord_list(3, i) = acos((coord_list(1, i)**2 + coord_list(3, i)**2 - coord_list(2, i)**2) / (2 * coord_list(1, i) * coord_list(3, i))) ! law of cosines
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! A shortcut for APH to Jacobi conversion.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_aph_to_jacobi(coord_list, mass) result(new_coord_list)
    real(real64), intent(in) :: coord_list(:, :)
    real(real64), intent(in) :: mass(:)
    real(real64) :: new_coord_list(size(coord_list, 1), size(coord_list, 2))

    new_coord_list = convert_aph_to_mass_jacobi(coord_list)
    new_coord_list = convert_mass_jacobi_to_jacobi(new_coord_list, mass)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! A shortcut for APH to Cartesian conversion.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_aph_to_cartesian(coord_list, mass) result(new_coord_list)
    real(real64), intent(in) :: coord_list(:, :)
    real(real64), intent(in) :: mass(:)
    real(real64) :: new_coord_list(size(coord_list, 1), size(coord_list, 2))

    new_coord_list = convert_aph_to_jacobi(coord_list, mass)
    new_coord_list = convert_jacobi_to_cartesian(new_coord_list, mass)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! A shortcut for APH to all bonds conversion.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_aph_to_all_bonds(coord_list, mass) result(new_coord_list)
    real(real64), intent(in) :: coord_list(:, :)
    real(real64), intent(in) :: mass(:)
    real(real64) :: new_coord_list(size(coord_list, 1), size(coord_list, 2))

    new_coord_list = convert_aph_to_cartesian(coord_list, mass)
    new_coord_list = convert_cartesian_to_all_bonds(new_coord_list)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! A shortcut for APH to internal conversion.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_aph_to_internal(coord_list, mass) result(new_coord_list)
    real(real64), intent(in) :: coord_list(:, :)
    real(real64), intent(in) :: mass(:)
    real(real64) :: new_coord_list(size(coord_list, 1), size(coord_list, 2))

    new_coord_list = convert_aph_to_all_bonds(coord_list, mass)
    new_coord_list = convert_all_bonds_to_internal(new_coord_list)
  end function

end module
