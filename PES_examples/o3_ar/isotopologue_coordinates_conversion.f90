module isotopologue_coordinates_conversion_mod
  use iso_fortran_env, only: real64
  use lapack_interface_mod
  implicit none

  ! Cartesian coordinates of ozone-666 atoms in a.u.
  ! Coordinates of each atom are given in columns in order x, y, z
  ! The coordinates are chosen to make 1-2 and 2-3 distances equal to 1.2717 Ang, 1-2-3 angle to 116.84 deg, as specified in the corresponding paper by Dawes and coworkers
  ! The center of mass is at the origin and axes coincide with the prinical axes of inertia
  ! x-axis coincides with the principle axis corresponding to the smallest moment of inertia
  ! z-axis coincides with the principle axis corresponding to the middle moment of inertia
  ! y-axis coincides with the principle axis corresponding to the largest moment of inertia
  real(real64), private, parameter :: ozone_cartesian(3, 3) = reshape([ 2.047279540451686d0, 0d0, -0.419503293481828d0, &
                                                                        0d0,                 0d0,  0.839006586963654d0, &
                                                                       -2.047279540451686d0, 0d0, -0.419503293481828d0], [3, 3])
  ! Masses of oxygen isotopes (16, 17, 18) in a.u.
  real(real64), private, parameter :: oxygen_masses(3) = [2.915694567487589d4, 3.098752142701935d4, 3.281046079389339d4]

  interface operator (-)
    module procedure :: matrix_minus_vector
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Subtracts a given vector from all columns of a given matrix.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function matrix_minus_vector(matrix, vector) result(res)
    real(real64), intent(in) :: matrix(:, :)
    real(real64), intent(in) :: vector(:)
    real(real64), allocatable :: res(:, :)
    integer :: j

    res = matrix
    do j = 1, size(res, 2)
      res(:, j) = res(:, j) - vector
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts isotopomer label to an array of atom masses in a.u.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_atom_masses(isotopomer) result(masses)
    character(3), intent(in) :: isotopomer
    real(real64) :: masses(3)
    integer :: i, atom_ind

    do i = 1, 3
      read(isotopomer(i:i), *) atom_ind
      atom_ind = atom_ind - 5
      masses(i) = oxygen_masses(atom_ind)
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates center of mass for a given matrix of cartesian coordinates (3xN, order: x, y, z) and corresponding particle masses (N).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_mass_center(coords_cartesian, masses) result(mass_center)
    real(real64), intent(in) :: coords_cartesian(:, :)
    real(real64), intent(in) :: masses(:)
    real(real64), allocatable :: mass_center(:)
    real(real64), allocatable :: weights(:)

    weights = masses / sum(masses)
    mass_center = matmul(coords_cartesian, weights)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates tensor of inertia.
! coords_cartesian: matrix of cartesian coordinates, 3xN, each point is written in column in order: x, y, z.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_inertia_tensor(coords_cartesian, masses) result(inertia_tensor)
    real(real64), intent(in) :: coords_cartesian(:, :)
    real(real64), intent(in) :: masses(:)
    real(real64) :: inertia_tensor(3, 3)

    inertia_tensor(1, 1) = dot_product(masses, coords_cartesian(2, :) ** 2 + coords_cartesian(3, :) ** 2)
    inertia_tensor(2, 2) = dot_product(masses, coords_cartesian(1, :) ** 2 + coords_cartesian(3, :) ** 2)
    inertia_tensor(3, 3) = dot_product(masses, coords_cartesian(1, :) ** 2 + coords_cartesian(2, :) ** 2)
    inertia_tensor(1, 2) = -dot_product(masses, coords_cartesian(1, :) * coords_cartesian(2, :))
    inertia_tensor(1, 3) = -dot_product(masses, coords_cartesian(1, :) * coords_cartesian(3, :))
    inertia_tensor(2, 3) = -dot_product(masses, coords_cartesian(2, :) * coords_cartesian(3, :))
    inertia_tensor(2, 1) = inertia_tensor(1, 2)
    inertia_tensor(3, 1) = inertia_tensor(1, 3)
    inertia_tensor(3, 2) = inertia_tensor(2, 3)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Rearranges order and direction of eigenvectors of inertia tensor to match the order (x, y, z) and direction of the reference isotopomer (666).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function rearrange_inertia_tensor(inertia_tensor) result(rearranged_tensor)
    real(real64), intent(in) :: inertia_tensor(3, 3)
    real(real64) :: rearranged_tensor(3, 3)
    integer :: j

    rearranged_tensor = inertia_tensor
    ! Rearrange axes into xyz order
    rearranged_tensor(:, 2) = inertia_tensor(:, 3)
    rearranged_tensor(:, 3) = inertia_tensor(:, 2)
    do j = 1, size(rearranged_tensor, 2)
      if (rearranged_tensor(j, j) < 0) then
        ! Enforce positive axes direction
        rearranged_tensor(:, j) = -rearranged_tensor(:, j)
      end if
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transforms spherical coordinates to cartesian.
! Spherical coordinates are given in order: R, theta (rad), phi (rad).
! Cartesian coordinates are given in order: x, y, z.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function spherical_to_cartesian(spherical) result(cartesian)
    real(real64), intent(in) :: spherical(3)
    real(real64) :: cartesian(3)

    cartesian(1) = spherical(1) * sin(spherical(2)) * cos(spherical(3))
    cartesian(2) = spherical(1) * sin(spherical(2)) * sin(spherical(3))
    cartesian(3) = spherical(1) * cos(spherical(2))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transforms cartesina coordinates to spherical.
! Cartesian coordinates are given in order: x, y, z.
! Spherical coordinates are given in order: R, theta (rad), phi (rad).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function cartesian_to_spherical(cartesian) result(spherical)
    real(real64), intent(in) :: cartesian(3)
    real(real64) :: spherical(3)

    spherical(1) = sqrt(sum(cartesian ** 2))
    spherical(2) = acos(cartesian(3) / spherical(1))
    spherical(3) = atan2(cartesian(2), cartesian(1))
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates rotation matrix and shift for a given isotopomer.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine calculate_isotopomer_shift_rotation(isotopomer, shift, rotation)
    character(3), intent(in) :: isotopomer
    real(real64), intent(out) :: shift(3)
    real(real64), intent(out) :: rotation(3, 3)
    real(real64) :: masses(3), inertia_moments(3)
    real(real64) :: isotopomer_cartesian(3, 3), inertia_tensor(3, 3)

    masses = get_atom_masses(isotopomer)
    shift = get_mass_center(ozone_cartesian, masses)
    isotopomer_cartesian = ozone_cartesian - shift
    inertia_tensor = get_inertia_tensor(isotopomer_cartesian, masses)
    call lapack_eigensolver_real_plain(inertia_tensor, inertia_moments)
    inertia_tensor = rearrange_inertia_tensor(inertia_tensor)
    rotation = transpose(lapack_matrix_inverse(inertia_tensor))
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Fast version of conversion that can reuse the same rotation matrix and shift for multiple points.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_coordinates_to_standard_ozone_fast(shift, rotation, coords_spherical) result(standard_coords_spherical)
    real(real64), intent(in) :: shift(3), coords_spherical(3)
    real(real64), intent(in) :: rotation(3, 3)
    real(real64) :: standard_coords_spherical(3)
    real(real64) :: coords_cartesian(3), standard_coords_cartesian(3)

    coords_cartesian = spherical_to_cartesian(coords_spherical)
    standard_coords_cartesian = matmul(rotation, coords_cartesian)
    standard_coords_cartesian = standard_coords_cartesian + shift
    standard_coords_spherical = cartesian_to_spherical(standard_coords_cartesian)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts spherical coordinates of a 3rd body atom (e.g. argon) given w.r.t. center of mass and principal axes of inertia of a given ozone isotopomer
! to corresponding coordinates in the frame of reference of the standard ozone isotopomer (666).
! The isotopomer in question is specified with a standard 3-letter shorthand notation (e.g. 668).
! The spherical coordinate are given in order: R (a.u.), theta (rad), phi(rad).
! R is the distance from the center of mass of the specified isotopomer.
! theta is the angle of deviation from the z-axis.
! phi is the angle of deviation from the x-axis.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function convert_coordinates_to_standard_ozone(isotopomer, coords_spherical) result(standard_coords_spherical)
    character(3), intent(in) :: isotopomer
    real(real64), intent(in) :: coords_spherical(3)
    real(real64) :: standard_coords_spherical(3)
    real(real64) :: shift(3)
    real(real64) :: rotation(3, 3)

    call calculate_isotopomer_shift_rotation(isotopomer, shift, rotation)
    standard_coords_spherical = convert_coordinates_to_standard_ozone_fast(shift, rotation, coords_spherical)
  end function

end module

program test
  use isotopologue_coordinates_conversion_mod
  real(real64) :: standard_coords(3)
  real(real64), parameter :: pi = acos(-1d0)

  standard_coords = convert_coordinates_to_standard_ozone('668', [5d0, pi/2d0, pi/2d0])
  print *, standard_coords
end program
