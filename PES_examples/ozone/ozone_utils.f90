module ozone_utils_mod
  use coordinate_coversion_mod
  use iso_fortran_env, only: real64
  use pes_utils_mod
  implicit none

  real(real64), parameter :: ozone_isotope_masses(3) = [15.99491461956d0, 16.99913170d0, 17.9991596129d0] * amu_to_aum ! aum; masses of isotopes of oxygen 16, 17 and 18 respectively
  ! As computed by program zpeO2
  real(real64), parameter :: zpe_66 = 791.6382691754641d0 / au_to_wn ! J = 0
  real(real64), parameter :: zpe_67 = 779.9050607081324d0 / au_to_wn ! J = 0
  real(real64), parameter :: zpe_68 = 769.3708543295361d0 / au_to_wn ! J = 0
  real(real64), parameter :: zpe_77 = 767.9904668774815d0 / au_to_wn ! J = 0
  real(real64), parameter :: zpe_78 = 757.2885872707559d0 / au_to_wn ! J = 0
  real(real64), parameter :: zpe_88 = 746.4315071358510d0 / au_to_wn ! J = 0

  ! Obtained as energy of dissociated O3 on ozone PES of Dawes for r1 = 900 Bohr, angle = 90 deg (does not matter which one), r2 is optimized to 10^-15 precision (see dawes_de.f90)
  real(real64), parameter :: De = 9274.99560025014d0 / au_to_wn

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Initializes masses and shifts. Isotopic composition of ozone has to be supplied via command line argument to determine appropriate shift.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine init_parameters(mass, shift)
    real(real64), intent(out) :: mass(3)
    real(real64), intent(out) :: shift
    integer :: i, atom_code
    character(3) :: ozone_isotopes

    ! Read command line arguments
    call get_command_argument(1, ozone_isotopes)
    if (len(trim(ozone_isotopes)) == 0) then
      stop 'Isotopic composition has to be specified'
    end if

    do i = 1, 3
      read(ozone_isotopes(i:i), *) atom_code
      mass(i) = ozone_isotope_masses(atom_code - 5)
    end do

    ! select corresponding ZPE
    if (ozone_isotopes == '666') then
      shift = -zpe_66
    else if (ozone_isotopes == '676') then
      shift = -zpe_67
    else if (ozone_isotopes == '686') then
      shift = -zpe_68
    else if (ozone_isotopes == '767' .or. ozone_isotopes == '777') then
      shift = -zpe_77
    else if (ozone_isotopes == '787') then
      shift = -zpe_78
    else if (ozone_isotopes == '868' .or. ozone_isotopes == '878' .or. ozone_isotopes == '888') then
      shift = -zpe_88
    else
      stop 'Unsupported isotopic composition'
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates ozone potential at given *APH coordinates*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function calc_potential_point(aph, mass, shift) result(potential)
    real(real64), intent(in) :: aph(3)
    real(real64), optional, intent(in) :: mass(3)
    real(real64), optional, intent(in) :: shift
    real(real64) :: potential
    real(real64) :: internal(3)
    real(real64) :: all_bonds(3, 1)
    external :: IMLS ! external PES procedure

    all_bonds = convert_aph_to_all_bonds(reshape(aph, [3, 1]), mass)
    ! convert all bonds to internal coordinates expected by IMLS
    internal(1) = all_bonds(1, 1)
    internal(2) = all_bonds(2, 1)
    internal(3) = acos((all_bonds(1, 1)**2 + all_bonds(2, 1)**2 - all_bonds(3, 1)**2) / (2 * all_bonds(1, 1) * all_bonds(2, 1))) * 180 / pi ! 866 angle

    call IMLS(internal, potential, 1)
    potential = potential / au_to_wn - De + shift
  end function

end module
