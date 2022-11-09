module constants_mod
  use iso_fortran_env, only: real64
  implicit none

  real(real64), parameter :: pi = acos(-1d0)
  real(real64), parameter :: rad_to_deg = 180 / pi
  ! The below constants are taken from: Thomas J. Bruno, Paris D. Svoronos CRC Handbook of Chemistry and Physics - 93rd edition (2012), page 1-13
  ! https://books.google.com/books?id=-BzP7Rkl7WkC&pg=SA8-PA31&lpg=SA8-PA31&dq=Thomas+J.+Bruno,+Paris+D.+Svoronos+CRC+Handbook+of+Chemistry+and+Physics+-+93rd+edition&source=bl&ots=kpIQrBNJuv&sig=MVecdbBZUQZGTiL4M5439oFC48Q&hl=en&sa=X&ved=0ahUKEwit2tOfxvnbAhVk2IMKHRbtCowQ6AEISzAE#v=onepage&q=Thomas%20J.%20Bruno%2C%20Paris%20D.%20Svoronos%20CRC%20Handbook%20of%20Chemistry%20and%20Physics%20-%2093rd%20edition&f=false
  real(real64), parameter :: au_to_wn = 219474.6313708d0 ! cm^-1 / Eh (wavenumbers per Hartree)
  real(real64), parameter :: amu_to_kg = 1.660538921d-27 ! kg / amu (kilograms per atomic mass unit)
  real(real64), parameter :: aum_to_kg = 9.10938291d-31 ! kg / aum (kilogram per electron (atomic unit of mass))
  real(real64), parameter :: amu_to_aum = amu_to_kg / aum_to_kg ! aum / amu
  real(real64), parameter :: hydrogen_masses(3) = [1.00782503207d0, 2.0141017778d0, 3.0160492777d0] * amu_to_aum ! H1-3
  real(real64), parameter :: helium_masses(2) = [3.0160293191d0, 4.00260325415d0] * amu_to_aum ! He3-4
  real(real64), parameter :: lithium_masses(2) = [6.015122795d0, 7.016003427d0] * amu_to_aum ! Li6-7
  real(real64), parameter :: beryllium_masses(1) = [9.01218305d0] * amu_to_aum ! Be9
  real(real64), parameter :: boron_masses(2) = [10.0129370d0, 11.0093054d0] * amu_to_aum ! B10-11
  real(real64), parameter :: carbon_masses(4) = [11.0114336d0, 12.0000000d0, 13.0033548378d0, 14.003241989d0] * amu_to_aum ! C11-14
  real(real64), parameter :: nitrogen_masses(2) = [14.0030740048d0, 15.0001088982d0] * amu_to_aum ! N14-15
  real(real64), parameter :: oxygen_masses(3) = [15.99491461956d0, 16.99913170d0, 17.9991596129d0] * amu_to_aum ! O16-18
  real(real64), parameter :: fluorine_masses(2) = [18.0009380d0, 18.99840322d0] * amu_to_aum ! F18-19
  real(real64), parameter :: neon_masses(3) = [19.9924401754d0, 20.99384668d0, 21.991385114d0] * amu_to_aum ! Ne20-22
  real(real64), parameter :: sodium_masses(3) = [21.9944364d0, 22.9897692809d0, 23.99096278d0] * amu_to_aum ! Na22-24
  real(real64), parameter :: magnesium_masses(3) = [23.985041700d0, 24.98583692d0, 25.982592929d0] * amu_to_aum ! Mg24-26
  real(real64), parameter :: aluminium_masses(1) = [26.98153863d0] * amu_to_aum ! Al27
  real(real64), parameter :: silicon_masses(3) = [27.9769265325d0, 28.976494700d0, 29.97377017d0] * amu_to_aum ! Si28-30
  real(real64), parameter :: phosphorus_masses(2) = [30.97376163d0, 31.97390727d0] * amu_to_aum ! P31-32
  real(real64), parameter :: sulfur_masses(5) = [31.97207100d0, 32.97145876d0, 33.96786690d0, 34.96903216d0, 35.96708076d0] * amu_to_aum ! S32-36
  real(real64), parameter :: chlorine_masses(2) = [34.96885268d0, 36.96590259d0] * amu_to_aum ! Cl35,37
  real(real64), parameter :: argon_masses(3) = [35.967545106d0, 37.9627324d0, 39.9623831225d0] * amu_to_aum ! Ar36,38,40
  
contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns mass of the atom (in a.u.) specified by its name followed by mass number. Returns 0 if atom_name is not recognized.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_atom_mass(atom_name) result(mass)
    character(*), intent(in) :: atom_name
    real(real64) :: mass

    select case (atom_name)
      case ('H1')
        mass = hydrogen_masses(1)
      case ('H2')
        mass = hydrogen_masses(2)
      case ('H3')
        mass = hydrogen_masses(3)
      case ('He3')
        mass = helium_masses(1)
      case ('He4')
        mass = helium_masses(2)
      case ('Li6')
        mass = lithium_masses(1)
      case ('Li7')
        mass = lithium_masses(2)
      case ('Be9')
        mass = beryllium_masses(1)
      case ('B10')
        mass = boron_masses(1)
      case ('B11')
        mass = boron_masses(2)
      case ('C11')
        mass = carbon_masses(1)
      case ('C12')
        mass = carbon_masses(2)
      case ('C13')
        mass = carbon_masses(3)
      case ('C14')
        mass = carbon_masses(4)
      case ('N14')
        mass = nitrogen_masses(1)
      case ('N15')
        mass = nitrogen_masses(2)
      case ('O16')
        mass = oxygen_masses(1)
      case ('O17')
        mass = oxygen_masses(2)
      case ('O18')
        mass = oxygen_masses(3)
      case ('F18')
        mass = fluorine_masses(1)
      case ('F19')
        mass = fluorine_masses(2)
      case ('Ne20')
        mass = neon_masses(1)
      case ('Ne21')
        mass = neon_masses(2)
      case ('Ne22')
        mass = neon_masses(3)
      case ('Na22')
        mass = sodium_masses(1)
      case ('Na23')
        mass = sodium_masses(2)
      case ('Na24')
        mass = sodium_masses(3)
      case ('Mg24')
        mass = magnesium_masses(1)
      case ('Mg25')
        mass = magnesium_masses(2)
      case ('Mg26')
        mass = magnesium_masses(3)
      case ('Al27')
        mass = aluminium_masses(1)
      case ('Si28')
        mass = silicon_masses(1)
      case ('Si29')
        mass = silicon_masses(2)
      case ('Si30')
        mass = silicon_masses(3)
      case ('P31')
        mass = phosphorus_masses(1)
      case ('P32')
        mass = phosphorus_masses(2)
      case ('S32')
        mass = sulfur_masses(1)
      case ('S33')
        mass = sulfur_masses(2)
      case ('S34')
        mass = sulfur_masses(3)
      case ('S35')
        mass = sulfur_masses(4)
      case ('S36')
        mass = sulfur_masses(5)
      case ('Cl35')
        mass = chlorine_masses(1)
      case ('Cl37')
        mass = chlorine_masses(2)
      case ('Ar36')
        mass = argon_masses(1)
      case ('Ar38')
        mass = argon_masses(2)
      case ('Ar40')
        mass = argon_masses(3)
      case default
        mass = 0
    end select
  end function

end module
