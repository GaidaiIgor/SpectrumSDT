module constants
  use iso_fortran_env, only: real64
  implicit none

  real(real64), parameter :: pi = acos(-1d0)
  ! The below constants are taken from: Thomas J. Bruno, Paris D. Svoronos CRC Handbook of Chemistry and Physics - 93rd edition (2012), page 1-13
  ! https://books.google.com/books?id=-BzP7Rkl7WkC&pg=SA8-PA31&lpg=SA8-PA31&dq=Thomas+J.+Bruno,+Paris+D.+Svoronos+CRC+Handbook+of+Chemistry+and+Physics+-+93rd+edition&source=bl&ots=kpIQrBNJuv&sig=MVecdbBZUQZGTiL4M5439oFC48Q&hl=en&sa=X&ved=0ahUKEwit2tOfxvnbAhVk2IMKHRbtCowQ6AEISzAE#v=onepage&q=Thomas%20J.%20Bruno%2C%20Paris%20D.%20Svoronos%20CRC%20Handbook%20of%20Chemistry%20and%20Physics%20-%2093rd%20edition&f=false
  real(real64), parameter :: amu = 1.660538921d-27 ! kg / amu (kilograms per atomic mass unit)
  real(real64), parameter :: aum = 9.10938291d-31 ! kg / aum (kilogram per electron (atomic unit of mass))
  real(real64), parameter :: amutoau = amu / aum ! aum / amu
  real(real64), parameter :: isomass(3) = [15.99491461956d0, 16.99913170d0, 17.9991596129d0] * amutoau ! aum; masses of isotopes of oxygen 16, 17 and 18 respectively
  real(real64), parameter :: autown = 219474.6313708d0 ! cm^-1 / Eh (wavenumbers per Hartree)
end module
