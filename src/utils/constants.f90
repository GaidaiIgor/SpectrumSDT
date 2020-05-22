module constants
  real*8, parameter :: pi = acos(-1d0), &
  ! The below constants are taken from: Thomas J. Bruno, Paris D. Svoronos CRC Handbook of Chemistry and Physics - 93rd edition (2012), page 1-13, https://books.google.com/books?id=-BzP7Rkl7WkC&pg=SA8-PA31&lpg=SA8-PA31&dq=Thomas+J.+Bruno,+Paris+D.+Svoronos+CRC+Handbook+of+Chemistry+and+Physics+-+93rd+edition&source=bl&ots=kpIQrBNJuv&sig=MVecdbBZUQZGTiL4M5439oFC48Q&hl=en&sa=X&ved=0ahUKEwit2tOfxvnbAhVk2IMKHRbtCowQ6AEISzAE#v=onepage&q=Thomas%20J.%20Bruno%2C%20Paris%20D.%20Svoronos%20CRC%20Handbook%20of%20Chemistry%20and%20Physics%20-%2093rd%20edition&f=false
                       amu = 1.660538921d-27, & ! kg (atomic mass unit)
                       aum = 9.10938291d-31, & ! kg (mass of electron)
                       amutoau = amu / aum, & ! AMU to atomic unit of mass (electron)
                       autown = 219474.6313708d0, & ! Hartree to cm-1
                       autoev = 27.21138505d0, & ! Hartree to eV
                       isomass(3) = [15.99491461956d0, 16.99913170d0, 17.9991596129d0] * amutoau, & ! masses of isotopes of oxygen 16, 17 and 18 respectively
                       ang2bohr = 0.52917721092d0, & ! Angstroms to Bohr
  ! Zero point energies of O2 are taken as a result of program zpeO2
                       zpe_66 = 791.633780122239d0 / autown, & ! J = 0
                                  ! 794.507661944641
                       ! zpe_66 = 794.5121300074147d0 / autown, & ! J = 1
                       zpe_68 = 769.3708063017872d0 / autown, &
                       zpe_88 = 746.427773959237d0 / autown, & ! J = 0
                       ! zpe_88 = 748.9861634170745d0 / autown, & ! J = 1
  ! Obtained as energy of dissociated O3 on ozone PES of Dawes for r1 = 900 Bohr, angle = 90 deg (does not matter which one), r2 is optimized to 10^-15 precision (see dawes_de.f90)
                       De_dawes = 9274.99560025014d0 / autown
end module
