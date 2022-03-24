program evaluate_pes
  use iso_fortran_env, only: real64
  implicit none
  
  real(real64), parameter :: pi = acos(-1d0)
  real(real64) :: potential ! in cm^-1
  real(real64) :: coords(3)
  character(:), allocatable :: data_path, pes_input_path_1, pes_data_path_1, pes_input_path_2, pes_data_path_2

  data_path = 'dawes/data'
  pes_input_path_1 = data_path // '/input-AUTOSURF-PES.dat'
  pes_data_path_1 = data_path // '/PES-O3-Ar-2112'
  pes_input_path_2 = data_path // '/input-AUTOSURF-PES2.dat'
  pes_data_path_2 = data_path // '/PES-O3-Ar-600'

  coords(1) = 3.295d0                      ! R: in Angstroms
  coords(2) = cos(110.15d0 * pi / 180)     ! theta: in cos(theta)
  coords(3) = 90                           ! phi: in degrees

  call PES_LR(coords, potential, pes_input_path_1, pes_data_path_1, pes_input_path_2, pes_data_path_2)
  print *, coords, potential
end program

