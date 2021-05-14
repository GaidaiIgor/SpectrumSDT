module algorithms_mod
  use algorithms_integer_mod
  use algorithms_real_mod
  implicit none

  interface findloc_all
    module procedure :: findloc_all_integer, findloc_all_real
  end interface

end module
