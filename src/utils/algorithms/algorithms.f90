module algorithms_mod
  use algorithms_integer_mod
  use algorithms_real_mod
  implicit none

  interface prefix_sum_exclusive
    module procedure :: prefix_sum_exclusive_integer, prefix_sum_exclusive_real
  end interface
  
  interface bubble_sort
    module procedure :: bubble_sort_integer, bubble_sort_real
  end interface

  interface get_unique
    module procedure :: get_unique_integer, get_unique_real
  end interface

  interface findloc_all
    module procedure :: findloc_all_integer, findloc_all_real
  end interface
end module
