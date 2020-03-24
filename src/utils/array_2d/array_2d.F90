module array_2d_mod
! aggregates modules for all template types
  use array_2d_complex_mod
  use array_2d_real_mod

  interface flatten_array_2d
    module procedure :: flatten_array_2d_complex, flatten_array_2d_real
  end interface
end module
