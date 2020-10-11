module debug_tools
  use general_utils
  use input_params_mod
  implicit none

  character(:), allocatable :: debug_mode
  character(:), allocatable :: test_mode
  integer :: debug_int_1
  integer :: parity
  integer :: signal = 0

contains

  subroutine init_debug(params)
    class(input_params), intent(in) :: params
    debug_mode = params % debug_mode
    test_mode = params % test_mode
    debug_int_1 = str2int(params % debug_param_1)
    parity = params % parity
  end subroutine
end module
