module debug_tools_mod
  use debug_tools_base_mod
  use input_params_mod
  use io_utils_mod
  implicit none

contains

  subroutine init_debug(params)
    class(input_params), intent(in) :: params
    debug_mode = params % debug % mode
    debug_int = params % debug % int_param
    debug_real = params % debug % real_param
  end subroutine
end module
