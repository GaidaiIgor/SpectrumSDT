module debug_tools
  use debug_tools_base
  use general_utils
  use input_params_mod
  implicit none

contains

  subroutine init_debug(params)
    class(input_params), intent(in) :: params
    debug_mode = params % debug_mode
    if (allocated(params % debug_param_1)) then
      debug_int_1 = str2int(params % debug_param_1)
    end if
  end subroutine
end module
