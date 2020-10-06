!-------------------------------------------------------------------------------------------------------------------------------------------
! System specific routines
!-------------------------------------------------------------------------------------------------------------------------------------------
module system_mod
  use general_utils
  use io_utils
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Executes specied shell command and returns command's output
! Uses provided file_name to store the results of execution (default is .command_output.txt)
!-------------------------------------------------------------------------------------------------------------------------------------------
  function execute_shell_command(command, file_name) result(output)
    character(*), intent(in) :: command
    character(*), optional, intent(in) :: file_name
    character(:), allocatable :: output
    character(:), allocatable :: command_act, file_name_act

    file_name_act = arg_or_default(file_name, '.command_output.txt')
    command_act = command // ' > ' // file_name_act
    call execute_command_line(command_act)
    output = read_file(file_name_act, delete = 1)
  end function
end module
