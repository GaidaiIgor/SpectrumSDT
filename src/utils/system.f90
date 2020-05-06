!-------------------------------------------------------------------------------------------------------------------------------------------
! System specific routines
!-------------------------------------------------------------------------------------------------------------------------------------------
module system_mod
  use general_utils
  use io_utils

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
    output = read_file(file_name_act, destruct = 1)
  end function

! !-----------------------------------------------------------------------
! ! List directory content
! ! max file name length = 255
! ! mask - file name pattern (like '*.txt')
! ! additional - additional arguments to find
! !-----------------------------------------------------------------------
!   subroutine list_directory(path, mask, result, additional)
!     character(*), intent(in) :: path, mask
!     character(:), allocatable :: command, tempFileName
!     character(255), allocatable :: result(:)
!     character(*), optional :: additional
!     integer :: fileState, fileUnit, linesAmount, i
!
!     tempFileName = '.searchResult'
!     command = 'find ' // path // ' -maxdepth 1 -name "' // mask // '" ! -name "' // tempFileName // '" -type f'
!     if (present(additional)) then
!       command = command // ' ' // additional
!     end if
!     command = command // ' -printf "%f\n" > ' // tempFileName
!     call system(command)
!     linesAmount = countLines(tempFileName)
!     allocate(result(linesAmount))
!
!     open(newunit = fileUnit, file = tempFileName)
!     i = 1
!     do
!       read(fileUnit, *, iostat = fileState) result(i)
!       if (fileState < 0) then
!         exit
!       end if
!       i = i + 1
!     end do
!     close(fileUnit, status = 'delete')
!   end subroutine
end module
