module path_resolution_mod
  use mpi
  implicit none

contains 

!-------------------------------------------------------------------------------------------------------------------------------------------
! Executes specied shell command and returns command's output. Parallel version.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function execute_shell_command(command) result(output)
    character(*), intent(in) :: command
    character(:), allocatable :: output
    integer :: proc_id, ierr, file_unit
    character(512) :: proc_id_str, command_result
    character(:), allocatable :: file_name, command_full

    ! The only way to get command output is to read it from file...
    call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
    write(proc_id_str, '(I0)') proc_id
    file_name = '.temp' // trim(proc_id_str)
    command_full = command // ' > ' // file_name
    call execute_command_line(command_full)

    ! Read the results
    open(newunit = file_unit, file = file_name)
    read(file_unit, '(A)') command_result
    close(file_unit, status = 'delete')
    output = trim(command_result)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Extracts everything before the last token in path.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_path_head(path) result(head)
    character(*), intent(in) :: path
    character(:), allocatable :: head
    integer :: ind_slash

    ind_slash = index(path, '/', .true.)
    head = path(:ind_slash-1)
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Resolves a path relative to location of executable file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function resolve_relative_exe_path(relative_exe_path) result(path)
    character(*), intent(in) :: relative_exe_path
    character(:), allocatable :: path
    character(512) :: path_arg

    call get_command_argument(0, path_arg)
    path = execute_shell_command('readlink -f $(which ' // trim(path_arg) // ')')
    path = get_path_head(path)
    ! only if path is non-empty
    if (len(path) > 0) then
      path = path // '/'
    end if
    path = path // relative_exe_path
  end function
  
end module
