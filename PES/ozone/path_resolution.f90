module path_resolution_mod
  use mpi
  implicit none

contains 

!-------------------------------------------------------------------------------------------------------------------------------------------
! Resolves a path relative to location of executable file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function resolve_relative_exe_path(relative_exe_path) result(path)
    character(*), intent(in) :: relative_exe_path
    character(:), allocatable :: path
    integer :: proc_id, ierr, script_unit, result_unit
    character(512) :: proc_id_str, program_name, command_result
    character(:), allocatable :: script_name, result_name

    call MPI_Comm_Rank(MPI_COMM_WORLD, proc_id, ierr)
    write(proc_id_str, '(I0)') proc_id
    script_name = '.temp_resolve_relative_exe_path_' // trim(proc_id_str) // '.py'
    call get_command_argument(0, program_name)

    open(newunit = script_unit, file = script_name)
    write(script_unit, '(A)') 'import shutil'
    write(script_unit, '(A)') 'from pathlib import Path'
    write(script_unit, '(A)') 'symlink_path = shutil.which("' // trim(program_name) // '")'
    write(script_unit, '(A)') 'target_path = Path(symlink_path).resolve().parent'
    write(script_unit, '(A)') 'relative_exe_path = target_path / "' // relative_exe_path // '"'
    write(script_unit, '(A)') 'print(relative_exe_path)'
    flush(script_unit)

    result_name = '.temp_resolve_relative_exe_path_result_' // trim(proc_id_str)
    call execute_command_line('python3 ' // script_name // ' > ' // result_name)

    open(newunit = result_unit, file = result_name)
    read(result_unit, '(A)') command_result
    path = trim(command_result)

    close(result_unit, status = 'delete')
    close(script_unit, status = 'delete')
  end function
  
end module
