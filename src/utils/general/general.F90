module general_utils
  use iso_fortran_env
  use ifcore

  implicit none

  interface arg_or_default
    module procedure :: arg_or_default_integer, arg_or_default_real, arg_or_default_real_array, arg_or_default_char_str
  end interface
  
  interface swap
    module procedure :: swap_integer, swap_real, swap_char_str
  end interface
  
  interface num2str
    module procedure :: num2str_integer, num2str_real
  end interface

  interface operator (.means.)
    procedure :: logical_consequence
  end interface

  contains
  
#include "general_integer.F90"
#include "general_real.F90"
#include "general_real_array.F90"
#include "general_char_str.F90"

!-----------------------------------------------------------------------
! Allocates array and compies the value of initialization array
!-----------------------------------------------------------------------
  subroutine initAllocatable(array, initializer)
    real*8, allocatable :: array(:)
    real*8 :: initializer(:)

    allocate(array(size(initializer)))
    array = initializer
  end subroutine

!-----------------------------------------------------------------------
! Compares given reals with specified precision
!-----------------------------------------------------------------------
  function compare_reals(a, b, comp_precision) result(res)
    real*8, intent(in) :: a, b
    real*8, optional :: comp_precision
    integer :: res
    real*8 :: difference, precision_act

    precision_act = arg_or_default(comp_precision, 1d-10)
    difference = a - b
    if (abs(difference) < precision_act) then
      res = 0
    else if (difference > 0) then
      res = 1
    else if (difference < 0) then
      res = -1
    end if
  end function

!-----------------------------------------------------------------------
! generates real grid using specified step
!-----------------------------------------------------------------------
  function generateRealRange(start, end, step) result(grid)
    real*8, intent(in) :: start, end, step
    real*8, allocatable :: grid(:)
    real*8 :: npointsf
    integer :: npoints, i

    npointsf = (end - start) / step
    npoints = npointsf
    ! If the last step is complete up to comparison accuracy
    if (compare_reals(npointsf - npoints, 1.d0) == 0) then
      npoints = npoints + 1
    end if
    
    ! +1 because of an extra point at the beginning of interval
    allocate(grid(npoints + 1))
    grid = [(start + i * step, i = 0,npoints)]
  end function
  
!-----------------------------------------------------------------------
! Generates equally spaced grid of points in the interval [start, end]
!-----------------------------------------------------------------------
  function linspace(start, end, npoints) result(grid)
    real*8, intent(in) :: start, end
    integer, intent(in) :: npoints
    real*8 :: grid(npoints)
    integer :: i
    real*8 :: step

    step = (end - start) / (npoints - 1)
    grid = [(start + i * step, i = 0,npoints - 1)]
  end function

!---------------------------------------------------------------------------------------------------------------------------------------------
! Prints progress %
! progress: current progress of some process (a number from 0 to 1)
! last_progress: the value of progress for the previous call. If it exceeds the value of progress by more than progress_step, progress is reported.
! progress_step: controls how often progress should be reported (default is every 0.1)
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine track_progress(progress, progress_step, reset)
    real*8, value :: progress
    real*8, optional :: progress_step
    integer, optional, intent(in) :: reset
    real*8, save :: last_progress
    real*8 :: progress_step_act

    if (present(reset)) then
      if (reset == 1) then
        last_progress = 0
      end if
    end if

    progress_step_act = merge(progress_step, 0.1d0, present(progress_step))
    if (progress > last_progress + progress_step_act) then
      last_progress = last_progress + progress_step_act
      print '(A,x,F6.2,A$)', '\r', last_progress * 100, '% done'
      if (compare_reals(last_progress, 1d0) == 0) then
        print *
      end if
    end if
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! A one-liner to check various constraints
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine assert(logical_statemenet, error_message)
    logical, value :: logical_statemenet
    character(*) :: error_message

    if (.not. logical_statemenet) then
      call TRACEBACKQQ(error_message)
      stop
    end if
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts integer to string
!-------------------------------------------------------------------------------------------------------------------------------------------
  function num2str_integer(num, format) result(str)
    integer, intent(in) :: num
    character(*), intent(in), optional :: format
    character(:), allocatable :: str
    character(256) :: buffer
    character(:), allocatable :: format_act
    
    format_act = arg_or_default(format, '(I0)')
    write(buffer, format_act) num
    str = trim(buffer)
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts real*8 to string
!-------------------------------------------------------------------------------------------------------------------------------------------
  function num2str_real(num, format) result(str)
    real*8, intent(in) :: num
    character(*), intent(in), optional :: format
    character(:), allocatable :: str
    character(256) :: buffer
    character(:), allocatable :: format_act
    
    format_act = arg_or_default(format, '(G10.2)')
    write(buffer, format_act) num
    str = trim(adjustl(buffer))
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to integer
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function str2int(str) result(int)
    character(*), intent(in) :: str
    integer :: int
    read(str, *) int
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts string to real*8
!-------------------------------------------------------------------------------------------------------------------------------------------
  elemental function str2real(str) result(real_num)
    character(*), intent(in) :: str
    real*8 :: real_num
    read(str, *) real_num
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Checks if a given strings starts with a given pattern
!-------------------------------------------------------------------------------------------------------------------------------------------
  function str_starts_with(str, pattern) result(res)
    character(*), intent(in) :: str, pattern
    logical :: res
    res = str(1:len(pattern)) == pattern
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Removes characters in *char_set* from the right end of *str* until first occurence of a character not from *char_set*
!-------------------------------------------------------------------------------------------------------------------------------------------
  function str_trim(str, char_set) result(res)
    character(*), intent(in) :: str
    character(*), optional, intent(in) :: char_set
    character(:), allocatable :: res
    integer :: pos
    character(:), allocatable :: char_set_act

    char_set_act = arg_or_default(char_set, ' ')
    pos = verify(str, char_set_act, .true.)
    res = str(1:pos)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns identity matrix of specified size
!-------------------------------------------------------------------------------------------------------------------------------------------
  function identity_matrix(size) result(matrix)
    integer, intent(in) :: size
    integer, allocatable :: matrix(:, :)
    integer :: i

    allocate(matrix(size, size))
    matrix = 0
    do i = 1, size
      matrix(i, i) = 1
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Tertiary if operator
!-------------------------------------------------------------------------------------------------------------------------------------------
  function iff(true_res, false_res, cond) result(res)
    character(*), intent(in) :: true_res, false_res
    logical, intent(in) :: cond
    character(:), allocatable :: res

    if (cond) then
      res = true_res
    else
      res = false_res
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Compares leading digits of given integers
!-------------------------------------------------------------------------------------------------------------------------------------------
  function compare_leading_digits(int1, int2) result(res)
    integer, intent(in) :: int1, int2
    logical :: res
    character(:), allocatable :: int1_str, int2_str

    int1_str = num2str(int1)
    int2_str = num2str(int2)
    res = str_starts_with(int1_str, int2_str)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Evaluates logical consequence function
!-------------------------------------------------------------------------------------------------------------------------------------------
  function logical_consequence(a, b) result(res)
    logical, intent(in) :: a, b
    logical :: res
    res = .not. (a .and. .not. b)
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Evaluates delta Kronecker function
!-------------------------------------------------------------------------------------------------------------------------------------------
  function delta(a, b) result(res)
    integer, intent(in) :: a, b
    integer :: res
    res = merge(1, 0, a == b)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Find _row_ and _col_ index of maximum element in the _matrix_
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine find_max_loc_matrix(matrix, row, col)
    real*8, intent(in) :: matrix(:, :)
    integer, intent(out) :: row, col
    integer, allocatable :: max_loc_col(:) ! row of maximum in each column
    real*8, allocatable :: max_val_col(:) ! value of maximum in each column
    
    max_loc_col = maxloc(matrix, 2)
    max_val_col = maxval(matrix, 2)
    col = maxloc(max_val_col, 1)
    row = max_loc_col(col)
  end subroutine
end module
