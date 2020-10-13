!-------------------------------------------------------------------------------------------------------------------------------------------
! This module contains functions that operate with system paths
!-------------------------------------------------------------------------------------------------------------------------------------------
module path_utils
  use general_utils
  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Creates specified folder
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine create_path(path)
    character(*), intent(in) :: path
    call execute_command_line('mkdir -p ' // path)
  end subroutine

!-----------------------------------------------------------------------
! extracts last token from path
!-----------------------------------------------------------------------
  function get_path_tail(path)
    character(*), intent(in) :: path
    character(:), allocatable :: get_path_tail
    integer :: ind_slash

    ind_slash = index(path, '/', .true.)
    get_path_tail = path(ind_slash+1:)
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
  
!-----------------------------------------------------------------------
! removes specified number of tokens from the end of path
!-----------------------------------------------------------------------
  function strip_path_tokens(path, n_tokens)
    character(*), intent(in) :: path
    integer, intent(in) :: n_tokens
    character(:), allocatable :: strip_path_tokens
    integer :: i
    
    strip_path_tokens = path
    do i = 1,n_tokens
      strip_path_tokens = get_path_head(strip_path_tokens)
    end do
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Appends a new token to a path
!-------------------------------------------------------------------------------------------------------------------------------------------
  function append_path_token(path, token) result(new_path)
    character(*), intent(in) :: path, token
    character(:), allocatable :: new_path, separator
    
    separator = '/'
    new_path = path
    if (new_path(len(new_path) - len(separator) + 1 : len(new_path)) /= separator) then
      new_path = new_path // separator
    end if
    new_path = new_path // token
  end function
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Appends a new token to a path
!-------------------------------------------------------------------------------------------------------------------------------------------
  function append_path_tokens(path, token1, token2, token3, token4) result(new_path)
    character(*), intent(in) :: path, token1
    character(*), optional, intent(in) :: token2, token3, token4
    character(:), allocatable :: new_path
    
    new_path = append_path_token(path, token1)
    if (present(token2)) then
      new_path = append_path_token(new_path, token2)
    end if
    if (present(token3)) then
      new_path = append_path_token(new_path, token3)
    end if
    if (present(token4)) then
      new_path = append_path_token(new_path, token4)
    end if
  end function

end module
