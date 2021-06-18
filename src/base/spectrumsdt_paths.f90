!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures that provide paths to SpectrumSDT files and folders. 
! Most procedures take either an instance of input parameters or `sym_path`, which is a path to a folder with calculations for a specific
! value of K and symmetry.
!-------------------------------------------------------------------------------------------------------------------------------------------
module spectrumsdt_paths_mod
  use general_utils_mod
  use input_params_mod
  use path_utils_mod
  use rovib_utils_base_mod
  implicit none

  interface get_k_folder_path
    module procedure :: get_k_folder_path_root, get_k_folder_path_root_range, get_k_folder_path_params
  end interface

  interface get_sym_path
    module procedure :: get_sym_path_params
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns name of rho info file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rho_info_name() result(res)
    character(:), allocatable :: res
    res = 'rho_info.txt'
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to rho grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rho_info_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, get_rho_info_name())
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns name of theta info file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_theta_info_name() result(res)
    character(:), allocatable :: res
    res = 'theta_info.txt'
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to theta grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_theta_info_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, get_theta_info_name())
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns name of phi info file.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_phi_info_name() result(res)
    character(:), allocatable :: res
    res = 'phi_info.txt'
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to phi grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_phi_info_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, get_phi_info_name())
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns name of the file with the requested values of PES.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_pes_request_path() result(res)
    character(:), allocatable :: res
    res = 'pes_in.txt'
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with calculated PES.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_pes_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, 'pes_out.txt')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for given Ks. K range version.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_folder_path_root_range(root_path, K, J, parity) result(res)
    character(*), intent(in) :: root_path
    integer, intent(in) :: K(2)
    integer, optional, intent(in) :: J, parity
    character(:), allocatable :: res
    character(:), allocatable :: k_folder_name
    
    if (present(J) .and. present(parity)) then
      if (K(1) == get_k_start(J, parity) .and. K(2) == J) then
        k_folder_name = 'K_all'
      end if
    end if

    if (.not. allocated(k_folder_name)) then
      if (K(1) == K(2)) then
        k_folder_name = 'K_' // num2str(K(1))
      else
        k_folder_name = 'K_' // num2str(K(1)) // '..' // num2str(K(2))
      end if
    end if

    res = append_path_tokens(root_path, k_folder_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for given Ks. Single K version.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_folder_path_root(root_path, K) result(res)
    character(*), intent(in) :: root_path
    integer, intent(in) :: K
    character(:), allocatable :: res
    res = get_k_folder_path_root_range(root_path, [K, K])
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for given Ks, taken from input *params*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_folder_path_params(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res

    if (params % use_rovib_coupling == 0) then
      res = get_k_folder_path_root_range(params % root_path, params % K)
    else
      res = get_k_folder_path_root_range(params % root_path, params % K, params % J, params % parity)
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for a given Ks and symmetry.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_sym_path_int(k_path, sym_code) result(res)
    character(*), intent(in) :: k_path
    integer, intent(in) :: sym_code
    character(:), allocatable :: res
    character(:), allocatable :: sym_name

    if (sym_code == 2) then
      sym_name = 'symmetry_all'
    else
      sym_name = 'symmetry_' // num2str(sym_code)
    end if
    res = append_path_token(k_path, sym_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for a given symmetry and Ks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_sym_path_root(root_path, K, sym_code) result(res)
    character(*), intent(in) :: root_path
    integer, intent(in) :: K, sym_code
    character(:), allocatable :: res
    character(:), allocatable :: k_path

    k_path = get_k_folder_path(root_path, K)
    res = get_sym_path_int(k_path, sym_code)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for symmetry and Ks, taken from input *params*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_sym_path_params(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    character(:), allocatable :: k_path

    k_path = get_k_folder_path(params)
    res = get_sym_path_int(k_path, params % basis % symmetry)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with basis calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_basis_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    res = append_path_tokens(sym_path, 'basis')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with basis results calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_basis_bin_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: basis_path

    basis_path = get_basis_path(sym_path)
    res = append_path_tokens(basis_path, 'bin')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with number of 1D basis functions in each rho-slice.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_basis_1D_summary_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: basis_path

    basis_path = get_basis_path(sym_path)
    res = append_path_tokens(basis_path, 'num_vectors_1d.fwc')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to k-block info (overlap structure).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_block_info_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: basis_path

    basis_path = get_basis_path(sym_path)
    res = append_path_tokens(basis_path, 'num_vectors_2d.fwc')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to a formatted file with 1D eigenvalues from all values of theta and rho.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_energies_1d_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    res = append_path_tokens(get_basis_path(sym_path), 'energies_1d.fwc')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to a formatted file with 2D eigenvalues from all rho slices.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_energies_2d_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    res = append_path_tokens(get_basis_path(sym_path), 'energies_2d.fwc')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to a binary file with 1D eigenvalues and eigenvectors from all theta slices in a specific rho slice.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_solutions_1d_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: basis_bin_path, file_name

    basis_bin_path = get_basis_bin_path(sym_path)
    file_name = 'basis_1d.' // num2str(slice_ind, '(I0)') // '.bin'
    res = append_path_tokens(basis_bin_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to a binary file with 2D eigenvalues and eigenvectors in a specific rho slice.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_solutions_2d_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: basis_bin_path, file_name

    basis_bin_path = get_basis_bin_path(sym_path)
    file_name = 'basis_2d.' // num2str(slice_ind, '(I0)') // '.bin'
    res = append_path_tokens(basis_bin_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with overlaps calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_overlaps_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    res = append_path_tokens(sym_path, 'overlaps')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with overlaps results calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_overlaps_bin_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_path

    overlaps_path = get_overlaps_path(sym_path)
    res = append_path_tokens(overlaps_path, 'bin')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with a regular overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_regular_overlap_path(sym_path, slice_ind_1, slice_ind_2) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind_1, slice_ind_2
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_bin_path, file_name

    overlaps_bin_path = get_overlaps_bin_path(sym_path)
    file_name = 'overlap.' // num2str(slice_ind_1, '(I0)') // '.' // num2str(slice_ind_2, '(I0)') // '.bin'
    res = append_path_tokens(overlaps_bin_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with a symmetric overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_symmetric_overlap_J_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_bin_path, file_name

    overlaps_bin_path = get_overlaps_bin_path(sym_path)
    file_name = 'sym_J.' // num2str(slice_ind) // '.bin'
    res = append_path_tokens(overlaps_bin_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with a symmetric overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_symmetric_overlap_K_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_bin_path, file_name

    overlaps_bin_path = get_overlaps_bin_path(sym_path)
    file_name = 'sym_K.' // num2str(slice_ind) // '.bin'
    res = append_path_tokens(overlaps_bin_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with a coriolis overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_coriolis_overlap_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_bin_path, file_name

    overlaps_bin_path = get_overlaps_bin_path(sym_path)
    file_name = 'coriolis.' // num2str(slice_ind) // '.bin'
    res = append_path_tokens(overlaps_bin_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with an asymmetric overlap block. Negative values of slice_ind indicate K=1 block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_asymmetric_overlap_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_bin_path, file_name

    overlaps_bin_path = get_overlaps_bin_path(sym_path)
    file_name = 'asym.' // num2str(slice_ind) // '.bin'
    res = append_path_tokens(overlaps_bin_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with an asymmetric overlap block. Negative values of slice_ind indicate K=1 block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_asymmetric_overlap_1_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_bin_path, file_name

    overlaps_bin_path = get_overlaps_bin_path(sym_path)
    file_name = 'asym_1.' // num2str(slice_ind) // '.bin'
    res = append_path_tokens(overlaps_bin_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with eigenpairs calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_eigensolve_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_tokens(get_sym_path(params), 'eigensolve')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to current stage.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_stage_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_tokens(get_sym_path(params), params % stage)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Appends parity folder to a given stage folder if parity matters.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_parity_path(params, stage_path) result(res)
    class(input_params), intent(in) :: params
    character(*), optional, intent(in) :: stage_path
    character(:), allocatable :: res
    character(:), allocatable :: stage_path_act

    stage_path_act = arg_or_default(stage_path, get_stage_path(params))
    if (params % use_rovib_coupling == 1 .and. params % K(1) <= 1) then
      res = append_path_tokens(stage_path_act, 'parity_' // num2str(params % parity))
    else
      res = stage_path_act
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with binary stage results.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_stage_binary_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    character(:), allocatable :: parity_path

    parity_path = get_parity_path(params)
    res = append_path_tokens(parity_path, 'bin')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with CAP.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_cap_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    character(:), allocatable :: parity_path

    parity_path = get_parity_path(params, get_eigensolve_path(params))
    res = append_path_tokens(parity_path, 'cap.fwc')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with eigenpairs results calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_eigensolve_bin_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    character(:), allocatable :: parity_path

    parity_path = get_parity_path(params, get_eigensolve_path(params))
    res = append_path_tokens(parity_path, 'bin')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with 3D expansion coefficients for a given state number k.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_solution_3d_path(params, k) result(res)
    class(input_params), intent(in) :: params
    integer, intent(in) :: k ! solution index
    character(:), allocatable :: res
    character(:), allocatable :: eigensolve_bin_folder, file_name

    eigensolve_bin_folder = get_eigensolve_bin_path(params)
    file_name = 'exp.' // num2str(k) // '.bin'
    res = append_path_tokens(eigensolve_bin_folder, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to computed spectrum.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_states_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    character(:), allocatable :: parity_path

    parity_path = get_parity_path(params, get_eigensolve_path(params))
    res = append_path_tokens(parity_path, 'states.fwc') ! fixed width columns
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with properties calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_properties_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_tokens(get_sym_path(params), 'properties')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with state properties calculation.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_state_properties_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    character(:), allocatable :: parity_path

    parity_path = get_parity_path(params, get_properties_path(params))
    res = append_path_tokens(parity_path, 'state_properties.fwc') ! fixed width columns
  end function

end module
