!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures that provide paths to SpectrumSDT files and folders. 
! Most procedures take either an instance of input parameters or `sym_path`, which is a path to a folder with calculations for a specific
! value of K and symmetry.
!-------------------------------------------------------------------------------------------------------------------------------------------
module spectrumsdt_paths_mod
  use general_utils
  use input_params_mod
  use path_utils
  use rovib_utils_base_mod
  implicit none

  interface get_k_folder_path
    module procedure :: get_k_folder_path_root, get_k_folder_path_root_range, get_k_folder_path_params
  end interface

  interface get_sym_path
    module procedure :: get_sym_path_str, get_sym_path_params
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to rho grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_rho_info_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, 'grid_rho.dat')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to theta grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_theta_info_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, 'grid_theta.dat')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to phi grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_phi_info_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, 'grid_phi.dat')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with calculated PES.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_pes_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, 'pes.out')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for given Ks. K range version.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_folder_path_root_range(root_path, K, J, parity, stage, rovib_coupling) result(res)
    character(*), intent(in) :: root_path
    integer, intent(in) :: K(2)
    integer, optional, intent(in) :: J, parity
    character(*), optional, intent(in) :: stage
    integer, optional, intent(in) :: rovib_coupling
    integer :: J_act, parity_act, rovib_coupling_act
    character(:), allocatable :: res
    character(:), allocatable :: stage_act, k_folder_name
    
    J_act = arg_or_default(J, -1)
    parity_act = arg_or_default(parity, -1)
    stage_act = arg_or_default(stage, '')
    rovib_coupling_act = arg_or_default(rovib_coupling, -1)
    if (rovib_coupling_act == 1 .and. (stage_act == 'eigensolve' .or. stage_act == 'properties') .and. &
        J_act /= -1 .and. parity_act /= -1 .and. K(1) == get_k_start(J_act, parity_act) .and. K(2) == J_act) then
      k_folder_name = 'K_all'
    else if (K(1) == K(2)) then
      k_folder_name = 'K_' // num2str(K(1))
    else
      k_folder_name = 'K_' // num2str(K(1)) // '..' // num2str(K(2))
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
    res = get_k_folder_path_root_range(params % root_path, params % K, params % J, params % parity, params % stage, params % use_rovib_coupling)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for a given symmetry and Ks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_sym_path_str(k_path, sym_name, parity) result(res)
    character(*), intent(in) :: k_path, sym_name
    integer, intent(in), optional :: parity
    character(:), allocatable :: res

    res = k_path
    if (present(parity)) then
      if (parity /= -1) then
        res = append_path_tokens(res, 'parity_' // num2str(parity)) 
      end if
    end if
    res = append_path_tokens(res, sym_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for a given symmetry and Ks.
! Parity is relevant for coupled eigensolves only.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_sym_path_int(k_path, sym_code, parity) result(res)
    character(*), intent(in) :: k_path
    integer, intent(in) :: sym_code
    integer, intent(in), optional :: parity
    character(:), allocatable :: res
    character(:), allocatable :: sym_name

    sym_name = 'symmetry_' // num2str(sym_code)
    res = get_sym_path_str(k_path, sym_name, parity)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for a given symmetry and Ks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_sym_path_root(root_path, K, sym_code, parity) result(res)
    character(*), intent(in) :: root_path
    integer, intent(in) :: K, sym_code
    integer, intent(in), optional :: parity
    character(:), allocatable :: res
    character(:), allocatable :: k_path

    k_path = get_k_folder_path(root_path, K)
    res = get_sym_path_int(k_path, sym_code, parity)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for symmetry and Ks, taken from input *params*.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_sym_path_params(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    character(:), allocatable :: k_path

    k_path = get_k_folder_path(params)
    if (params % use_rovib_coupling == 1 .and. (params % stage == 'eigensolve' .or. params % stage == 'properties')) then
      res = get_sym_path_int(k_path, params % basis % symmetry, params % parity)
    else
      res = get_sym_path_int(k_path, params % basis % symmetry)
    end if
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
  function get_basis_results_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: basis_path

    basis_path = get_basis_path(sym_path)
    res = append_path_tokens(basis_path, 'out_basis')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with number of 1D basis functions in each rho-slice.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_basis_1D_summary_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: basis_path

    basis_path = get_basis_path(sym_path)
    res = append_path_tokens(basis_path, 'num_vectors_1d.dat')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to k-block info (overlap structure).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_block_info_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: basis_path

    basis_path = get_basis_path(sym_path)
    res = append_path_tokens(basis_path, 'num_vectors_2d.dat')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to a formatted file with 1D eigenvalues from all values of theta and rho.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_energies_1d_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    res = append_path_tokens(get_basis_path(sym_path), 'energies_1d.txt')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to a formatted file with 2D eigenvalues from all rho slices.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_energies_2d_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    res = append_path_tokens(get_basis_path(sym_path), 'energies_2d.txt')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to a binary file with 1D eigenvalues and eigenvectors from all theta slices in a specific rho slice.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_solutions_1d_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: basis_results_path, file_name

    basis_results_path = get_basis_results_path(sym_path)
    file_name = 'bas1.' // num2str(slice_ind, '(I0)') // '.bin.out'
    res = append_path_tokens(basis_results_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to a binary file with 2D eigenvalues and eigenvectors in a specific rho slice.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_solutions_2d_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: basis_results_path, file_name

    basis_results_path = get_basis_results_path(sym_path)
    file_name = 'bas2.' // num2str(slice_ind, '(I0)') // '.bin.out'
    res = append_path_tokens(basis_results_path, file_name)
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
  function get_overlaps_results_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_path

    overlaps_path = get_overlaps_path(sym_path)
    res = append_path_tokens(overlaps_path, 'out_overlaps')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with a regular overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_regular_overlap_file_path(sym_path, slice_ind_1, slice_ind_2) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind_1, slice_ind_2
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_results_path, file_name

    overlaps_results_path = get_overlaps_results_path(sym_path)
    file_name = 'overlap.' // num2str(slice_ind_1, '(I0)') // '.' // num2str(slice_ind_2, '(I0)') // '.bin.out'
    res = append_path_tokens(overlaps_results_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with a symmetric overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_symmetric_overlap_J_file_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_results_path, file_name

    overlaps_results_path = get_overlaps_results_path(sym_path)
    file_name = 'sym_J.' // num2str(slice_ind) // '.bin.out'
    res = append_path_tokens(overlaps_results_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with a symmetric overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_symmetric_overlap_K_file_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_results_path, file_name

    overlaps_results_path = get_overlaps_results_path(sym_path)
    file_name = 'sym_K.' // num2str(slice_ind) // '.bin.out'
    res = append_path_tokens(overlaps_results_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with a coriolis overlap block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_coriolis_overlap_file_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_results_path, file_name

    overlaps_results_path = get_overlaps_results_path(sym_path)
    file_name = 'coriolis.' // num2str(slice_ind) // '.bin.out'
    res = append_path_tokens(overlaps_results_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with an asymmetric overlap block. Negative values of slice_ind indicate K=1 block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_asymmetric_overlap_file_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_results_path, file_name

    overlaps_results_path = get_overlaps_results_path(sym_path)
    file_name = 'asym.' // num2str(slice_ind) // '.bin.out'
    res = append_path_tokens(overlaps_results_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with an asymmetric overlap block. Negative values of slice_ind indicate K=1 block.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_asymmetric_overlap_file_1_path(sym_path, slice_ind) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: slice_ind
    character(:), allocatable :: res
    character(:), allocatable :: overlaps_results_path, file_name

    overlaps_results_path = get_overlaps_results_path(sym_path)
    file_name = 'asym_1.' // num2str(slice_ind) // '.bin.out'
    res = append_path_tokens(overlaps_results_path, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with eigenpairs calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_eigensolve_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    res = append_path_tokens(sym_path, 'eigensolve')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with eigenpairs results calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_eigensolve_results_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: eigensolve_path

    eigensolve_path = get_eigensolve_path(sym_path)
    res = append_path_tokens(eigensolve_path, 'out_eigensolve')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with 3D expansion coefficients for a given state number k.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_solution_3d_path(sym_path, k) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: k ! solution index
    character(:), allocatable :: res
    character(:), allocatable :: eigensolve_results_folder, file_name

    eigensolve_results_folder = get_eigensolve_results_path(sym_path)
    file_name = 'exp.' // num2str(k, '(I0)') // '.bin.out'
    res = append_path_tokens(eigensolve_results_folder, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to computed spectrum.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_spectrum_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: eigensolve_path

    eigensolve_path = get_eigensolve_path(sym_path)
    res = append_path_tokens(eigensolve_path, 'states.fwc') ! Fixed width columns
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with properties calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_properties_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    res = append_path_tokens(sym_path, 'properties')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with properties results calculation.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_properties_result_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: properties_path

    properties_path = get_properties_path(sym_path)
    res = append_path_tokens(properties_path, 'states.ssdtp') ! SpectrumSDT properties
  end function

end module
