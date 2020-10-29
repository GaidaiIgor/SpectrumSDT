!-------------------------------------------------------------------------------------------------------------------------------------------
! Procedures that provide paths to SpectrumSDT files and folders. 
! Most procedures take either an instance of input parameters or `sym_path`, which is a path to a folder with calculations for a specific
! value of K and symmetry.
!-------------------------------------------------------------------------------------------------------------------------------------------
module spectrumsdt_paths_mod
  use general_utils
  use input_params_mod
  use path_utils
  implicit none

  interface get_k_folder_path
    module procedure :: get_k_folder_path_root, get_k_folder_path_params
  end interface

  interface get_sym_path
    module procedure :: get_sym_path_str, get_sym_path_params
  end interface

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to rho grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_grid_rho_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, 'grid_rho.dat')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to theta grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_grid_theta_path(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    res = append_path_token(params % grid_path, 'grid_theta.dat')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to phi grid.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_grid_phi_path(params) result(res)
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
! Generates path to folder with calculation results for given Ks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_folder_path_root(root_path, K) result(res)
    character(*), intent(in) :: root_path
    integer, intent(in) :: K
    character(:), allocatable :: res
    character(:), allocatable :: k_folder_name
    
    k_folder_name = iff(K == -1, 'K_all', 'K_' // num2str(K))
    res = append_path_tokens(root_path, k_folder_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with calculation results for given Ks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_k_folder_path_params(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    integer :: K

    K = merge(-1, params % K(1), params % use_rovib_coupling == 1 .and. params % stage /= 'overlaps')
    res = get_k_folder_path_root(params % root_path, K)
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
! Parity is relevant for coupled eigencalcs only.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_sym_path_int(k_path, sym_code, parity) result(res)
    character(*), intent(in) :: k_path
    integer, intent(in) :: sym_code
    integer, intent(in), optional :: parity
    character(:), allocatable :: res
    character(:), allocatable :: sym_name

    sym_name = iff(sym_code == 0, 'even', 'odd')
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
! Generates path to folder with calculation results for a given symmetry and Ks.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_sym_path_params(params) result(res)
    class(input_params), intent(in) :: params
    character(:), allocatable :: res
    integer :: parity
    character(:), allocatable :: k_path

    k_path = get_k_folder_path_params(params)
    parity = merge(params % parity, -1, params % use_rovib_coupling == 1 .and. params % stage /= 'overlaps')
    res = get_sym_path_int(k_path, params % symmetry, parity)
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
! Generates path to k-block info (overlap structure).
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_block_info_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: basis_results_path

    basis_results_path = get_basis_results_path(sym_path)
    res = append_path_tokens(basis_results_path, 'nvec2.dat')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to the file with 2D energies from all slices.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_2d_energies_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: basis_results_path

    basis_results_path = get_basis_results_path(sym_path)
    res = append_path_tokens(basis_results_path, 'val2.out')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with 1D eigenvalues and eigenvectors from all theta slices in a specific rho slice.
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
! Generates path to file with 2D eigenvalues and eigenvectors in a specific rho slice.
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
  function get_eigencalc_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    res = append_path_tokens(sym_path, 'eigencalc')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to folder with eigenpairs results calculations.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_eigencalc_results_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: eigencalc_path

    eigencalc_path = get_eigencalc_path(sym_path)
    res = append_path_tokens(eigencalc_path, 'out_eigencalc')
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to file with 3D expansion coefficients for a given state number k.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_solution_3d_path(sym_path, k) result(res)
    character(*), intent(in) :: sym_path
    integer, intent(in) :: k ! solution index
    character(:), allocatable :: res
    character(:), allocatable :: eigencalc_results_folder, file_name

    eigencalc_results_folder = get_eigencalc_results_path(sym_path)
    file_name = 'exp.' // num2str(k, '(I0)') // '.bin.out'
    res = append_path_tokens(eigencalc_results_folder, file_name)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Generates path to computed spectrum.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_spectrum_path(sym_path) result(res)
    character(*), intent(in) :: sym_path
    character(:), allocatable :: res
    character(:), allocatable :: eigencalc_path

    eigencalc_path = get_eigencalc_path(sym_path)
    res = append_path_tokens(eigencalc_path, 'states.fwc')
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
    res = append_path_tokens(properties_path, 'states.fwc') ! Fixed width columns
  end function

end module
