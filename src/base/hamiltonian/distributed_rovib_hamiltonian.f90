!-------------------------------------------------------------------------------------------------------------------------------------------
! Distributed Hamiltonian class.
!-------------------------------------------------------------------------------------------------------------------------------------------
module distributed_rovib_hamiltonian_mod
  use block_borders_mod
  use general_utils_mod
  use formulas_mod
  use input_params_mod
  use iso_fortran_env, only: real64
  use k_block_info_mod
  use matrix_block_info_mod
  use parallel_utils_mod
  use rovib_utils_mod
  use spectrumsdt_io_mod
  use spectrumsdt_paths_mod
  implicit none

  private
  public :: distributed_rovib_hamiltonian, load_overlap_block

  type :: distributed_rovib_hamiltonian
    type(matrix_block_info), pointer :: local_chunk_info => null() ! contains structure info on this chunk without regard for other chunks
    type(matrix_block_info), pointer :: global_chunk_info => null() ! contains structure info on thus chunk with regard for other chunks (location, borders, etc.)
    complex(real64), allocatable :: proc_chunk(:, :)
    integer, allocatable :: all_counts(:) ! how many rows each processor has
    integer, allocatable :: all_shifts(:) ! prefix sum of all_counts
    integer :: compression

  contains

    private
    final :: finalize_distributed_rovib_hamiltonian
    procedure :: load_chunk_info
    procedure :: load_k_block_overlaps
    procedure :: load_chunk_overlaps
    procedure :: include_kinetic_energy
    procedure :: load_k_block_eivals_2d
    procedure :: load_chunk_eivals_2d
    procedure :: include_cap
    procedure :: load_k_block_sym_term
    procedure :: load_chunk_sym_term
    procedure :: load_k_block_coriolis_term
    procedure :: load_chunk_coriolis_term
    procedure :: load_k_block_asym_term
    procedure :: load_chunk_asym_term
    procedure, public :: build
  end type
  
contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! Destructor.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine finalize_distributed_rovib_hamiltonian(this)
    type(distributed_rovib_hamiltonian), intent(inout) :: this
    call deallocate_recursive(this % local_chunk_info)
    call deallocate_recursive(this % global_chunk_info)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads K block sizes on the main diagonal of hamiltonian (based on number of solutions).
! Offdiagonal blocks are ignored. Blocks' positions and borders are ignored, only sizes are initialized here.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_block_sizes_diag(params, ham_info)
    class(input_params), intent(in) :: params
    type(matrix_block_info), intent(out) :: ham_info
    integer :: K, K_load, K_ind, next_sym, total_rows, n_k_blocks
    character(:), allocatable :: root_path, sym_folder, block_info_path
    type(matrix_block_info), pointer :: k_blocks_info(:, :)

    n_k_blocks = params % K(2) - params % K(1) + 1
    allocate(k_blocks_info(n_k_blocks, n_k_blocks))
    next_sym = get_k_symmetry(params % K(1), params % basis % symmetry)
    total_rows = 0
    K_ind = 1

    do K = params % K(1), params % K(2)
      K_load = merge(params % basis % fixed % K, K, params % basis % fixed % enabled == 1)
      root_path = iff(params % basis % fixed % enabled == 1, params % basis % fixed % root_path, params % root_path)
      sym_folder = get_sym_path_root(root_path, K_load, next_sym)
      block_info_path = get_block_info_path(sym_folder)
      call load_k_subblock_sizes_diag(block_info_path, k_blocks_info(K_ind, K_ind))
      total_rows = total_rows + k_blocks_info(K_ind, K_ind) % rows
      next_sym = 1 - next_sym
      K_ind = K_ind + 1
    end do

    ham_info = matrix_block_info(1, 1, total_rows, total_rows)
    ham_info % subblocks => k_blocks_info
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads blocks information for the whole rovib-enabled hamiltonian matrix.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_rovib_hamiltonian_info(params, ham_info)
    class(input_params), intent(in) :: params
    type(matrix_block_info), intent(out) :: ham_info ! enveloping block for all Ks

    call load_block_sizes_diag(params, ham_info) ! Loads sizes of diagonal blocks (down to n-blocks)
    call ham_info % compute_offdiag_subblock_sizes_diag() ! Calculates sizes of offdiagonal blocks. All block sizes are known after this.
    call ham_info % shift_to(1, 1) ! Recursively establishes blocks' borders
    call ham_info % update_block_indexes() ! Recursively establishes blocks' indexes (all blocks know their indexes within the parental block)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transforms K-block info to correspond to compressed offdiagonal K-block representation.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine compress_offdiag_K_block(global_K_block_info, compressed_info)
    type(matrix_block_info), intent(in) :: global_K_block_info
    type(matrix_block_info), intent(out) :: compressed_info
    integer :: n_row_ind, global_n_row_ind

    call compressed_info % copy_without_subblocks(global_K_block_info)
    allocate(compressed_info % subblocks(size(global_K_block_info % subblocks, 1), 1))

    do n_row_ind = 1, size(global_K_block_info % subblocks, 1)
      global_n_row_ind = global_K_block_info % subblocks(n_row_ind, 1) % block_row_ind
      call compressed_info % subblocks(n_row_ind, 1) % copy(global_K_block_info % subblocks(n_row_ind, global_n_row_ind))
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transforms chunk info to correspond to compressed matrix representation.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine compress_chunk_info(global_chunk_info)
    type(matrix_block_info), pointer, intent(inout) :: global_chunk_info
    type(matrix_block_info), pointer :: compressed_info
    integer :: K_row_ind, K_col_ind, global_K_row_ind, global_K_col_ind_left, global_K_col_ind_right, first_row, first_col, rows, cols

    allocate(compressed_info)
    call compressed_info % copy_without_subblocks(global_chunk_info)
    allocate(compressed_info % subblocks(size(global_chunk_info % subblocks, 1), 5))
    ! Copy necessary blocks
    do K_row_ind = 1, size(global_chunk_info % subblocks, 1)
      global_K_row_ind = global_chunk_info % subblocks(K_row_ind, 1) % block_row_ind
      global_K_col_ind_left = max(global_K_row_ind - 2, 1)
      global_K_col_ind_right = min(global_K_row_ind + 2, size(global_chunk_info % subblocks, 2))

      ! Copy non-empty blocks
      do K_col_ind = 1, global_K_col_ind_right - global_K_col_ind_left + 1
        ! Compress offdiagonal blocks
        if (global_K_col_ind_left + K_col_ind - 1 /= global_K_row_ind) then
          call compress_offdiag_K_block(global_chunk_info % subblocks(K_row_ind, global_K_col_ind_left + K_col_ind - 1), compressed_info % subblocks(K_row_ind, K_col_ind))
        else
          call compressed_info % subblocks(K_row_ind, K_col_ind) % copy(global_chunk_info % subblocks(K_row_ind, global_K_col_ind_left + K_col_ind - 1))
        end if
      end do

      ! Add empty blocks for rows where number of blocks is < 5
      do K_col_ind = K_col_ind, 5
        first_row = compressed_info % subblocks(K_row_ind, K_col_ind - 1) % borders % top
        first_col = compressed_info % subblocks(K_row_ind, K_col_ind - 1) % borders % right + 1
        rows = compressed_info % subblocks(K_row_ind, K_col_ind - 1) % borders % bottom - first_row + 1
        cols = 0
        compressed_info % subblocks(K_row_ind, K_col_ind) = matrix_block_info(first_row, first_col, rows, cols)
        ! The empty block does not correspond to an actual block, so its indexes are invalidated to stress that
        compressed_info % subblocks(K_row_ind, K_col_ind) % block_row_ind = -1
        compressed_info % subblocks(K_row_ind, K_col_ind) % block_col_ind = -1
      end do
    end do

    ! Point to the compressed object
    call deallocate_recursive(global_chunk_info)
    global_chunk_info => compressed_info
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Determines part of the overall hamiltonian available to the current process (chunk).
! Saves two chunks: global chunk retains all indexes and borders exactly same as in the global info.
! Local chunk recalculates borders and indexes so that they correspond to the actual hamiltonian chunk stored in the current process.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_chunk_info(this, ham_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    type(matrix_block_info), intent(in) :: ham_info
    integer :: first_row, proc_rows

    call get_proc_elem_range(ham_info % rows, first_row, proc_rows, this % all_counts, this % all_shifts) ! determines first_row and proc_rows values
    allocate(this % global_chunk_info)
    call this % global_chunk_info % copy(ham_info)
    call this % global_chunk_info % cut_rows_between(first_row, first_row + proc_rows - 1)

    if (this % compression == 1) then
      call compress_chunk_info(this % global_chunk_info)
    end if

    allocate(this % local_chunk_info)
    call this % local_chunk_info % copy(this % global_chunk_info)
    if (this % compression == 1) then
      call this % local_chunk_info % calculate_compressed_columns()
    end if
    call this % local_chunk_info % shift_to(1, 1) ! update block borders to correspond to local chunk
    call this % local_chunk_info % update_block_indexes()
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! General wrapper for loading overlap blocks. Converts real blocks to complex.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function load_overlap_block(root_path, K_row, K_col, overlap_type, K_row_sym, slice_ind_row, slice_ind_col, rows, columns, is_file_real) result(block)
    character(*), intent(in) :: root_path
    integer, intent(in) :: K_row, K_col, overlap_type, K_row_sym, slice_ind_row, slice_ind_col, rows, columns
    logical, intent(in) :: is_file_real
    complex(real64), allocatable :: block(:, :)

    if (is_file_real) then
      block = load_overlap_block_real(root_path, K_row, K_col, overlap_type, K_row_sym, slice_ind_row, slice_ind_col, rows, columns)
    else
      block = load_overlap_block_complex(root_path, K_row, K_col, overlap_type, K_row_sym, slice_ind_row, slice_ind_col, rows, columns)
    end if
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads regular overlaps for a diagonal K-block from disk. The block is described by the given block infos.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_k_block_overlaps(this, params, local_k_block_info, global_k_block_info, full_ham_k_block_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), target, intent(in) :: local_k_block_info, global_k_block_info, full_ham_k_block_info
    logical :: is_file_real
    integer :: ir, ic, K, K_load, K_ind, K_sym, slice_ind_row, slice_ind_col, rows, columns, first_row, last_row
    complex(real64), allocatable :: overlap_block(:, :)
    complex(real64), pointer :: local_overlap_info_data(:, :)
    character(:), allocatable :: root_path
    type(matrix_block_info), pointer :: global_overlap_info, full_ham_overlap_info

    K_ind = global_k_block_info % block_row_ind
    K = k_ind_to_k(K_ind, params % K(1))
    K_sym = get_k_symmetry(K, params % basis % symmetry)
    is_file_real = params % use_geometric_phase == 0

    do ir = 1, size(local_k_block_info % subblocks, 1)
      do ic = 1, size(local_k_block_info % subblocks, 2)
        global_overlap_info => global_k_block_info % subblocks(ir, ic)
        slice_ind_row = global_overlap_info % block_row_ind
        slice_ind_col = global_overlap_info % block_col_ind
        full_ham_overlap_info => full_ham_k_block_info % subblocks(slice_ind_row, slice_ind_col)

        if (full_ham_overlap_info % is_empty()) then
          cycle
        end if

        rows = full_ham_overlap_info % rows
        columns = full_ham_overlap_info % columns
        K_load = merge(params % basis % fixed % K, K, params % basis % fixed % enabled == 1)
        root_path = iff(params % basis % fixed % enabled == 1, params % basis % fixed % root_path, params % root_path)
        overlap_block = load_overlap_block(root_path, K_load, K_load, 0, K_sym, slice_ind_row, slice_ind_col, rows, columns, is_file_real)
        ! Determine which part of overlaps block should be stored in the current chunk
        first_row = global_overlap_info % borders % top - full_ham_overlap_info % borders % top + 1
        last_row = size(overlap_block, 1) - (full_ham_overlap_info % borders % bottom - global_overlap_info % borders % bottom)

        local_overlap_info_data => local_k_block_info % subblocks(ir, ic) % extract_from_matrix(this % proc_chunk)
        local_overlap_info_data = overlap_block(first_row : last_row, :)
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads overlaps for this processor's chunk of hamiltonian.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_chunk_overlaps(this, params, ham_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), intent(in) :: ham_info
    integer :: ir, ic, K1_ind, K2_ind

    do ir = 1, size(this % global_chunk_info % subblocks, 1)
      do ic = 1, size(this % global_chunk_info % subblocks, 2)
        K1_ind = this % global_chunk_info % subblocks(ir, ic) % block_row_ind
        K2_ind = this % global_chunk_info % subblocks(ir, ic) % block_col_ind

        ! Only diagonal K blocks have pure overlaps
        if (K1_ind /= K2_ind .or. K1_ind == -1) then
          cycle
        end if

        call this % load_k_block_overlaps(params, this % local_chunk_info % subblocks(ir, ic), this % global_chunk_info % subblocks(ir, ic), ham_info % subblocks(K1_ind, K2_ind))
      end do
    end do
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Multiplies this processor's chunk with kinetic energy matrix elements.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine include_kinetic_energy(this, kinetic)
    class(distributed_rovib_hamiltonian), target, intent(inout) :: this
    complex(real64), intent(in) :: kinetic(:, :) ! Kinetic energy matrix for rho grid
    complex(real64), pointer :: overlap_block(:, :)
    integer :: i, j, n, m, K1_ind, K2_ind, slice_ind_row, slice_ind_col
    type(matrix_block_info), pointer :: local_k_block_info, global_k_block_info, local_overlap_block_info, global_overlap_block_info

    ! iterate through K-blocks
    do i = 1, size(this % local_chunk_info % subblocks, 1)
      do j = 1, size(this % local_chunk_info % subblocks, 2)
        local_k_block_info => this % local_chunk_info % subblocks(i, j)
        global_k_block_info => this % global_chunk_info % subblocks(i, j)
        K1_ind = global_k_block_info % block_row_ind
        K2_ind = global_k_block_info % block_col_ind

        ! only diagonal K blocks need to include kinetic energy
        if (K1_ind /= K2_ind .or. K1_ind == -1) then
          cycle
        end if
        
        ! iterate thorugh overlap blocks
        do n = 1, size(local_k_block_info % subblocks, 1)
          do m = 1, size(local_k_block_info % subblocks, 2)
            local_overlap_block_info => local_k_block_info % subblocks(n, m)
            global_overlap_block_info => global_k_block_info % subblocks(n, m)
            if (local_overlap_block_info % is_empty()) then
              cycle
            end if

            slice_ind_row = global_overlap_block_info % block_row_ind
            slice_ind_col = global_overlap_block_info % block_col_ind
            overlap_block => local_overlap_block_info % extract_from_matrix(this % proc_chunk)
            overlap_block = overlap_block * kinetic(slice_ind_row, slice_ind_col)
          end do
        end do
      end do
    end do
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads 2d eigenvalues from file for given K, symmetry and slice index.
!-------------------------------------------------------------------------------------------------------------------------------------------
  function load_eivals_2d(root_path, K, K_sym, slice_ind) result(eivals)
    character(*), intent(in) :: root_path
    integer, intent(in) :: K, K_sym, slice_ind
    real(real64), allocatable :: eivals(:)
    integer :: file_unit, n_eivals
    character(:), allocatable :: k_path, sym_path, eivals_path

    k_path = get_k_folder_path(root_path, K)
    sym_path = get_sym_path_int(k_path, K_sym)
    eivals_path = get_solutions_2d_path(sym_path, slice_ind)
    open(newunit = file_unit, file = eivals_path, form = 'unformatted')
    read(file_unit) n_eivals
    allocate(eivals(n_eivals))
    read(file_unit) eivals
    close(file_unit)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Adds 2d energies to the portion of main diagonal belonging to this process's chunk.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_k_block_eivals_2d(this, params, local_k_block_info, global_k_block_info, full_ham_k_block_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), target, intent(in) :: local_k_block_info, global_k_block_info, full_ham_k_block_info
    integer :: K_ind, K, K_sym, n, m, K_load, slice_ind_row, slice_ind_col, col_shift, row, col
    real(real64), allocatable :: eivals(:)
    complex(real64), pointer :: overlap_block(:, :)
    character(:), allocatable :: root_path
    type(matrix_block_info), pointer :: local_overlap_block_info, global_overlap_block_info, full_ham_overlap_block_info

    K_ind = global_k_block_info % block_row_ind
    K = k_ind_to_k(K_ind, params % K(1))
    K_sym = get_k_symmetry(K, params % basis % symmetry)

    ! Iterate thorugh n-blocks
    do n = 1, size(local_k_block_info % subblocks, 1)
      do m = 1, size(local_k_block_info % subblocks, 2)
        local_overlap_block_info => local_k_block_info % subblocks(n, m)
        global_overlap_block_info => global_k_block_info % subblocks(n, m)

        if (local_overlap_block_info % is_empty()) then
          cycle
        end if

        slice_ind_row = global_overlap_block_info % block_row_ind
        slice_ind_col = global_overlap_block_info % block_col_ind

        if (slice_ind_row /= slice_ind_col) then
          cycle
        end if

        full_ham_overlap_block_info => full_ham_k_block_info % subblocks(slice_ind_row, slice_ind_col)
        ! First column of the main diagonal may not be first in this process chunk, so we need to make corresponding adjustments
        col_shift = global_overlap_block_info % borders % top - full_ham_overlap_block_info % borders % top

        K_load = merge(params % basis % fixed % K, K, params % basis % fixed % enabled == 1)
        root_path = iff(params % basis % fixed % enabled == 1, params % basis % fixed % root_path, params % root_path)
        eivals = load_eivals_2d(root_path, K_load, K_sym, slice_ind_row)
        overlap_block => local_overlap_block_info % extract_from_matrix(this % proc_chunk)
        do row = 1, size(overlap_block, 1)
          col = row + col_shift
          overlap_block(row, col) = overlap_block(row, col) + eivals(col)
        end do
      end do
    end do
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Adds 2d energies to the portion of main diagonal belonging to this process's chunk.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_chunk_eivals_2d(this, params, ham_info)
    class(distributed_rovib_hamiltonian), target, intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), intent(in) :: ham_info
    integer :: ir, ic, K_row_ind, K_col_ind
    type(matrix_block_info), pointer :: local_k_block_info, global_k_block_info

    ! Iterate through all K-blocks because we don't know where the diagonal blocks are in an arbitrary chunk
    do ir = 1, size(this % local_chunk_info % subblocks, 1)
      do ic = 1, size(this % local_chunk_info % subblocks, 2)
        local_k_block_info => this % local_chunk_info % subblocks(ir, ic)
        global_k_block_info => this % global_chunk_info % subblocks(ir, ic)
        K_row_ind = global_k_block_info % block_row_ind
        K_col_ind = global_k_block_info % block_col_ind

        ! Only diagonal K blocks need to include 2D eigenvalues
        if (K_row_ind /= K_col_ind .or. K_row_ind == -1) then
          cycle
        end if

        call this % load_k_block_eivals_2d(params, local_k_block_info, global_k_block_info, ham_info % subblocks(K_row_ind, K_col_ind))
      end do
    end do
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Adds complex absorbing potential to the portion of main diagonal belonging to this process's chunk.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine include_cap(this, cap, ham_info)
    class(distributed_rovib_hamiltonian), target, intent(inout) :: this
    complex(real64), intent(in) :: cap(:)
    class(matrix_block_info), target, intent(in) :: ham_info
    integer :: ir, ic, n, m, K1_ind, K2_ind, slice_ind_row, slice_ind_col, col_shift, row, col
    complex(real64), pointer :: overlap_block(:, :)
    type(matrix_block_info), pointer :: local_k_block_info, global_k_block_info, local_overlap_block_info, global_overlap_block_info, full_ham_overlap_block_info

    ! Iterate through all K-blocks because we don't know where the diagonal blocks are in an arbitrary chunk
    do ir = 1, size(this % local_chunk_info % subblocks, 1)
      do ic = 1, size(this % local_chunk_info % subblocks, 2)
        local_k_block_info => this % local_chunk_info % subblocks(ir, ic)
        global_k_block_info => this % global_chunk_info % subblocks(ir, ic)
        K1_ind = global_k_block_info % block_row_ind
        K2_ind = global_k_block_info % block_col_ind

        ! Only diagonal K blocks need to include 2D eigenvalues
        if (K1_ind /= K2_ind .or. K1_ind == -1) then
          cycle
        end if
        
        ! Iterate thorugh n-blocks
        do n = 1, size(local_k_block_info % subblocks, 1)
          do m = 1, size(local_k_block_info % subblocks, 2)
            local_overlap_block_info => local_k_block_info % subblocks(n, m)
            global_overlap_block_info => global_k_block_info % subblocks(n, m)

            if (local_overlap_block_info % is_empty()) then
              cycle
            end if
            slice_ind_row = global_overlap_block_info % block_row_ind
            slice_ind_col = global_overlap_block_info % block_col_ind

            ! Only the main diagonal includes CAP
            if (slice_ind_row /= slice_ind_col) then
              cycle
            end if
            overlap_block => local_overlap_block_info % extract_from_matrix(this % proc_chunk)

            ! First column of the main diagonal may not be first in this process chunk, so we need to make corresponding adjustments
            full_ham_overlap_block_info => ham_info % subblocks(K1_ind, K2_ind) % subblocks(slice_ind_row, slice_ind_col)
            col_shift = global_overlap_block_info % borders % top - full_ham_overlap_block_info % borders % top
            do row = 1, size(overlap_block, 1)
              col = row + col_shift
              overlap_block(row, col) = overlap_block(row, col) + cap(slice_ind_row)
            end do
          end do
        end do
      end do
    end do
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads symmetric term contribution for a specific K-block. The block is described by the given block infos.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_k_block_sym_term(this, params, local_k_block_info, global_k_block_info, full_ham_k_block_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), target, intent(in) :: local_k_block_info, global_k_block_info, full_ham_k_block_info
    logical :: is_file_real
    integer :: K_ind, K, K_sym, ir, ic, K_load, slice_row_ind, slice_col_ind, rows, first_row, last_row, J_shift, K_shift, J_factor, K_factor
    complex(real64), allocatable :: overlap_block_J(:, :), overlap_block_K(:, :), overlap_block_slice(:, :)
    complex(real64), pointer :: existing_overlap_block(:, :)
    character(:), allocatable :: root_path
    type(matrix_block_info), pointer :: local_overlap_info, global_overlap_info, full_ham_overlap_info

    K_ind = global_k_block_info % block_row_ind
    K = k_ind_to_k(K_ind, params % K(1))
    K_sym = get_k_symmetry(K, params % basis % symmetry)
    is_file_real = params % use_geometric_phase == 0

    ! Iterate over n-blocks
    do ir = 1, size(local_k_block_info % subblocks, 1)
      do ic = 1, size(local_k_block_info % subblocks, 2)
        local_overlap_info => local_k_block_info % subblocks(ir, ic)
        global_overlap_info => global_k_block_info % subblocks(ir, ic)

        slice_row_ind = global_overlap_info % block_row_ind
        slice_col_ind = global_overlap_info % block_col_ind
        full_ham_overlap_info => full_ham_k_block_info % subblocks(slice_row_ind, slice_col_ind)

        ! Only diagonal n-blocks have this term
        if (slice_row_ind /= slice_col_ind .or. full_ham_overlap_info % is_empty()) then
          cycle
        end if

        call assert(full_ham_overlap_info % rows == full_ham_overlap_info % columns, 'Assertion error: diagonal n-block is not square')
        rows = full_ham_overlap_info % rows

        K_load = merge(params % basis % fixed % K, K, params % basis % fixed % enabled == 1)
        root_path = iff(params % basis % fixed % enabled == 1, params % basis % fixed % root_path, params % root_path)
        overlap_block_J = load_overlap_block(root_path, K_load, K_load, 10, K_sym, slice_row_ind, slice_col_ind, rows, rows, is_file_real)
        overlap_block_K = load_overlap_block(root_path, K_load, K_load, 11, K_sym, slice_row_ind, slice_col_ind, rows, rows, is_file_real)

        ! Determine which part of overlaps block should be stored in the current chunk
        first_row = global_overlap_info % borders % top - full_ham_overlap_info % borders % top + 1
        last_row = size(overlap_block_J, 1) - (full_ham_overlap_info % borders % bottom - global_overlap_info % borders % bottom)
        ! Slice and combine the blocks
        J_shift = merge(params % basis % fixed % J, 0, params % basis % fixed % enabled == 1)
        K_shift = merge(params % basis % fixed % K, 0, params % basis % fixed % enabled == 1)
        J_factor = (params % J * (params % J + 1) - J_shift * (J_shift + 1)) / 2
        K_factor = K ** 2 - K_shift ** 2
        overlap_block_slice = J_factor * overlap_block_J(first_row : last_row, :) + K_factor * overlap_block_K(first_row : last_row, :)

        ! Add to existing regular overlaps
        existing_overlap_block => local_overlap_info % extract_from_matrix(this % proc_chunk)
        existing_overlap_block = existing_overlap_block + overlap_block_slice
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads sym term overlaps for this processor's chunk of hamiltonian.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_chunk_sym_term(this, params, ham_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), intent(in) :: ham_info
    integer :: ir, ic, K_row_ind, K_col_ind

    do ir = 1, size(this % global_chunk_info % subblocks, 1)
      do ic = 1, size(this % global_chunk_info % subblocks, 2)
        K_row_ind = this % global_chunk_info % subblocks(ir, ic) % block_row_ind
        K_col_ind = this % global_chunk_info % subblocks(ir, ic) % block_col_ind

        if (K_row_ind /= K_col_ind .or. K_row_ind == -1) then
          cycle
        end if

        call this % load_k_block_sym_term(params, this % local_chunk_info % subblocks(ir, ic), this % global_chunk_info % subblocks(ir, ic), &
            ham_info % subblocks(K_row_ind, K_col_ind))
      end do
    end do
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads coriolis term contribution for a specific K-block. The block is described by the given block infos.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_k_block_coriolis_term(this, params, local_k_block_info, global_k_block_info, full_ham_k_block_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), target, intent(in) :: local_k_block_info, global_k_block_info, full_ham_k_block_info
    logical :: is_file_real
    integer :: K_row_ind, K_col_ind, K_row, K_col, K_row_sym, ir, ic, K_row_load, K_col_load, slice_ind_row, slice_ind_col, rows, columns, first_row, last_row
    real(real64) :: W
    complex(real64), allocatable :: overlap_block(:, :)
    complex(real64), pointer :: local_overlap_info_data(:, :)
    character(:), allocatable :: root_path
    type(matrix_block_info), pointer :: global_overlap_info, full_ham_overlap_info

    K_row_ind = global_k_block_info % block_row_ind
    K_col_ind = global_k_block_info % block_col_ind
    K_row = k_ind_to_k(K_row_ind, params % K(1))
    K_col = k_ind_to_k(K_col_ind, params % K(1))
    K_row_sym = get_k_symmetry(K_row, params % basis % symmetry)
    is_file_real = params % use_geometric_phase == 0

    if (params % basis % fixed % enabled == 1) then
      ! All K-blocks are the same, so we only need to keep K_row_load different from K_col_load for the correct block transposition 
      ! (alghough even that should not matter since offdiagonal n-blocks are 0)
      K_row_load = merge(params % basis % fixed % K, params % basis % fixed % K + 1, K_row < K_col) ! +1 even if fixed_basis % K == fixed_basis % J, since only lower value is used for loading
      K_col_load = merge(params % basis % fixed % K + 1, params % basis % fixed % K, K_row < K_col)
      root_path = params % basis % fixed % root_path
    else
      K_row_load = K_row
      K_col_load = K_col
      root_path = params % root_path
    end if

    do ir = 1, size(local_k_block_info % subblocks, 1)
      do ic = 1, size(local_k_block_info % subblocks, 2)
        global_overlap_info => global_k_block_info % subblocks(ir, ic)
        slice_ind_row = global_overlap_info % block_row_ind
        slice_ind_col = global_overlap_info % block_col_ind
        full_ham_overlap_info => full_ham_k_block_info % subblocks(slice_ind_row, slice_ind_col)

        ! only diagonal n-blocks have this term
        if (slice_ind_row /= slice_ind_col .or. full_ham_overlap_info % is_empty()) then
          cycle
        end if

        ! Since matrix is symmetric, lower-diagonal blocks are not actually stored on disk and upper diagonal blocks are read and transposed instead
        ! To correct the sign of the block, W-factor should also be taken from the upper-diagonal block even if the current block is lower-diagonal
        if (K_row < K_col) then
          W = calculate_W(params % J, K_row, K_col, params % parity)
        else
          W = calculate_W(params % J, K_col, K_row, params % parity)
        end if

        ! Load block
        rows = full_ham_overlap_info % rows
        columns = full_ham_overlap_info % columns
        overlap_block = load_overlap_block(root_path, K_row_load, K_col_load, 2, K_row_sym, slice_ind_row, slice_ind_col, rows, columns, is_file_real)
        ! Determine which part of overlaps block should be stored in the current chunk
        first_row = global_overlap_info % borders % top - full_ham_overlap_info % borders % top + 1
        last_row = size(overlap_block, 1) - (full_ham_overlap_info % borders % bottom - global_overlap_info % borders % bottom)

        ! Update Hamiltonian
        local_overlap_info_data => local_k_block_info % subblocks(ir, ic) % extract_from_matrix(this % proc_chunk)
        local_overlap_info_data = W * overlap_block(first_row : last_row, :)
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads coriolis overlaps for this processor's chunk of hamiltonian.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_chunk_coriolis_term(this, params, ham_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), intent(in) :: ham_info
    integer :: ir, ic, K_row_ind, K_col_ind

    do ir = 1, size(this % global_chunk_info % subblocks, 1)
      do ic = 1, size(this % global_chunk_info % subblocks, 2)
        K_row_ind = this % global_chunk_info % subblocks(ir, ic) % block_row_ind
        K_col_ind = this % global_chunk_info % subblocks(ir, ic) % block_col_ind

        if (abs(K_row_ind - K_col_ind) /= 1) then
          cycle
        end if
        
        call this % load_k_block_coriolis_term(params, this % local_chunk_info % subblocks(ir, ic), this % global_chunk_info % subblocks(ir, ic), &
            ham_info % subblocks(K_row_ind, K_col_ind))
      end do
    end do
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads asymmetric term contribution for a specific K-block. The block is described by the given block infos.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_k_block_asym_term(this, params, local_k_block_info, global_k_block_info, full_ham_k_block_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), target, intent(in) :: local_k_block_info, global_k_block_info, full_ham_k_block_info
    logical :: is_file_real
    integer :: K_row_ind, K_col_ind, K_row, K_col, K_row_sym, K_row_load, K_col_load, ir, ic, slice_ind_row, slice_ind_col, rows, columns, first_row, last_row
    complex(real64), allocatable :: overlap_block(:, :), overlap_block_slice(:, :)
    complex(real64), pointer :: existing_overlap_block(:, :)
    character(:), allocatable :: root_path
    type(matrix_block_info), pointer :: local_overlap_info, global_overlap_info, full_ham_overlap_info

    K_row_ind = global_k_block_info % block_row_ind
    K_col_ind = global_k_block_info % block_col_ind
    K_row = k_ind_to_k(K_row_ind, params % K(1))
    K_col = k_ind_to_k(K_col_ind, params % K(1))
    K_row_sym = get_k_symmetry(K_row, params % basis % symmetry)
    is_file_real = params % use_geometric_phase == 0

    if (params % basis % fixed % enabled == 1) then
      K_row_load = merge(params % basis % fixed % K, params % basis % fixed % K + 2, K_row <= K_col)
      K_col_load = merge(params % basis % fixed % K + 2, params % basis % fixed % K, K_row <= K_col)
      root_path = params % basis % fixed % root_path
    else
      K_row_load = K_row
      K_col_load = K_col
      root_path = params % root_path
    end if

    do ir = 1, size(local_k_block_info % subblocks, 1)
      do ic = 1, size(local_k_block_info % subblocks, 2)
        local_overlap_info => local_k_block_info % subblocks(ir, ic)
        global_overlap_info => global_k_block_info % subblocks(ir, ic)

        slice_ind_row = global_overlap_info % block_row_ind
        slice_ind_col = global_overlap_info % block_col_ind
        full_ham_overlap_info => full_ham_k_block_info % subblocks(slice_ind_row, slice_ind_col)

        ! only diagonal n-blocks have this term
        if (slice_ind_row /= slice_ind_col .or. full_ham_overlap_info % is_empty()) then
          cycle
        end if

        ! load block
        rows = full_ham_overlap_info % rows
        columns = full_ham_overlap_info % columns
        overlap_block = load_overlap_block(root_path, K_row_load, K_col_load, 3, K_row_sym, slice_ind_row, slice_ind_col, rows, columns, is_file_real)

        ! Determine which part of overlaps block should be stored in the current chunk
        first_row = global_overlap_info % borders % top - full_ham_overlap_info % borders % top + 1
        last_row = size(overlap_block, 1) - (full_ham_overlap_info % borders % bottom - global_overlap_info % borders % bottom)
        overlap_block_slice = overlap_block(first_row : last_row, :)
        overlap_block_slice = calculate_U(params % J, K_row, K_col, params % parity) * overlap_block_slice

        existing_overlap_block => local_overlap_info % extract_from_matrix(this % proc_chunk)
        ! add existing regular overlaps for K = 1
        if (K_row == 1 .and. K_col == 1) then
          existing_overlap_block = existing_overlap_block + overlap_block_slice
        else
          existing_overlap_block = overlap_block_slice
        end if
      end do
    end do
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Loads asym term overlaps for this processor's chunk of hamiltonian.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine load_chunk_asym_term(this, params, ham_info)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    type(matrix_block_info), intent(in) :: ham_info
    integer :: ir, ic, K_row_ind, K_col_ind, K_row, K_col

    do ir = 1, size(this % global_chunk_info % subblocks, 1)
      do ic = 1, size(this % global_chunk_info % subblocks, 2)
        K_row_ind = this % global_chunk_info % subblocks(ir, ic) % block_row_ind
        K_col_ind = this % global_chunk_info % subblocks(ir, ic) % block_col_ind
        K_row = k_ind_to_k(K_row_ind, params % K(1))
        K_col = k_ind_to_k(K_col_ind, params % K(1))

        if (.not. (abs(K_row - K_col) == 2 .or. K_row == 1 .and. K_col == 1)) then
          cycle
        end if

        call this % load_k_block_asym_term(params, this % local_chunk_info % subblocks(ir, ic), this % global_chunk_info % subblocks(ir, ic), &
            ham_info % subblocks(K_row_ind, K_col_ind))
      end do
    end do
  end subroutine
  
!-------------------------------------------------------------------------------------------------------------------------------------------
! Builds a Hamiltonian matrix.
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine build(this, params, kinetic, cap)
    class(distributed_rovib_hamiltonian), intent(inout) :: this
    class(input_params), intent(in) :: params
    complex(real64), intent(in) :: kinetic(:, :) ! Kinetic energy for rho grid
    complex(real64), intent(in), optional :: cap(:) ! Complex absorbing potential for rho
    type(matrix_block_info) :: ham_info

    call load_rovib_hamiltonian_info(params, ham_info) ! loads blocks info for the whole hamiltonian (sizes, positions, borders)
    call this % load_chunk_info(ham_info) ! loads blocks info for the local chunk of the hamiltonian available to this process
    allocate(this % proc_chunk(this % local_chunk_info % rows, this % local_chunk_info % columns)) ! allocate proc chunk
    this % proc_chunk = 0
    call print_parallel('Hamiltonian is allocated')

    call this % load_chunk_overlaps(params, ham_info)
    call print_parallel('Overlaps are loaded')

    call this % include_kinetic_energy(kinetic) ! has to be called when only overlaps are loaded
    call this % load_chunk_eivals_2d(params, ham_info)

    if (params % cap % type /= 'none') then
      call assert(present(cap), 'Error: cap has to be provided if cap type is not none')
      call this % include_cap(cap, ham_info)
    end if

    if (params % basis % fixed % enabled == 1) then
      call this % load_chunk_sym_term(params, ham_info)
      call print_parallel('Sym term is loaded')
    end if

    if (params % use_rovib_coupling == 0) then
      return
    end if

    if (params % debug % enable_terms(1) == 0) then
      call print_parallel('Coriolis term is disabled')
    else
      call this % load_chunk_coriolis_term(params, ham_info)
      call print_parallel('Coriolis term is loaded')
    end if

    if (params % debug % enable_terms(2) == 0) then
      call print_parallel('Asymmetric term is disabled')
    else
      call this % load_chunk_asym_term(params, ham_info)
      call print_parallel('Asymmetric term is loaded')
    end if
  end subroutine

end module
