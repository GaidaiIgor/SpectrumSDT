!-----------------------------------------------------------------------
!  PESAPH (dawes/pesinterface.f)
!  Containes all things related to Dawes PES
!  Author: Alexander Teplukhin, Igor Gayday
!-----------------------------------------------------------------------
 module pesinterface_mod
   use constants
   use input_params_mod
   use mpi_f08
   use parallel_utils
   use pesgeneral

   implicit none
   real*8 :: pes_mass(3)

   private
   public :: pes_mass
   public :: init_pots, calc_potvib, calc_pots

 contains

!-----------------------------------------------------------------------
!  Initialization.
!-----------------------------------------------------------------------
  subroutine init_pots(params)
    implicit none
    class(input_params), intent(in) :: params
    
    call init_pots_general(params)
    pes_mass(1) = m0
    pes_mass(2) = m1
    pes_mass(3) = m2
  end subroutine

!-----------------------------------------------------------------------
! Calculates vibrational potential at given point.
! init_pots has to be called first to set pes_mass
!-----------------------------------------------------------------------
  real*8 function calc_potvib(rho,tet,phi)
    implicit none
    real*8 rho,tet,phi
    real*8 vpot,r(3),r2(3),cart(9)
    r2(1) = rho
    r2(2) = tet
    r2(3) = phi

    call INT_Cart(r2,cart,pes_mass,4)
    call cart_int(r2,cart,pes_mass,2)
    
    r(1)=min(r2(1),r2(2))
    r(2)=max(r2(1),r2(2))
    r(3)=180d0*acos(r2(3))/pi

    call IMLS(r, vpot, 1)
    vpot = vpot / autown - De_dawes
    calc_potvib = vpot + shift
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Converts a 1D index of 1D-representation into 3 indexes corresponding to it in 3D representation
! Assumes the 3D array was flattened using the following order of dimenisions: 3, 2, 1 (3rd coordinate is changing most frequently)
! n1, n2, n3 - sizes of the 3D array
!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine convert_1d_ind_to_3d(ind_1d, n1, n2, n3, i1_3d, i2_3d, i3_3d)
    integer, intent(in) :: ind_1d, n1, n2, n3
    integer, intent(out) :: i1_3d, i2_3d, i3_3d
    integer :: ind_1d_0, i1_3d_0, i2_3d_0, i3_3d_0

    ! Convert using 0-based indexing
    ind_1d_0 = ind_1d - 1
    i1_3d_0 = ind_1d_0 / (n2 * n3)
    i2_3d_0 = (ind_1d_0 - i1_3d_0 * n2 * n3) / n3
    i3_3d_0 = ind_1d_0 - i1_3d_0 * n2 * n3 - i2_3d_0 * n3

    ! Shift back to 1-based indexing
    i1_3d = i1_3d_0 + 1
    i2_3d = i2_3d_0 + 1
    i3_3d = i3_3d_0 + 1
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Calculates potential energy surface
!-------------------------------------------------------------------------------------------------------------------------------------------
   subroutine calc_pots(grid1, grid2, grid3)
     real*8 :: grid1(:), grid2(:), grid3(:)
     integer :: i1, i2, i3, n1, n2, n3, proc_k, first_k, k, proc_points, total_points, proc_id
     integer, allocatable :: proc_counts(:), proc_shifts(:)
     real*8, allocatable :: proc_pot_vib_1d(:), global_pot_vib_1d(:)

     n1 = size(grid1)
     n2 = size(grid2)
     n3 = size(grid3)
     total_points = n1 * n2 * n3
     call get_proc_elem_range(total_points, first_k, proc_points, proc_counts, proc_shifts)
     allocate(proc_pot_vib_1d(proc_points))

     ! Calculate a chunk
     do proc_k = 1, proc_points
       if (get_proc_id() == 0) then
         call track_progress(proc_k * 1d0 / proc_points, 0.01d0)
       end if

       k = first_k + proc_k - 1
       call convert_1d_ind_to_3d(k, n1, n2, n3, i1, i2, i3)
       proc_pot_vib_1d(proc_k) = calc_potvib(grid1(i1), grid2(i2), grid3(i3))
     end do

     ! Allocate global storage on 0th proc
     proc_id = get_proc_id()
     if (proc_id == 0) then
       allocate(global_pot_vib_1d(total_points))
       allocate(potvib(n3, n2, n1))
     end if

     call print_parallel('Exchanging data...')
     call MPI_Gatherv(proc_pot_vib_1d, size(proc_pot_vib_1d), MPI_DOUBLE_PRECISION, global_pot_vib_1d, proc_counts, proc_shifts, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD)

     if (proc_id == 0) then
       k = 1
       do i1 = 1, n1
         do i2 = 1, n2
           do i3 = 1, n3
             potvib(i3, i2, i1) = global_pot_vib_1d(k)
             k = k + 1
           end do
         end do
       end do
     end if
   end subroutine

end module
