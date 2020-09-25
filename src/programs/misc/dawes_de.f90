program main
  use constants
  use general_utils
  implicit none
  real*8 :: coords(3)
  real*8 :: min_point(3), diss_min_point(2)
  real*8 :: min_search_limits(2, 2), diss_search_limits(2)
  real*8 :: De, energy
  
  min_search_limits = reshape([2.4d0, 116d0, 2.5d0, 117d0], [2, 2])
  min_point = find_surface_min_energy(min_search_limits, 5, 1d-15)
  print *, 'Minimum (bond length, angle, energy)'
  print *, min_point
  
  ! coords = [2.40315416155719852042693d0, 2.40315416155719852042693d0, 116.835751764854677503535d0] ! O3 -0.0053200768 84798630    -5.3200768 E-3
  ! ! coords = [900d0, 2.2824356078d0, 90d0] ! O2 9274.99028017325
  ! call IMLS([coords(1), coords(1), coords(3)], energy, 1)
  ! ! print *, min_point(1) - coords(1), min_point(2) - coords(3)
  ! print *, energy
  
  coords = [900d0, 2.5d0, 90d0]
  diss_search_limits = [2.2d0, 2.3d0]
  diss_min_point = find_surface_diss_energy(diss_search_limits, coords, 5, 1d-15)
  print *, 'Dissociation (bond length in bohr, energy)'
  print *, diss_min_point
  
  De = diss_min_point(2) - min_point(3)
  print *, 'De:'
  print *, De

contains

  !-------------------------------------
  ! optimizes value of coord(2) (O2 bond length) on Dawes surface
  !-------------------------------------
  function find_surface_diss_energy(search_limits, coords, npoints, eps) result(min_point)
    real*8, value :: search_limits(2), coords(3)
    integer, intent(in) :: npoints
    real*8, intent(in) :: eps
    real*8 :: min_point(2)
    real*8 :: energy(npoints)
    real*8 :: grid(npoints)
    real*8 :: search_length
    integer :: min_index
    integer :: i, from_index, to_index
    
    do while (.true.)
      ! generate grid
      grid(:) = linspace(search_limits(1), search_limits(2), npoints)
      
      ! compute potential
      do i = 1, npoints
        coords(2) = grid(i)
        call IMLS(coords, energy(i), 1)
      end do
      
      min_index = minloc(energy, 1)
      from_index = max(min_index - 1, 1)
      to_index = min(min_index + 1, size(grid))
      search_limits(1) = grid(from_index)
      search_limits(2) = grid(to_index)
      search_length = search_limits(2) - search_limits(1)
      
      if (search_length < eps) then
        min_point(1) = grid(min_index)
        min_point(2) = energy(min_index)
        return
      end if
    end do
  end function

  !-------------------------------------
  ! finds miniumum on Dawes surface
  !-------------------------------------
  function find_surface_min_energy(search_limits, npoints, eps) result(min_point)
    real*8, value :: search_limits(2, 2)
    integer, intent(in) :: npoints
    real*8, intent(in) :: eps
    real*8 :: min_point(3)
    real*8 :: energy(npoints, npoints)
    real*8 :: grids(2, npoints)
    real*8 :: coords(3), search_length(2)
    integer :: min_index(2)
    integer :: i, j, from_index, to_index
    
    do while (.true.)
      ! generate grids for each dimension
      do i = 1, 2
        grids(i, :) = linspace(search_limits(i, 1), search_limits(i, 2), npoints)
      end do
      
      ! compute 2D potential
      do i = 1, npoints
        coords(1) = grids(1, i)
        coords(2) = grids(1, i)
        do j = 1,npoints
          coords(3) = grids(2, j)
          call IMLS(coords, energy(i, j), 1)
        end do
      end do
      
      min_index = minloc(energy)
      do i = 1, 2
        from_index = max(min_index(i) - 1, 1)
        to_index = min(min_index(i) + 1, size(grids, 2))
        search_limits(i, 1) = grids(i, from_index)
        search_limits(i, 2) = grids(i, to_index)
        search_length(i) = search_limits(i, 2) - search_limits(i, 1)
      end do
      
      ! print *, 'Next iteration'
      ! call print_real_matrix(grids)
      ! print *
      ! call print_real_matrix(energy)
      ! print *
      ! print *, min_index
      
      ! do i = 1, 2
        ! min_point(i) = grids(i, min_index(i))
      ! end do
      ! min_point(3) = energy(min_index(1), min_index(2))
      ! print *, min_point
      
      if (all(search_length < eps)) then
        do i = 1, 2
          min_point(i) = grids(i, min_index(i))
        end do
        min_point(3) = energy(min_index(1), min_index(2))
        return
      end if
    end do
  end function

  !-------------------------------------
  ! finds miniumum on 3d surface
  !-------------------------------------
  function find_surface_min_3d(search_limits, npoints, eps) result(min_point)
    real*8, value :: search_limits(3, 2)
    integer, intent(in) :: npoints
    real*8, intent(in) :: eps
    real*8 :: min_point(4)
    real*8 :: energy(npoints, npoints, npoints)
    real*8 :: grids(3, npoints)
    real*8 :: coords(3), search_length(3)
    integer :: min_index(3)
    integer :: i, j, k
    
    do while (.true.)
      ! generate grids for each dimension
      do i = 1, 3
        grids(i, :) = linspace(search_limits(i, 1), search_limits(i, 2), npoints)
      end do
      
      ! compute 3D potential
      do i = 1, npoints
        coords(1) = grids(1, i)
        do j = 1,npoints
          coords(2) = grids(2, j)
          do k = 1,npoints
            coords(3) = grids(3, k)
            call IMLS(coords, energy(i, j, k), 1)
          end do
        end do
      end do
      
      min_index = minloc(energy)
      do i = 1, 3
        search_limits(i, 1) = grids(i, min_index(i) - 1)
        search_limits(i, 2) = grids(i, min_index(i) + 1)
        search_length(i) = search_limits(i, 2) - search_limits(i, 1)
      end do
      
      if (all(search_length < eps)) then
        do i = 1, 3
          min_point(i) = grids(i, min_index(i))
        end do
        min_point(4) = energy(min_index(1), min_index(2), min_index(3))
        return
      end if
    end do
  end function
  
end program
