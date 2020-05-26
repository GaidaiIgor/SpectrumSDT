!-------------------------------------------------------------------------------------------------------------------------------------------
! Data and procedures related to complex absorbing potential (CAP)
!-------------------------------------------------------------------------------------------------------------------------------------------
module cap_mod
  use constants
  use general_vars
  use input_params_mod
  use parallel_utils

  integer, parameter :: ncap = 3 ! Number of CAPs
  character(*), parameter :: capdir = 'caps'
  integer :: capid      ! CAP id
  integer :: ndbw       ! Number of de Broglie waves
  real*8 :: capebar     ! Barrier energy used for CAPs
  real*8 :: dac         ! Delta for strength
  real*8 :: dwc         ! Delta for exponential
  real*8 :: drc         ! Delta for position
  real*8 :: emin        ! Emin for Manolopoulos CAP
  real*8, allocatable :: all_caps(:, :) ! CAPs on grid

contains

  !-----------------------------------------------------------------------
  !  Calculates CAP on a given grid for rho
  !  Three types: 1. Adopted from dimentionally reduced
  !               2. Balint-Kurti
  !               3. Manolopoulos
  !-----------------------------------------------------------------------
  subroutine calc_cap(nr,g,u)
    implicit none
    ! CAP #1
    real*8 ac,wc,rc
    real*8 d0
    ! CAP #2
    real*8,parameter :: strs(*) = (/ 1.92, 1.88, 1.85 /)
    ! CAP #3
    real*8,parameter :: aM = 0.112449d0
    real*8,parameter :: bM = 0.00828735d0
    real*8,parameter :: cM = 2.62206d0
    real*8 dM
    real*8 xM
    ! Common
    real*8,parameter :: eabsmin = 1 / autown
    real*8 g(nr)     ! Grid
    real*8 r         ! Running grid point
    real*8 eabs      ! Absolute energy
    real*8 dbl       ! De Broglie wavelength
    real*8 dampstr   ! Damping strength
    real*8 damplen   ! Damping length
    integer ir,nr,u

    ! Allocate array
    allocate(all_caps(nr,ncap))
    all_caps = 0

    ! Setup parameters for default CAP
    d0 = sqrt((m0/mu) * (1.0d0 - m0/mtot))
    ac = (10d0 + dac) * 1000 / autown
    wc = 6d0 * d0 + dwc
    rc = 7d0 * d0 + drc

    ! Fill in default CAP
    do ir=1,nr
      r = g(ir)
      if(r > rc)then
        all_caps(ir,1) = ac * exp ( - wc / ( r - rc ) )
      end if
    end do

    ! Setup parameters for Balint-Kurti CAP
    eabs = capebar
    if(eabs < eabsmin)then
      eabs = eabsmin
    end if
    dbl     = 2 * pi / sqrt(2 * mu * eabs)
    damplen = ndbw       * dbl
    dampstr = strs(ndbw) * Eabs
    rc = g(nr) + (g(nr)-g(nr-1))/2 - damplen + drc

    ! Fill in Balint-Kurti CAP
    do ir=1,nr
      r = g(ir)
      if(r > rc)then
        all_caps(ir,2) = dampstr * 13.22d0 * exp(-2 * damplen/(r-rc))
      end if
    end do

    ! Setup parameters for Manolopoulos CAP
    eabs = emin
    if(eabs < eabsmin)then
      eabs = eabsmin
    end if
    dM   = cM / (4 * pi)
    damplen = cM / (2 * dM * sqrt(2 * mu * eabs))
    rc = g(nr) + (g(nr)-g(nr-1))/2 - damplen + drc

    ! Fill in Manolopoulos CAP
    do ir=1,nr
      r = g(ir)
      if(r > rc)then
        xM = 2 * dM * sqrt(2 * mu * eabs) * (r - rc)
        all_caps(ir,3) = eabs * (aM * xM - bM * xM**3 + 4/(cM-xM)**2 - 4/(cM+xM)**2 )
        if(r > rc + damplen)then
          all_caps(ir,3) = all_caps(ir-1,3)
        end if
      end if
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Prints CAP.
  !-----------------------------------------------------------------------
  subroutine prnt_cap(fn)
    implicit none
    integer i1,id
    character(*)::fn
    open(1,file=fn)
    do i1=1,size(all_caps,1)
      write(1,'(I3)',advance='no')i1
      do id=1,ncap
        write(1,'(F25.17)',advance='no') all_caps(i1,id) * autown
      end do
      write(1,*)
    end do
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  ! Initializes CAPs (complex absorbing potential)
  !-----------------------------------------------------------------------
  subroutine init_caps(params, capebarin)
    implicit none
    class(input_params), intent(in) :: params
    real*8 :: capebarin
    character(256) :: fn

    if (params % cap_type == 'none') then
      capid = 0
    else if (params % cap_type == 'Manolopoulos') then
      capid = 3
    end if

    ! CAP params are hardcoded for now
    dac = 0
    dwc = 0
    drc = 0
    ndbw = 3
    emin = 7 / autown
    capebar = capebarin

    call calc_cap(n1,g1,LG)
    write(fn,'(4A,I5.5,A)')outdir,'/',capdir,'/cap',myid+1,'.out'

    if (get_proc_id() == 0 .and. params % mode == 'diagonalization') then
      call prnt_cap(fn)
    end if
  end subroutine

!-----------------------------------------------------------------------
! Writes active (as specified by global capid) CAP/-i into *p*
!-----------------------------------------------------------------------
  subroutine get_cap(p)
    implicit none
    real*8 :: p(n1)
    p = all_caps(:, capid)
  end subroutine

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns active (as specified by global capid) CAP/-i as real*8
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_real_cap() result(cap)
    real*8, allocatable :: cap(:)
    call assert(allocated(all_caps), 'Error: caps are not initialized')
    allocate(cap(n1))
    call get_cap(cap)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Returns active (as specified by global capid) CAP as complex*16
!-------------------------------------------------------------------------------------------------------------------------------------------
  function get_complex_cap() result(cap)
    complex*16, allocatable :: cap(:)
    real*8, allocatable :: imag_part(:)

    call assert(allocated(all_caps), 'Error: caps are not initialized')
    allocate(cap(n1), imag_part(n1))
    call get_cap(imag_part)
    cap = (0, -1d0) * imag_part
  end function

end module
