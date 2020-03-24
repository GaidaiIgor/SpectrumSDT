!-----------------------------------------------------------------------
!  PESGeneral (pesgeneral.f)
!  Contains common things used by different PESes
!  Author: Alexander Teplukhin
!-----------------------------------------------------------------------
module pesgeneral
  use constants
  use formulas_mod
  use input_params_mod
  use general_utils
  
  implicit none
  real*8, parameter :: r0 = 2.2819d0
  ! Potentials on 3D grid
  real*8, allocatable :: potvib(:,:,:) ! Vibrational
  real*8, allocatable :: pottot(:,:,:) ! Total
  ! Masses
  integer mol        ! Molecule code, for example, 686
  real*8 m0          ! Mass of central atom
  real*8 m1          ! Mass of 1st terminal atom
  real*8 m2          ! Mass of 2nd terminal atom
  real*8 mtot        ! Total mass
  real*8 mu          ! Three-particle reduced mass
  real*8 mu0         ! Reduced mass of m1 and m2
  real*8 mu1         ! Reduced mass of m0 and m2
  real*8 mu2         ! Reduced mass of m0 and m1
  ! Rotation
  real*8 jx,jy,jz    ! Components of total angular momentum J
  integer jlarge     ! Absolute value of J
  integer klarge     ! Projection onto Z axis
  ! CAP
  integer,parameter   :: ncap = 3 ! Number of CAPs
  real*8, allocatable :: cap(:,:) ! CAPs on grid
  real*8 capebar     ! Barrier energy used for CAPs
  real*8 dac         ! Delta for strength
  real*8 dwc         ! Delta for exponential
  real*8 drc         ! Delta for position
  integer ndbw       ! Number of de Broglie waves
  real*8 emin        ! Emin for Manolopoulos CAP
  integer capid      ! CAP id
  ! Miscellaneous
  real*8 shift       ! Used for shifting PES down
  real*8 threshold   ! Threshold in the lowest channel

contains
  
!-----------------------------------------------------------------------
!  Initialization of common PES variables.
!  Mass is read as 3-digit value, like "686".
!  Each particular PES must call this subroutine from init_pots.
!-----------------------------------------------------------------------
  subroutine init_pots_general(params)
    implicit none
    class(input_params), intent(in) :: params
    integer :: a1, a2, a3, total_code
    real*8 maxmu

    ! setup global variables
    mol = str2int(params % molecule)
    jlarge = params % J
    klarge = params % K(1)
    if (params % cap_type == 'none') then
      capid = 0
    else if (params % cap_type == 'Manolopoulos') then
      capid = 3
    end if

    ! Other parameters are hardcoded for now. Can be made a part of params if necessary.
    dac = 0
    dwc = 0
    drc = 0
    ndbw = 3
    emin = 7 / autown

    ! Set masses: split the mass number into digits
    a1 = mod(mol/100, 10)
    a2 = mod(mol/10,  10)
    a3 = mod(mol,     10)
    if( a1<6.or.a1>8 .or. a2<6.or.a2>8 .or. a3<6.or.a3>8 ) stop 'Wrong molecule'

    m1 = isomass(a1-5)
    m0 = isomass(a2-5)
    m2 = isomass(a3-5)
    mtot = m0 + m1 + m2
    mu = sqrt(m0 * m1 * m2 / mtot)
    mu0 = m1 * m2 / (m1 + m2)
    mu1 = m0 * m2 / (m0 + m2)
    mu2 = m0 * m1 / (m0 + m1)

    ! select corresponding ZPE. Usage of * over + is to distinguish 677 from 668, etc.
    total_code = a1 * a2 * a3
    if (total_code == 216) then
      shift = -zpe_66
    else if (total_code == 288) then
      shift = -zpe_68
    else if (total_code == 384) then
      shift = -zpe_88
    end if

    ! Set threshold ?
    ! Large reduced mass gives both small ZPE and rotpot in channel
    maxmu = max(mu0,mu1,mu2)
    threshold = klarge**2 / (2*maxmu*r0**2) ! TODO: what is that for?
  end subroutine

  !-----------------------------------------------------------------------
  !  Calculates rotational potential in symmetric top approximation.
  !-----------------------------------------------------------------------
  function calc_potrot(rho, theta, J, K) result(res)
    implicit none
    real*8 :: rho, theta
    integer, optional, intent(in) :: J, K
    real*8 :: res
    res = calc_potrotj2(rho, theta, J) + calc_potrotk2(rho, theta, K)
  end function

  !-----------------------------------------------------------------------
  !  Calculates J**2 component of rotational potential.
  !-----------------------------------------------------------------------
  function calc_potrotj2(rho, theta, J) result(res)
    implicit none
    real*8 :: rho, theta
    integer, optional, intent(in) :: J
    real*8 :: res
    integer :: J_act
    real*8 :: a, b

    a = get_rotational_a(mu, rho, theta)
    b = get_rotational_b(mu, rho, theta)
    J_act = merge(J, jlarge, present(J))
    res = (a + b)/2 * J_act*(J_act + 1)
  end function

  !-----------------------------------------------------------------------
  !  Calculates K**2 component of rotational potential.
  !-----------------------------------------------------------------------
  function calc_potrotk2(rho, theta, K) result(res)
    implicit none
    real*8 :: rho, theta
    integer, optional, intent(in) :: K
    real*8 :: res
    integer :: K_act
    real*8 :: a, b, c

    a = get_rotational_a(mu, rho, theta)
    b = get_rotational_b(mu, rho, theta)
    c = get_rotational_c(mu, rho, theta)
    K_act = merge(K, klarge, present(K))
    res = (c - (a + b)/2) * K_act**2
  end function

  !-----------------------------------------------------------------------
  !  Calculates asymmetric component of rotational potential.
  !-----------------------------------------------------------------------
  real*8 function calc_potrotasym(rho,tet)
    implicit none
    real*8 rho, tet
    real*8 :: a, b
    a = get_rotational_a(mu, rho, tet)
    b = get_rotational_b(mu, rho, tet)
    calc_potrotasym = (a - b) / 2
  end function

  !-----------------------------------------------------------------------
  !  Calculates rotational potential at a given point.
  !-----------------------------------------------------------------------
  real*8 function calc_potrotfull(rho,tet)
    implicit none
    real*8 rho,tet
    calc_potrotfull = ( jx**2/(1-sin(tet)) + jy**2/(1+sin(tet)) + jz**2/(2*sin(tet)**2) ) / (mu*rho**2)
  end function

  !-----------------------------------------------------------------------
  !  Calculates extra potential term at a given point.
  !-----------------------------------------------------------------------
  real*8 function calc_potxtr(rho,tet)
    implicit none
    real*8 rho,tet
    calc_potxtr = -(0.25d0+4/(sin(2*tet)**2))/(2d0*mu*rho**2)
  end function

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
    allocate(cap(nr,ncap))
    cap = 0

    ! Setup parameters for default CAP
    d0 = sqrt( (m0/mu) * (1.0d0 - m0/mtot) )
    ac = (10d0 + dac) * 1000 / autown
    wc = 6d0 * d0 + dwc
    rc = 7d0 * d0 + drc
    write(u,*)'Default CAP start:  ', rc
    write(u,*)'Default CAP length: ', wc

    ! Fill in default CAP
    do ir=1,nr
      r = g(ir)
      if(r > rc)then
        cap(ir,1) = ac * exp ( - wc / ( r - rc ) )
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
    write(u,*)'Balint-Kurti CAP start:  ', rc
    write(u,*)'Balint-Kurti CAP length: ', damplen

    ! Check if grid CAP is too close to covalent well
    if(rc < 5)then
      write(u,*)'Balint-Kurti CAP is too close to covalent well'
      write(u,*)'Required grid end: ', 5 + damplen
    end if

    ! Fill in Balint-Kurti CAP
    do ir=1,nr
      r = g(ir)
      if(r > rc)then
        cap(ir,2) = dampstr * 13.22d0 * exp(-2 * damplen/(r-rc))
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
    write(u,*)'Manolopoulos CAP start:  ', rc
    write(u,*)'Manolopoulos CAP length: ', damplen

    ! Check if grid CAP is too close to covalent well
    if(rc < 5)then
      write(u,*)'Manolopoulos CAP is too close to covalent well'
      write(u,*)'Required grid end: ', 5 + damplen
    end if

    ! Fill in Manolopoulos CAP
    do ir=1,nr
      r = g(ir)
      if(r > rc)then
        xM = 2 * dM * sqrt(2 * mu * eabs) * (r - rc)
        cap(ir,3) = eabs * (aM * xM - bM * xM**3 + 4/(cM-xM)**2 - 4/(cM+xM)**2 )
        if(r > rc + damplen)then
          cap(ir,3) = cap(ir-1,3)
        end if
      end if
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Prints vibrational potential on 3D grid.
  !-----------------------------------------------------------------------
  subroutine prnt_potvib
    implicit none
    integer i1,i2,i3
    open(1,file='potvib.dat',form='unformatted')
    write(1)potvib
    close(1)
    open(1,file='potvib.out')
    do i1=1,size(potvib,3)
      do i2=1,size(potvib,2)
        do i3=1,size(potvib,1)
          write(1,*)potvib(i3,i2,i1)*autown
        end do
      end do
    end do
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  !  Prints total potential on 3D grid.
  !-----------------------------------------------------------------------
  subroutine prnt_pottot
    implicit none
    integer i1,i2,i3
    open(1,file='pottot.dat',form='unformatted')
    write(1)pottot
    close(1)
    open(1,file='pottot.out')
    do i1=1,size(pottot,3)
      do i2=1,size(pottot,2)
        do i3=1,size(pottot,1)
          write(1,*)pottot(i3,i2,i1)*autown
        end do
      end do
    end do
    close(1)
  end subroutine

  !-----------------------------------------------------------------------
  !  Prints CAP.
  !-----------------------------------------------------------------------
  subroutine prnt_cap(fn)
    implicit none
    integer i1,id
    character(*)::fn
    open(1,file=fn)
    do i1=1,size(cap,1)
      write(1,'(I3)',advance='no')i1
      do id=1,ncap
        write(1,'(F25.17)',advance='no')cap(i1,id) * autown
      end do
      write(1,*)
    end do
    close(1)
  end subroutine
end module
