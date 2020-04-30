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
  use general_vars
  
  implicit none
  real*8, parameter :: r0 = 2.2819d0
  ! Potentials on 3D grid
  real*8, allocatable :: potvib(:,:,:) ! Vibrational
  real*8, allocatable :: pottot(:,:,:) ! Total
  ! Rotation
  real*8 jx,jy,jz    ! Components of total angular momentum J
  integer jlarge     ! Absolute value of J
  integer klarge     ! Projection onto Z axis
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
end module
