!-----------------------------------------------------------------------
!  ZPE of O2 (zpeO2.f)
!  Calculates a vibrational spectrum of O2
!  The coordinate is an interatomic distance
!  Author: Alexander Teplukhin
!-----------------------------------------------------------------------
       use config
       use constants
       use general_vars
       use input_params_mod
       use mkl
       use pesinterface_mod
       use lapack95
       implicit none
       real*8,allocatable::g(:),freq(:),pot(:)
       real*8,allocatable::eivals(:),eivecs(:,:),ham(:,:),w(:)
       real*8 rmin,rmax,rl,js
       real*8 potmin
       integer n,nstate,costheta
       type(input_params) :: params

       call init_parameters_zpe_o2(params)
       write(*,*)'Total potential optimization'
       call find_minimum(calc_ptot)
       write(*,*)'Vibrational potential optimization'
       call find_minimum(calc_pvib)
       call load_grid_pottot
       call init_matrix
       allocate(eivals(nstate),eivecs(n,nstate),w(n))
       call syev(ham,w,'V')
       eivals = w(1:nstate)
       eivecs = ham(:,1:nstate)
       call print_pot
       call print_eivals
       call print_eivecs

       contains

!-----------------------------------------------------------------------
!  Input parameters.
!-----------------------------------------------------------------------
       subroutine init_parameters_zpe_o2(params)
       implicit none
       type(input_params), intent(inout) :: params
       
       ! open(1,file='zpeO2.config')
       params = process_user_settings('spectrumsdt.config')
       call init_pots(params)
       n = 64
       nstate = 16
       rl = 60
       rmin = 1.972625
       rmax = 2.778875
       costheta = 0
       js = params % J
       end subroutine

!-----------------------------------------------------------------------
!  Finds PES minimum using Newton-Raphson method.
!-----------------------------------------------------------------------
       subroutine find_minimum(potfunc)
       implicit none
       real*8 r1,r2,r3,r4
       real*8 f1,f2,f3,df,ds
       real*8 potfunc
       r2 = 2.282d0
       r1 = r2 - 1d-1
       r3 = r2 + 1d-1
       do
         f1 = potfunc(r1,rl,costheta)
         f2 = potfunc(r2,rl,costheta)
         f3 = potfunc(r3,rl,costheta)
         df = (f3 - f1) / (r3 - r1)
         ds = (f3 - 2*f2 + f1) / (r2 - r1)**2
         r4 = r2 - df / ds
         potmin = potfunc(r4,rl,costheta)
         write(*,*)r4,potmin*autown
         if(abs(r4-r2)<1E-10)exit
         r1 = r2
         r2 = r4
         r3 = 2*r2-r1
       enddo
       write(*,*)'Minimum: ',r4,potmin*autown
       end subroutine

!-----------------------------------------------------------------------
!  Loads grid and potential.
!-----------------------------------------------------------------------
       subroutine load_grid_pottot
       implicit none
       real*8 l
       integer i
       allocate(g(n),pot(n))
       l = rmax - rmin
       do i=1,n
         g(i) = rmin + l/n * (i-0.5d0)
         pot(i) = calc_ptot(g(i),rl,costheta) - potmin
       enddo
       close(1)
       end subroutine

!-----------------------------------------------------------------------
!  Calculates vibraional potential.
!-----------------------------------------------------------------------
       real*8 function calc_pvib(rs,rl,costheta)
       implicit none
       real*8 r(3),r2(3),cart(9)
       real*8 rs,rl,vpot
       integer costheta
       r2(1) = rs
       r2(2) = rl
       r2(3) = costheta*1d0 ! cos(0 or pi/2)
       call INT_Cart(r2,cart,pes_mass,1)
       call cart_int(r2,cart,pes_mass,2)
       r(1)=min(r2(1),r2(2))
       r(2)=max(r2(1),r2(2))
       r(3)=180d0*acos(r2(3))/pi
       call IMLS(r,vpot,1)
       calc_pvib = vpot / autown
       end function

!-----------------------------------------------------------------------
!  Calculates rotational potential.
!-----------------------------------------------------------------------
       real*8 function calc_prot(rs)
       implicit none
       real*8 rs
       calc_prot = js*(js+1) / (2d0*mu0*rs**2)
       end function

!-----------------------------------------------------------------------
!  Calculates total potential.
!-----------------------------------------------------------------------
       real*8 function calc_ptot(rs,rl,costheta)
       implicit none
       real*8 rs,rl
       integer costheta
       calc_ptot = calc_pvib(rs,rl,costheta) + calc_prot(rs)
       end function

!-----------------------------------------------------------------------
!  Computes the Hamiltonian matrix.
!-----------------------------------------------------------------------
       subroutine init_matrix
       implicit none
       real*8 psi(n),m
       integer i
       ! Init kinetic operator
       allocate(freq(n))
       call init_derivd(2,n,rmax-rmin,freq)
       ! Init matrix
       allocate(ham(n,n))
       do i=1,n
         psi = 0d0
         psi(i) = 1d0
         call calc_derivd(2,1,n,freq,psi)
         psi = - psi / (2d0*mu0)
         psi(i) = psi(i) + pot(i)
         ham(:,i) = psi
       enddo
       end subroutine

!-----------------------------------------------------------------------
!  Prints computed potential
!-----------------------------------------------------------------------
       subroutine print_pot
       implicit none
       integer i
       open(1,file='pot.out')
       do i=1,n
         write(1,'(2F25.17)')g(i),pot(i)*autown
       enddo
       close(1)
       end subroutine

!-----------------------------------------------------------------------
!  Prints computed spectrum with symmetries.
!-----------------------------------------------------------------------
       subroutine print_eivals
       implicit none
       integer k
       open(1,file='eivals.out')
       do k=1,nstate
         write(1,'(I4,2F25.17)')k,eivals(k)*autown,symmetry(eivecs(:,k))
       enddo
       close(1)
       end subroutine

!-----------------------------------------------------------------------
!  Prints computed wave functions, squared value of the wave function.
!-----------------------------------------------------------------------
       subroutine print_eivecs
       implicit none
       integer k,istate
       character*32 fn
       do istate=1,nstate
         write(fn,'(A5,I3.3,A4)')'state',istate,'.out'
         open(1,file=fn)
         do k=1,n
           write(1,'(F25.16)')eivecs(k,istate)
         enddo
         close(1)
       enddo
       end subroutine

!-----------------------------------------------------------------------
!  Computes the symmetry of the wave function.
!-----------------------------------------------------------------------
       real*8 function symmetry(psi)
       implicit none
       real*8 psi(n),s
       integer i
       s = 0
       do i=1,n
         s = s + psi(i)*psi(n-i+1)
       enddo
       symmetry = s
       end function

       end

