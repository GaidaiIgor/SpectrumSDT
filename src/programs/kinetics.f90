!-----------------------------------------------------------------------
! Kinetics (kinetics.f)
! Kinetic model of ozone recombination.
! Calculation of isotopic effects.
! Author: Alexander Teplukhin
!-----------------------------------------------------------------------
program kinetics
  use pesgeneral
  implicit none
  real*8,parameter :: temp = 298
  real*8,parameter :: kltown = 0.6950345740d0
  real*8,parameter :: kltojo = 1.3806485279d-23
  real*8,parameter :: kltoha = 3.166810518d-6
  real*8,parameter :: ktwn = temp * kltown
  real*8,parameter :: ktjo = temp * kltojo
  real*8,parameter :: ktha = temp * kltoha
  real*8,parameter :: autocm = 0.52917721067d-8
  real*8,parameter :: mconc = 6.53E+18
  real*8,parameter :: gatos1 = 1 / 2.418884326505d-17

  ! Output file unit
  integer,parameter  :: OU = 9

  ! Number of channels wihtin group
  integer, parameter :: nch = 5

  ! Config
  character(:),allocatable :: chpath ! Channels path
  character(:),allocatable :: sppath ! Spectrum path
  integer jlmin,jlmax,jlstep ! J data
  integer klmax,klstep       ! K data
  integer level              ! Level
  logical loadkstab          ! Loads kstab from config
  logical gafilter           ! Gamma filter
  real*8  pwmin              ! pwmin
  real*8  kmin               ! kmin
  integer usech              ! Use only specific channel or all
  integer enfilter           ! Energy filter type: 1 - de, 2 - en
  real*8  enmax              ! Maximum energy

  ! J&K variables
  integer njl,ijl,jl
  integer nkl,ikl,kl

  ! Computed constants
  real*8 erca                            ! Energy of reactants, A
  real*8 ercb                            ! Energy of reactants, B
  real*8 qrca                            ! P-func of reactants, A
  real*8 qrcb                            ! P-func of reactants, B
  real*8 kex                             ! K exchange
  real*8 kbfactor                        ! Kex or 1/Kex

  ! Channel data
  real*8,  allocatable :: chb(:,:,:,:,:) ! Channel barrier position
  real*8,  allocatable :: che(:,:,:,:,:) ! Channel barrier energy
  real*8,  allocatable :: chc(:,:,:,:,:) ! Channel barrier coeff
  real*8,  allocatable :: chi(:,:,:,:,:) ! Channel numbers

  ! Common data
  real*8,  allocatable :: ets(:,:,:,:) ! Energy of TS
  real*8,  allocatable :: edd(:,:,:,:) ! Energy double dagger
  real*8,  allocatable :: qts(:,:,:,:) ! P-func of TS
  real*8,  allocatable :: kre(:,:,:,:) ! krec
  real*8                  rsf(nch)     ! Resonances factors for Qts

  ! Groups
  integer,   parameter :: ngr = 3
  character, parameter :: grname(*) = (/ 'B', 'A', 'S' /)
  integer    igr

  ! Symmetries
  integer,   parameter :: nsy = 2
  character, parameter :: sys(*) = (/ 'S', 'A' /)
  character  sy
  integer    isy

  ! Running channel number
  integer ich
  integer ichi

  ! Dissociation channels in      686          868
  integer dclo      ! Lowest     A (2)   or   B (1)
  integer dchi      ! Highest    B (1)   or   A (2)

  ! Miscellaneous
  character(*),parameter   :: outdir = 'kinetics'
  character(:),allocatable :: dir
  character(256) str1
  character(256) str2
  real*8  kstab
  real*8  kdiss
  real*8  sstab
  logical ex

  ! Setup directory
  inquire(directory=outdir,exist=ex)
  if(ex)then
   stop 'Output directory already exists'
  else
   call execute_command_line('mkdir ' // outdir)
  endif

  ! Load config
  open(1,file='kinetics.config')
  read(1,*)jlmin,jlmax,jlstep
  read(1,*)klmax,klstep
  read(1,*)loadkstab,kstab
  read(1,*)level
  read(1,'(A)')str1
  read(1,'(A)')str2
  read(1,*)gafilter,pwmin,kmin,usech
  read(1,*)enfilter,enmax
  close(1)
  njl = (jlmax - jlmin) / jlstep + 1
  nkl = klmax / klstep + 1
  chpath = trim(str1)
  sppath = trim(str2)

  ! The code only works for klstep = 2
  if(klstep /= 2)stop 'klstep must be 2'

  ! Allocate JK arrays
  allocate(ets(nkl,njl,ngr,nsy), edd(nkl,njl,ngr,nsy), qts(nkl,njl,ngr,nsy), kre(nkl,njl,ngr,nsy))

  ! Set default resonance factors
  rsf = 1

  ! Open output file for summary
  open(OU,file=outdir//'/summary.out')

  ! Calculate constants
  write(*,*)'Calculating constants'
  call calc_const
  write(*,*)

  ! Collect channels
  write(*,*)'Collecting channels'
  call coll_chan
  write(*,*)

  ! Calculate kstab
  if(loadkstab)then
    write(*,*)'Skipping kstab calculation'
    call calc_kdiss
  else
    write(*,*)'Calculating kstab'
    call calc_kstab
  endif
  write(*,*)

  ! Calculate effects
  write(*,*)'Calculating effects'
  call calc_eff
  contains

!-----------------------------------------------------------------------
! Returns directory name.
!-----------------------------------------------------------------------
  subroutine get_dir
    character(256) str
    logical ex

    ! Get directory name
    write(str,'(A,I2.2,A,I2.2,A)')'J',jl,'/K',kl,sy
    dir = trim(str)
    write(*,'(A)',advance='no')dir

    ! Check if directory exists
    inquire(directory=dir,exist=ex)
    if(ex)then
      write(*,*)
    else
      write(*,*)'Skipped'
      dir = ''
    endif
  end subroutine

!-----------------------------------------------------------------------
! Calculates constants.
!-----------------------------------------------------------------------
  subroutine calc_const
    real*8,parameter :: dzpeO2(*) = [22.26698417, 22.93902199] ! 686, 868

    ! Variables for Qtr
    real*8 muOO2
    real*8 dbw
    real*8 qtra
    real*8 qtrb

    ! Variables for Qrt
    integer j ! j small running
    real*8  e ! Energy
    real*8  q ! Partition function
    real*8  qrta
    real*8  qrtb

    ! Get directory to read masses from
    sy = sys(1)
    jl = jlmin
    kl = 0
    call get_dir
    if(len(dir) == 0) stop 'No directory to read masses'

    ! Read masses
    open(1,file=dir//'/spectrumsdt.config')
    call init_pots_general(1)
    close(1)

    ! Translational partition function of O2 in channel A
    muOO2 = (m0 + m2) * m1 / mtot
    dbw   = sqrt(2 * pi / (muOO2 * ktha)) * autocm
    qtra  = 1 / dbw**3
    write(OU,'(A,ES35.17E3)')'Qtra:',qtra

    ! Translational partition function of O2 in channel B
    muOO2 = (m1 + m2) * m0 / mtot
    dbw   = sqrt(2 * pi / (muOO2 * ktha)) * autocm
    qtrb  = 1 / dbw**3
    write(OU,'(A,ES35.17E3)')'Qtrb:',qtrb

    ! Rotational partition function of O2 in channel A
    qrta = 0
    do j=0,100,1
      e = rcerot(murc(),j)
      qrta = qrta + (2 * j + 1) * exp(-e / ktwn)
    enddo
    write(OU,'(A,ES35.17E3)')'Qrta:',qrta

    ! Rotational partition function of O2 in channel B
    qrtb = 0
    do j=1,100,2
      e = rcerot(murc(),j) - rcerot(murc(),1)
      qrtb = qrtb + (2 * j + 1) * exp(-e / ktwn)
    enddo
    write(OU,'(A,ES35.17E3)')'Qrtb:',qrtb

    ! Qrca
    qrca = qtra * qrta
    write(OU,'(A,ES35.17E3)')'Qrca:',qrca

    ! Qrcb
    qrcb = qtrb * qrtb
    write(OU,'(A,ES35.17E3)')'Qrcb:',qrcb

    ! Set energy of reactants, threshold
    if(mol == 686)then
      erca = 0
      ercb = dzpeO2(1) + rcerot(murc(),1)
    else if (mol == 868)then
      erca = dzpeO2(2) - rcerot(murc(),1)
      ercb = 0
    endif
    write(OU,'(A,F35.17)')'Erca:',erca
    write(OU,'(A,F35.17)')'Ercb:',ercb

    ! Calculate Kex
    if(mol == 686)then
      kex = exp(-ercb / ktwn) * qrcb / qrca
    else if(mol == 868)then
      kex = exp(-erca / ktwn) * qrca / qrcb
    endif
    write(OU,'(A,F35.17)')'Kex:',kex

    ! Set lowest and highest diss channels
    if(mol == 686)then
      dclo = 2
      dchi = 1
    else if (mol == 868)then
      dclo = 1
      dchi = 2
    endif

    ! Calculate kB factor, Kex or 1/Kex
    kbfactor = kex**((-1)**dclo)
  end subroutine

!-----------------------------------------------------------------------
! Collects channel data.
!-----------------------------------------------------------------------
  subroutine coll_chan
    ! Channel properties
    integer ichan
    real*8  barps
    real*8  baren
    real*8  barom

    ! Miscellaneous
    integer ncht(3)
    integer i
    integer ios
    integer ist
    character(256) fn

    ! Allocate arrays
    allocate(chb(nch,nkl,njl,ngr,nsy), che(nch,nkl,njl,ngr,nsy), chc(nch,nkl,njl,ngr,nsy), chi(nch,nkl,njl,ngr,nsy))
    chb = 0
    che = 0
    chc = 0
    chi = 0

    ! Loop over symmetries, Js and Ks
    do isy=1,nsy
      sy = sys(isy)
      do ijl=1,njl
        jl = jlmin + (ijl - 1) * jlstep
        do ikl=1,nkl
          kl = (ikl - 1) * klstep

          ! Skip K > J
          if (kl > jl) cycle

          ! Get directory name or cycle if there is no such directory
          call get_dir
          if(len(dir) == 0)cycle

          ! Setup channel counter
          ncht = 0

          ! Loop over channels
          open(1,file=chpath//'/'//dir//'/chrecog/channels.dat')
          do
            ! Read channel
            read(1,'(2I35,2F35.17,70X,F35.17)',iostat=ios) ichan,igr,barps,baren,barom
            if(ios /= 0)exit

            ! Skip if group is full
            if( ncht(igr) == nch ) cycle

            ! Increment number of channels
            ncht(igr) = ncht(igr) + 1

            ! Convert parabola coeff to barrier omega in wn
            barom = sqrt(-2 * barom / (autown * mu)) * autown

            ! Save properties
            ist = ncht(igr)
            chb(ist,ikl,ijl,igr,isy) = barps
            che(ist,ikl,ijl,igr,isy) = baren
            chc(ist,ikl,ijl,igr,isy) = barom
            chi(ist,ikl,ijl,igr,isy) = ichan
          enddo
          close(1)

          ! Stop if not enough channels
          if( sum(ncht) /= 3 * nch ) stop 'Not enough channels'
        enddo
      enddo
    enddo
  end subroutine

!-----------------------------------------------------------------------
! Calculates kstab and sstab.
!-----------------------------------------------------------------------
  subroutine calc_kstab
    real*8,parameter :: mM = 39.962383123d0 * amutoau
    real*8,parameter :: autokg = 9.10938291d-31             
    real*8 krec
    real*8 krel
    real*8 mustab
    real*8 vel

    ! Set parameters
    igr = 3
    kstab = 1

    ! Nullify arrays
    ets = 0
    edd = 0
    qts = 0
    krec = 0

    ! Loop over symmetries, Js and Ks
    do isy=1,nsy
      sy = sys(isy)
      do ijl=1,njl
        jl = jlmin + (ijl - 1) * jlstep
        do ikl=1,nkl
          kl = (ikl - 1) * klstep

          ! Skip K > J
          if (kl > jl) cycle

          ! Get directory name or cycle if there is no such directory
          call get_dir
          if(len(dir) == 0)cycle

          ! Calculate individual properties
          ets(ikl,ijl,igr,isy) = calc_ets()
          edd(ikl,ijl,igr,isy) = calc_edd()
          qts(ikl,ijl,igr,isy) = calc_qts()

          ! Add contribution to krec
          krec = krec + calc_kre()
        enddo
      enddo
    enddo

    ! Multiply krec by JK box size
    krec = krec * jlstep * klstep

    ! Velocity
    mustab = mtot * mM / (mtot + mM)
    vel = sqrt( 8 * ktjo / (pi * mustab * autokg) ) * 1d2 ! cm / s

    ! Calculate effective kstab and sstab
    if(mol == 686)then
      krel = 1.007d0
    else
      krel = 1.04d0
    endif
    kstab = 60d-35 * krel / (2 * krec)
    sstab = kstab / (autocm**2 * vel)
    call calc_kdiss
    write(OU,'(A,ES35.17E3)')'kstab:',kstab
    write(OU,'(A,ES35.17E3)')'kdiss:',kdiss
    write(OU,'(A,ES35.17E3)')'sstab:',sstab
  end subroutine

!-----------------------------------------------------------------------
! Calculates kdiss.
!-----------------------------------------------------------------------
  subroutine calc_kdiss
    kdiss = kstab * 9.59d0 / 123.92d0
  end subroutine

!-----------------------------------------------------------------------
! Calculates effects.
!-----------------------------------------------------------------------
  subroutine calc_eff
    ! Effects distributions
    real*8, allocatable :: dze(:,:)     ! Delta ZPE effect
    real*8, allocatable :: eta(:,:)     ! Eta effect
    real*8, allocatable :: dze_an(:,:)  ! Delta ZPE effect, analytic
    real*8, allocatable :: eta_an(:,:)  ! Eta effect, analytic
    real*8, allocatable :: rte(:,:,:)   ! Ratios of exp(ets)
    real*8, allocatable :: rtq(:,:,:)   ! Ratios of Q

    ! Total effects
    real*8 total_dze
    real*8 total_eta
    real*8 dkl
    real*8 krec(5,3) ! 1 - B
                     ! 2 - A
                     ! 3 - S
                     ! 4 - dchi *   Kex for dZPE-effect
                     ! 5 - B * or / Kex for  eta-effect
                     ! Second index: sym, asym and total

    ! Allocate JK arrays
    allocate(dze(nkl,njl), eta(nkl,njl), dze_an(nkl,njl), eta_an(nkl,njl), rte(nkl,njl,ngr), rtq(nkl,njl,ngr))

    ! Nullify arrays
    ets = 0
    edd = 0
    qts = 0
    kre = 0
    dze = 0
    eta = 0
    dze_an = 0
    eta_an = 0
    rte = 0
    rtq = 0
    krec = 0

    ! Preparation for levels /= 1
    if(level /= 1)then
      open(11,file=outdir//'/res.B.out')
      open(12,file=outdir//'/res.A.out')
      open(13,file=outdir//'/res.S.out')
    endif

    ! Extra preparation for level 4
    if(level == 4)then
      open(20,file=outdir//'/demo.zpe.out')
      open(21,file=outdir//'/demo.zpe.a1.out')
      open(22,file=outdir//'/demo.zpe.b1.out')
      open(30,file=outdir//'/mol.sym.out')
      open(31,file=outdir//'/mol.sym.a1.out')
      open(32,file=outdir//'/mol.sym.b1.out')
      open(40,file=outdir//'/mol.asym.out')
      open(41,file=outdir//'/mol.asym.a1.out')
      open(42,file=outdir//'/mol.asym.b1.out')
    endif

    ! Loop Js and Ks
    do ijl=1,njl
      jl = jlmin + (ijl - 1) * jlstep
      do ikl=1,nkl
        kl = (ikl - 1) * klstep

        ! Skip K > J
        if (kl > jl) cycle

        ! Loop over symmetries
        do isy=1,nsy
          sy = sys(isy)

          ! Determine dkl
          if(kl == jl)then
            dkl = 1.5d0
          elseif(kl == 0)then
            if(isy == 1)then
              dkl = 1.5d0
            else
              dkl = 0.5d0
            endif
          else
            dkl = 2d0
          endif

          ! Get directory name or cycle if there is no such directory
          call get_dir
          if(len(dir) == 0)cycle

          ! Loop over groups
          do igr=1,ngr
            ! Process resonances
            call proc_res()

            ! Calculate individual properties
            ets(ikl,ijl,igr,isy) = calc_ets()
            edd(ikl,ijl,igr,isy) = calc_edd()
            qts(ikl,ijl,igr,isy) = calc_qts()
            kre(ikl,ijl,igr,isy) = calc_kre()

            ! Add contribution to krec
            krec(igr,isy) = krec(igr,isy) + dkl * kre(ikl,ijl,igr,isy)
          enddo

          ! Add contribution to other krecs
          krec(4,isy) = krec(4,isy) + dkl * kre(ikl,ijl,dchi,isy)*kex
          krec(5,isy) = krec(5,isy) + dkl * kre(ikl,ijl,1,isy)*kbfactor

          ! Extra output for level 4
          if(level == 4)call level4_extra()
        enddo

        ! Distributions of isotope effects for levels 1 and 4
        if(level == 1 .or. level == 4)then
          dze(ikl,ijl) = calc_dze()
          eta(ikl,ijl) = calc_eta()
        endif

        ! Analytic isotope effects for level 1
        if(level == 1)then
          dze_an(ikl,ijl) = calc_dze_an()
          eta_an(ikl,ijl) = calc_eta_an()
        endif

        ! Edd ratios for level 1
        if(level == 1)then
          do igr=1,ngr
            rte(ikl,ijl,igr) = calc_rte()
          enddo
        endif

        ! Q ratios for levels 1 and 4
        if(level == 1 .or. level == 4)then
          do igr=1,ngr
            rtq(ikl,ijl,igr) = calc_rtq()
          enddo
        endif
      enddo
    enddo

    ! Multiply by jlstep
    krec = krec * jlstep

    ! Sum krec over symmetries
    krec(:,3) = krec(:,1) + krec(:,2)

    ! Close resonance files
    if(level /= 1)then
      close(11)
      close(12)
      close(13)
    endif

    ! Close etra level 4 files
    if(level == 4)then
      close(20)
      close(21)
      close(22)
      close(30)
      close(31)
      close(32)
      close(40)
      close(41)
      close(42)
    endif

    ! Write jk distributions
    open(1,file=outdir//'/distr.out')
    call write_jkdistr(ets)
    call write_jkdistr(edd)
    call write_jkdistr(qts,.true.)
    call write_jkdistr(kre,.true.)
    if(level == 1 .or. level == 4)then
      call write_jkdistr_nosy_nogr(dze)
    endif
    if(level == 1)then
      call write_jkdistr_nosy_nogr(dze_an)
    endif
    if(level == 1 .or. level == 4)then
      call write_jkdistr_nosy_nogr(eta)
    endif
    if(level == 1)then
      call write_jkdistr_nosy_nogr(eta_an)
      call write_jkdistr_nosy(rte)
    endif
    if(level == 1 .or. level == 4)then
      call write_jkdistr_nosy(rtq)
    endif
    close(1)

    ! Write total results
    do isy=1,nsy+1
      ! Calculate total effects
      total_dze = krec(dclo,isy) / krec(4,isy)
      total_eta = (krec(2,isy) + krec(5,isy)) / (2 * krec(3,isy))

      ! Write
      write(OU,'(5F25.17)') krec(2,isy) * 1d35, krec(1,isy) * 1d35, krec(3,isy) * 1d35, total_dze, total_eta
    enddo
    close(OU)
  end subroutine

!-----------------------------------------------------------------------
! Processes resonances.
!-----------------------------------------------------------------------
  subroutine proc_res
    select case(level)
    ! On level 2 process all channels or process specific channel
    case(2)
      if(usech == 0)then
        do ich=1,nch
          rsf(ich) = proc_res_chan(ich)
        enddo
      else
        rsf = 0
        rsf(usech) = proc_res_chan(usech)
      endif

     ! On level 3 and 4 the lowest channel matters only
     ! On level 3 skip pathways A and B
    case(3:)
      rsf = 0
      if(level == 3 .and. igr /= 3)return
      rsf(1) = proc_res_chan(0)
    end select
  end subroutine

!-----------------------------------------------------------------------
! Processes resonances in one channel. Return resonance factor.
!-----------------------------------------------------------------------
  function proc_res_chan(ich) result(factor)
    integer ich
    real*8  factor

    ! Channel properties
    integer ichan
    real*8  baren
    real*8  barom

    ! State properties
    real*8  en
    real*8  de
    real*8  ga
    real*8  gatot
    real*8  pw
    real*8  w
    real*8  expde
    real*8  kr
    real*8  tc

    ! Extra state properties for level 4
    real*8  de4(ngr)
    real*8  ga4(2)
    real*8  pw4(2)

    ! Miscellaneous
    integer i
    integer ios
    integer ist
    real*8  thr
    character(:),allocatable :: spdir
    character(512) specstr
    character(512) resstr
    character(256) fn

    ! Setup resonance factor
    factor = 0

    ! Setup threshold
    thr = get_erc()

    ! Setup parameters depending on level
    if(level >= 3)then
      spdir = '3dsdt'
      if(level == 4)spdir = '3dsdtp'
      baren = che(1,ikl,ijl,igr,isy)
      barom = chc(1,ikl,ijl,igr,isy)
      ichan = chi(1,ikl,ijl,igr,isy)
      write(fn,'(6A)') sppath,'/',dir,'/',spdir,'/spec.out'
    else
      spdir = 'chdiag'
      baren = che(ich,ikl,ijl,igr,isy)
      barom = chc(ich,ikl,ijl,igr,isy)
      ichan = chi(ich,ikl,ijl,igr,isy)
      write(fn,'(6A,I4.4,A)') sppath,'/',dir,'/',spdir,'/spec.',ichan,'.out'
    endif

    ! Loop over states
    open(2,file=trim(fn))
    ist = 0
    do
      ! Read state properties
      read(2,'(A)',iostat=ios)specstr
      if(ios /= 0)exit

      ! Extract state number for levels 3 and 4 spectrum files
      if(level >= 3)specstr = specstr(5:)

      ! Parse spectrum line
      if(level == 4)then
        read(specstr,*)en,de4,gatot,ga4,pw4
        de = de4(igr)
        if(igr == 1)then
          ga = ga4(1)
        else
          ga = ga4(2)
        endif
        if(igr == 3)then
          pw = pw4(1)
        else
          pw = pw4(2)
        endif
      else
        read(specstr,*)en,de,gatot,pw
        ga = gatot
      endif
      ist = ist + 1

      ! Default filter
      if( en < thr .or. gatot > 300d0 )cycle

      ! Gamma filter
      if( gafilter .and. ga > max(gamod(de),10d0) )cycle

      ! Energy filter
      if( enfilter == 1 .and. de > enmax )cycle
      if( enfilter == 2 .and. en > enmax )cycle

      ! Pw filter
      if(pw < pwmin) cycle

      ! Calculate weight
      w = ga / autown * gatos1 / (gatot / autown * gatos1 + (kstab * pw + kdiss) * mconc)

      ! Calculate expde
      expde = exp( - de / ktwn )

      ! Calculate rate coefficient
      kr = (2 * jl + 1) * exp( - (baren - thr) / ktwn ) * expde * w * pw * kstab / get_qrc()

      ! kmin filter
      if(kr < kmin) cycle

      ! Calculate transmission coefficient
      tc = 1 / (1 + exp( - 2 * pi * de / barom))

      ! Write resonance
      write(resstr,'(2I5,A5,2I5,A,4F30.17)') jl,kl,sy,ichan,ist,trim(specstr),w,kr*1d35,barom,tc
      write(10+igr,'(A)')trim(resstr)

      ! Update resonance factor
      factor = factor + expde * w * pw
    enddo
    close(2)
  end function

!-----------------------------------------------------------------------
! Returns model gamma.
!-----------------------------------------------------------------------
  real*8 function gamod(de)
    real*8 de
    real*8,parameter :: a = 2d0
    gamod = sqrt(2 * de / (autown * mu)) / a * autown
  end function

!-----------------------------------------------------------------------
! Plots extra results for level 4.
!-----------------------------------------------------------------------
  subroutine level4_extra()
    ! State properties
    real*8  en
    real*8  de
    real*8  ga
    real*8  gatot

    ! Extra state properties for level 4
    real*8  de4(ngr)
    real*8  ga4(2)
    real*8  pw4(2)

    ! Miscellaneous
    integer i
    integer ios
    integer ist
    real*8  thr
    real*8  e
    real*8  g
    real*8  ea
    real*8  eb
    character(:),allocatable :: spdir
    character(512) specstr
    character(512) resstr
    character(256) fn

    ! Setup threshold
    igr = dclo
    thr = get_erc()

    ! Loop over states
    write(fn,'(4A)')sppath,'/',dir,'/3dsdtp/spec.out'
    open(2,file=trim(fn))
    ist = 0
    do
      ! Read state properties
      read(2,'(A)',iostat=ios)specstr
      if(ios /= 0)exit

      ! Extract state number
      specstr = specstr(5:)

      ! Parse spectrum line
      read(specstr,*)en,de4,gatot,ga4,pw4
      ist = ist + 1

      ! Default filter
      if( en < thr .or. gatot > 300d0 )cycle

      ! Write demo zpe
      if(pw4(2) > pwmin)then
        ! Calculate e and g
        eb = che(1,ikl,ijl,1,isy)
        ea = che(1,ikl,ijl,2,isy)
        if(mol == 686)then
          e = (en-ea) / (eb-ea)
          g = (ga4(2) - ga4(1)) / gatot
        else
          e = (en-eb) / (ea-eb)
          g = (ga4(1) - ga4(2)) / gatot
        endif

        ! Write resonance
        write(resstr,'(2I5,A5,2I5,A,2F30.17)') jl,kl,sy,-1,ist,trim(specstr),e,g
        write(20,'(A)')trim(resstr)

        ! Write resonance by symmetry
        write(resstr,'(2I5,A5,2I5,A,2F30.17)') jl,kl,sy,-1,ist,trim(specstr),e,g
        write(20+isy,'(A)')trim(resstr)
      endif

      ! Write resonances of symmetric molecule
      if(pw4(1) > pwmin)then
        ! Write resonance
        write(resstr,'(2I5,A5,2I5,A)') jl,kl,sy,-1,ist,trim(specstr)
        write(30,'(A)')trim(resstr)

        ! Write resonance by symmetry
        write(resstr,'(2I5,A5,2I5,A)') jl,kl,sy,-1,ist,trim(specstr)
        write(30+isy,'(A)')trim(resstr)
      endif

      ! Write resonances of asymmetric molecule
      if(pw4(2) > pwmin)then
        ! Write resonance
        write(resstr,'(2I5,A5,2I5,A)') jl,kl,sy,-1,ist,trim(specstr)
        write(40,'(A)')trim(resstr)

        ! Write resonance by symmetry
        write(resstr,'(2I5,A5,2I5,A)') jl,kl,sy,-1,ist,trim(specstr)
        write(40+isy,'(A)')trim(resstr)
      endif
    enddo
    close(2)
  end subroutine

!-----------------------------------------------------------------------
! Writes JK distribution to file.
!-----------------------------------------------------------------------
  subroutine write_jkdistr(a,sciarg)
    real*8  a(:,:,:,:)
    logical, optional :: sciarg
    logical :: sci

    ! Check optional argument
    if(present(sciarg))then
      sci = sciarg
    else
      sci = .false.
    endif

    ! Write JK distribution
    do ikl=nkl,1,-1
      kl = (ikl - 1) * klstep
      do isy=1,nsy
        do igr=1,ngr
          write(1,'(I30)',advance='no')kl
          do ijl=1,njl
            jl = jlmin + (ijl - 1) * jlstep
            if (kl > jl) then
              write(1,'(A30)',advance='no')'N'
            else
              if(sci)then
                write(1,'(ES30.17E3)',advance='no')a(ikl,ijl,igr,isy)
              else
                write(1,'(F30.17)',advance='no')a(ikl,ijl,igr,isy)
              endif
            endif
          enddo
          write(1,'(A30)',advance='no')'N'
        enddo
      enddo
      write(1,*)
    enddo

    ! Write J line
    do isy=1,nsy
      do igr=1,ngr
        write(1,'(A30)',advance='no')'N'
        do ijl=1,njl
          jl = jlmin + (ijl - 1) * jlstep
          write(1,'(I30)',advance='no')jl
        enddo
        write(1,'(A30)',advance='no')'N'
      enddo
    enddo
    write(1,*)
    write(1,*)
  end subroutine

!-----------------------------------------------------------------------
! Writes JK distribution to file. No symmetries.
!-----------------------------------------------------------------------
  subroutine write_jkdistr_nosy(a)
    real*8  a(:,:,:)

    ! Write JK distribution
    do ikl=nkl,1,-1
      kl = (ikl - 1) * klstep
      do igr=1,ngr
        write(1,'(I30)',advance='no')kl
        do ijl=1,njl
          jl = jlmin + (ijl - 1) * jlstep
          if (kl > jl) then
            write(1,'(A30)',advance='no')'N'
          else
            write(1,'(F30.17)',advance='no')a(ikl,ijl,igr)
          endif
        enddo
        write(1,'(A30)',advance='no')'N'
      enddo
      write(1,*)
    enddo

      ! Write J line
    do igr=1,ngr
      write(1,'(A30)',advance='no')'N'
      do ijl=1,njl
        jl = jlmin + (ijl - 1) * jlstep
        write(1,'(I30)',advance='no')jl
      enddo
      write(1,'(A30)',advance='no')'N'
    enddo
    write(1,*)
    write(1,*)
  end subroutine

!-----------------------------------------------------------------------
! Writes JK distribution to file. No symmetries. No groups.
!-----------------------------------------------------------------------
  subroutine write_jkdistr_nosy_nogr(a)
    real*8  a(:,:)

    ! Write JK distribution
    do ikl=nkl,1,-1
      kl = (ikl - 1) * klstep
      write(1,'(I30)',advance='no')kl
      do ijl=1,njl
        jl = jlmin + (ijl - 1) * jlstep
        if (kl > jl) then
          write(1,'(A30)',advance='no')'N'
        else
          write(1,'(F30.17)',advance='no')a(ikl,ijl)
        endif
      enddo
      write(1,*)
    enddo

    ! Write J line
    write(1,'(A30)',advance='no')'N'
    do ijl=1,njl
      jl = jlmin + (ijl - 1) * jlstep
      write(1,'(I30)',advance='no')jl
    enddo
    write(1,*)
    write(1,*)
  end subroutine

!-----------------------------------------------------------------------
! Returns energy of reactants.
!-----------------------------------------------------------------------
  real*8 function get_erc()
    if(igr == 1)then
      get_erc = ercb
    else
      get_erc = erca
    endif
  end function

!-----------------------------------------------------------------------
! Returns Q of reactants.
!-----------------------------------------------------------------------
  real*8 function get_qrc()
    if(igr == 1)then
      get_qrc = qrcb
    else
      get_qrc = qrca
    endif
  end function

!-----------------------------------------------------------------------
! Calculates energy of TS.
!-----------------------------------------------------------------------
  real*8 function calc_ets()
    real*8  thr
    real*8  en

    ! On level 1 return first channel above threshold
    if(level == 1)then
      ! Get threshold
      thr = get_erc()

      ! Find the first state above threshold
      do ich=1,nch
        en = che(ich,ikl,ijl,igr,isy)
        if( en > thr )then
          ichi = ich
          calc_ets = en
          exit
        endif
      enddo

      ! Stop if all are below threshold
      if(ich == nch+1) stop 'All channels are bound'

     ! On other levels return first channel
    else
      calc_ets = che(1,ikl,ijl,igr,isy)
      ichi = 1
    endif
  end function

!-----------------------------------------------------------------------
! Calculates energy double dagger.
!-----------------------------------------------------------------------
  real*8 function calc_edd()
    calc_edd = ets(ikl,ijl,igr,isy) - get_erc()
  end function

!-----------------------------------------------------------------------
! Calculates Q of transition state.
!-----------------------------------------------------------------------
  real*8 function calc_qts()
    real*8  ei
    real*8  e0
    real*8  q
    e0 = ets(ikl,ijl,igr,isy)
    q  = 0
    do ich=ichi,nch
      ei = che(ich,ikl,ijl,igr,isy)
      q  = q + exp( - (ei - e0) / ktwn ) * rsf(ich)
    enddo
    calc_qts = q
  end function

!-----------------------------------------------------------------------
! Returns reduced mass of O2 in group.
!-----------------------------------------------------------------------
  real*8 function murc()
    integer js
    if(igr == 1)then
      murc = mu0
    else
      murc = mu1
    endif
  end function

!-----------------------------------------------------------------------
! Calculates rotational energy of reactants.
!-----------------------------------------------------------------------
  real*8 function rcerot(murc,js)
    real*8  murc
    integer js
    rcerot = js * (js + 1) / (2*murc*r0**2) * autown
  end function

!-----------------------------------------------------------------------
! Calculates krec.
!-----------------------------------------------------------------------
  real*8 function calc_kre()
    ! Get J
    jl = jlmin + (ijl - 1) * jlstep

    ! Calculate rate coefficient
    calc_kre = (2 * jl + 1) * kstab * exp( - edd(ikl,ijl,igr,isy) / ktwn ) * qts(ikl,ijl,igr,isy) / get_qrc()
  end function

!-----------------------------------------------------------------------
! Calculates delta ZPE effect.
!-----------------------------------------------------------------------
  real*8 function calc_dze()
    real*8 klo
    real*8 khi

    ! Calculate terms
    klo = kre(ikl,ijl,dclo,1) + kre(ikl,ijl,dclo,2)
    khi = kre(ikl,ijl,dchi,1) + kre(ikl,ijl,dchi,2)

    ! Calculate ratio
    calc_dze = klo / ( khi * kex )
  end function

!-----------------------------------------------------------------------
! Calculates eta effect.
!-----------------------------------------------------------------------
  real*8 function calc_eta()
    real*8 ka
    real*8 kb
    real*8 ks

    ! Calculate terms
    ka = kre(ikl,ijl,2,1) + kre(ikl,ijl,2,2)
    kb = kre(ikl,ijl,1,1) + kre(ikl,ijl,1,2)
    ks = kre(ikl,ijl,3,1) + kre(ikl,ijl,3,2)

    ! Calculate ratio
    calc_eta = (ka + kb * kbfactor) / (2 * ks)
  end function

!-----------------------------------------------------------------------
! Calculates delta ZPE effect, analytic.
!-----------------------------------------------------------------------
  real*8 function calc_dze_an()
    real*8,parameter :: dzpeTS(*) = [31.95011991293968378, 32.91642837772187136 ] ! 686, 868
    real*8 mul
    real*8 ka
    real*8 kb

    ! Calculate A
    mul = m1 * (m0 + m2) / mtot
    ka  = exp( - vrot(mul,mu1) / ktha)

    ! Calculate B
    mul = m0 * (m1 + m2) / mtot
    kb  = exp( - vrot(mul,mu0) / ktha)

    ! Calculate ratio
    if(mol == 686)then
      calc_dze_an = ka / kb * exp(dzpeTS(1) / ktwn)
    else
      calc_dze_an = kb / ka * exp(dzpeTS(2) / ktwn)
    endif
  end function

!-----------------------------------------------------------------------
! Calculates eta effect, analytic.
!-----------------------------------------------------------------------
  real*8 function calc_eta_an()
    if(mol == 686)then
      calc_eta_an = ( 1 + 1 / calc_dze_an() ) / 2
    else
      calc_eta_an = ( 1 + calc_dze_an() ) / 2
    endif
  end function

!-----------------------------------------------------------------------
! Calculates rotational potential in Jacobi coordinates.
!-----------------------------------------------------------------------
  real*8 function vrot(mul,mus)
    real*8,parameter :: rl = 3.5d0
    real*8 mul
    real*8 mus

    ! Get J and K
    jl = jlmin + (ijl - 1) * jlstep
    kl = (ikl - 1) * klstep

    ! Add J component
    vrot = jl * (jl + 1) / (2 * mul * rl**2) + kl * (kl + 1) / (2 * mus * r0**2)
  end function

!-----------------------------------------------------------------------
! Calculates ratios between exp(edd).
!-----------------------------------------------------------------------
  real*8 function calc_rte()
    real*8 k1
    real*8 k2
    integer igr1
    integer igr2

    ! Set group numbers
    select case(igr)
    case(1)
      if(mol == 686)then
        igr1 = 2
        igr2 = 1
      else
        igr1 = 1
        igr2 = 2
      endif
    case(2)
      igr1 = 2
      igr2 = 3
    case(3)
      igr1 = 1
      igr2 = 3
    end select

    ! Calculate 1
    k1 = exp( - edd(ikl,ijl,igr1,1) / ktwn )  + exp( - edd(ikl,ijl,igr1,2) / ktwn )

    ! Calculate 2
    k2 = exp( - edd(ikl,ijl,igr2,1) / ktwn )  + exp( - edd(ikl,ijl,igr2,2) / ktwn )

    ! Calculate ratio
    calc_rte = k1 / k2
  end function

!-----------------------------------------------------------------------
! Calculates ratios between Qts.
!-----------------------------------------------------------------------
  real*8 function calc_rtq()
    real*8 k1
    real*8 k2
    integer igr1
    integer igr2

    ! Set group numbers
    select case(igr)
    case(1)
      if(mol == 686)then
        igr1 = 2
        igr2 = 1
      else
        igr1 = 1
        igr2 = 2
      endif
    case(2)
      igr1 = 2
      igr2 = 3
    case(3)
      igr1 = 1
      igr2 = 3
    end select

    ! Calculate 1
    k1 = qts(ikl,ijl,igr1,1) + qts(ikl,ijl,igr1,2)

    ! Calculate 2
    k2 = qts(ikl,ijl,igr2,1) + qts(ikl,ijl,igr2,2)

    ! Calculate ratio
    calc_rtq = k1 / k2
  end function
end

