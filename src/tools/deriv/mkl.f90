module mkl
  implicit none

contains
  !-----------------------------------------------------------------------
  !  Derivatives calculation using MKL library (mkl.f).
  !  Not tested for odd number of points.
  !  Author: Alexander Teplukhin
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !  Real version. Initialization. 
  !  TODO: Only works for t == 2. Remove t from args?
  !-----------------------------------------------------------------------
  subroutine init_derivd(t,nk,l,freq)
    use constants
    integer nk,k,t
    
    real*8 freq(nk)
    real*8 f,l
    f=2*pi/l
    freq = 0
    do k=2,nk/2
      freq(k)=(k-1)*f
    end do
    if(t.eq.2)freq = -freq**2
  end subroutine

  !-----------------------------------------------------------------------
  !  Real version. Calculation.
  !-----------------------------------------------------------------------
  subroutine calc_derivd(t,nf,nk,freq,psi)
    use MKL_DFTI
    integer nk,nf,st,i,k,t
    real*8 psi(nk,nf),freq(nk),a,b
    type(DFTI_DESCRIPTOR),POINTER :: My_Desc1_Handle, My_Desc2_Handle

    st = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_REAL, 1, nk )
    st = DftiSetValue( My_Desc1_Handle, DFTI_NUMBER_OF_TRANSFORMS,nf)
    st = DftiSetValue( My_Desc1_Handle, DFTI_INPUT_DISTANCE, nk)
    st = DftiSetValue( My_Desc1_Handle, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT)
    st = DftiCommitDescriptor( My_Desc1_Handle )
    st = DftiComputeForward( My_Desc1_Handle, psi(:,1))
    st = DftiFreeDescriptor(My_Desc1_Handle)

    do i=1,nf
      psi(1, i) = 0
      do k=2,nk/2
        a = psi(2*(k-1),i)
        b = psi(2*(k-1)+1,i)
        if (t.eq.1) then
          psi(2*(k-1),i)   = - b * freq(k)
          psi(2*(k-1)+1,i) =   a * freq(k)
        else
          psi(2*(k-1),i)   =   a * freq(k)
          psi(2*(k-1)+1,i) =   b * freq(k)
        end if
      end do
    end do

    st = DftiCreateDescriptor(My_Desc2_Handle, DFTI_DOUBLE, DFTI_REAL, 1, nk)
    st = DftiSetValue( My_Desc2_Handle, DFTI_NUMBER_OF_TRANSFORMS,nf)
    st = DftiSetValue( My_Desc2_Handle, DFTI_INPUT_DISTANCE, nk)
    st = DftiSetValue( My_Desc2_Handle, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT)
    st = DftiSetValue( My_Desc2_Handle, DFTI_BACKWARD_SCALE, 1d0/nk)
    st = DftiCommitDescriptor(My_Desc2_Handle)
    st = DftiComputeBackward(My_Desc2_Handle, psi(:,1))
    st = DftiFreeDescriptor(My_Desc2_Handle)
  end subroutine

  !-----------------------------------------------------------------------
  !  Complex version. Initialization.
  !-----------------------------------------------------------------------
  subroutine init_derivz(nk,l,freq)
    use constants
    integer nk,k
    complex*16 f, freq(nk)
    real*8 l
    
    f=dcmplx(0,2*pi/l) 
    freq = 0
    do k=2,nk/2+1
      freq(k)=(k-1)*f
    end do
    do k=nk/2+2,nk
      freq(k)=(k-nk-1)*f
    end do
  end subroutine

  !-----------------------------------------------------------------------
  !  Complex version. Calculation.
  !-----------------------------------------------------------------------
  subroutine calc_derivz(nf,nk,freq,psi)
    use MKL_DFTI
    integer nk,nf,st,i,k
    complex*16 psi(nk,nf),freq(nk)
    type(DFTI_DESCRIPTOR),POINTER :: My_Desc1_Handle, My_Desc2_Handle

    st = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nk )
    st = DftiSetValue( My_Desc1_Handle, DFTI_NUMBER_OF_TRANSFORMS,nf)
    st = DftiSetValue( My_Desc1_Handle, DFTI_INPUT_DISTANCE, nk)
    st = DftiCommitDescriptor( My_Desc1_Handle )
    st = DftiComputeForward( My_Desc1_Handle, psi(:,1))
    st = DftiFreeDescriptor(My_Desc1_Handle)

    do i=1,nf
      do k=1,nk
        psi(k,i) = psi(k,i) * freq(k)
      end do
    end do

    st = DftiCreateDescriptor(My_Desc2_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, nk)
    st = DftiSetValue( My_Desc2_Handle, DFTI_NUMBER_OF_TRANSFORMS,nf)
    st = DftiSetValue( My_Desc2_Handle, DFTI_INPUT_DISTANCE, nk)
    st = DftiSetValue( My_Desc2_Handle, DFTI_BACKWARD_SCALE, 1d0/nk)
    st = DftiCommitDescriptor(My_Desc2_Handle)
    st = DftiComputeBackward(My_Desc2_Handle, psi(:,1))
    st = DftiFreeDescriptor(My_Desc2_Handle)
  end subroutine
end module
