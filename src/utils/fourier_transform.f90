module fourier_transform_mod
!-------------------------------------------------------------------------------------------------------------------------------------------
! Contains procedures related to Fourier Transform
!-------------------------------------------------------------------------------------------------------------------------------------------
  use constants, only : pi
  use general_utils, only : linspace
  use iso_fortran_env, only : real64

  implicit none

contains

!-------------------------------------------------------------------------------------------------------------------------------------------
! The simplest implementation of Discrete Fourier Transform (by definition).
! Any continuous periodic (or defined in a finite domain) function can be expanded (represented) in a basis of complex exponents. The 
! frequencies of these exponents are continuous too, so such basis set is infinite. The dependency of the expansion coefficiens of the 
! exponents in this expansion on their frequency is a representation of the original function in the frequency domain, 
! aka Fourier Transform (Continuous).
! If the continuous function is not known, but equidistant discrete simples of it are known, then it is possible to compute equidistant
! discrete samples of its Fourier Transform (Discrete) and this is what this procedure does.
! In other words, this procedure finds linear combination coefficients of complex exponents of specific frequencies, such that this linear
! combination passes through the given set of points.
! time_samples - equidistant discrete samples (values) of a periodic function at different times (or whatever arguments it depends on). 
! The values of argument (time) do not matter since the samples are equidistant.
! freq_samples - equidistnat frequency samples of Fourier Transform given on a grid -F/2 : F/2 - f_step : f_step for even number of points
! and -F/2 + f_step/2 : F/2 - f_step/2 : f_step for odd number of points. The frequencies in the return value of this procedure are 
! cyclically shifted (naturally) to align 0 frequency with the 1st element of freq_samples. This way the frequency
! increases in the positive index direction (freq_samples(2) corresponds to frequency = f_step, freq_samples(3) to 2*f_step, etc.) and 
! decreases in the negative index direction (freq_samples(n) corresponds to frequency = -f_step, freq_samples(n - 1) to -2*f_step, etc.).
! The second half of the array (negative frequencies) is actually redundant since the values there are complex conjugates of the values in
! the first (positive) half, but this way interpretation of the answer and further processing is simpler, so I prefer that over neglegible 
! memory efficiency.
! Here F = 2*pi/t_step and f_step = 2*pi/T, where T is the period of the original function (equivalently - range of definition), and t_step 
! is a sampling rate of the original function. This procedure does not need to know explicitly the values of t_step and T, but the 
! interpretation of the answer implicitly depdends on them. 
!-------------------------------------------------------------------------------------------------------------------------------------------
  function dft(time_samples) result(freq_samples)
    complex(real64), intent(in) :: time_samples(:)
    complex(real64) :: freq_samples(size(time_samples))
    integer :: k, n, num_points
    complex(real64) :: i

    i = cmplx(0, 1, kind = real64)
    num_points = size(time_samples)
    freq_samples = 0d0
    do k = 0, num_points - 1
      do n = 0, num_points - 1
        freq_samples(k + 1) = freq_samples(k + 1) + time_samples(n + 1) * exp(-i*2*pi*k*n / num_points)
      end do
    end do
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Transforms the frequency domain representation, obtained by discrete_fourier_transform back to the time domain
!-------------------------------------------------------------------------------------------------------------------------------------------
 function inverse_dft(freq_samples) result(time_samples)
    complex(real64), intent(in) :: freq_samples(:)
    complex(real64) :: time_samples(size(freq_samples))
    integer :: k, n, num_points
    complex(real64) :: i

    i = cmplx(0, 1, kind = real64)
    num_points = size(freq_samples)
    time_samples = 0d0
    do n = 0, num_points - 1
      do k = 0, num_points - 1
        time_samples(n + 1) = time_samples(n + 1) + freq_samples(k + 1) * exp(i*2*pi*k*n / num_points)
      end do
    end do
    time_samples = time_samples / num_points
 end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Essentially, fits *time_samples* with linear combination of complex exponents and returns derivative of the fit function of specified
! *order* at the same points. 
! *period* is the period of the given periodic function (length of definition domain)
!-------------------------------------------------------------------------------------------------------------------------------------------
  function dft_derivative(time_samples, period, order) result(deriv_time_samples)
    complex(real64), intent(in) :: time_samples(:)
    real(real64), intent(in) :: period
    integer, intent(in) :: order
    complex(real64) :: deriv_time_samples(size(time_samples))
    integer :: num_points
    real(real64) :: freq_step, freq_max
    real(real64) :: frequencies(size(time_samples))
    complex(real64) :: i
    complex(real64) :: freq_samples(size(time_samples)), deriv_freq_samples(size(time_samples))

    i = cmplx(0, 1, kind = real64)
    num_points = size(time_samples)
    freq_step = 2*pi / period
    freq_max = freq_step * num_points
    if (mod(num_points, 2) == 0) then
      frequencies = linspace(-freq_max / 2, freq_max / 2 - freq_step, num_points)
      frequencies = cshift(frequencies, num_points / 2)
    else
      frequencies = linspace(-freq_max / 2 + freq_step / 2, freq_max / 2 - freq_step / 2, num_points)
      frequencies = cshift(frequencies, (num_points - 1) / 2)
    end if

    freq_samples = dft(time_samples)
    deriv_freq_samples = freq_samples * (i*frequencies) ** order
    deriv_time_samples = inverse_dft(deriv_freq_samples)
  end function

!-------------------------------------------------------------------------------------------------------------------------------------------
! Applies dft_derivative to all columns of a given matrix
!-------------------------------------------------------------------------------------------------------------------------------------------
  function dft_derivative_all(time_samples, period, order) result(deriv_time_samples)
    complex(real64), intent(in) :: time_samples(:, :)
    real(real64), intent(in) :: period
    integer, intent(in) :: order
    complex(real64) :: deriv_time_samples(size(time_samples, 1), size(time_samples, 2))
    integer :: n

    do n = 1, size(time_samples, 2)
      deriv_time_samples(:, n) = dft_derivative(time_samples(:, n), period, order)
    end do
  end function
end module
