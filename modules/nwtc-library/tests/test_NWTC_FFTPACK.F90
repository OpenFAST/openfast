module test_NWTC_FFTPACK

use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_FFTPACK
use NWTC_Library
use nwtc_library_test_tools

implicit none

private
public :: test_NWTC_FFTPACK_suite

real(SiKi), parameter :: tol = 1.0e-4_SiKi

contains

subroutine test_NWTC_FFTPACK_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      new_unittest("FFTPACK_real_kind", test_fftpack_real_kind), &
      new_unittest("FFT_roundtrip", test_fft_roundtrip), &
      new_unittest("FFT_forward_known", test_fft_forward_known), &
      new_unittest("FFT_forward_sign", test_fft_forward_sign), &
      new_unittest("FFT_backward_sign", test_fft_backward_sign), &
      new_unittest("FFT_cx_sign", test_fft_cx_sign), &
      new_unittest("FFT_cx_vs_real", test_fft_cx_vs_real), &
      new_unittest("CFFT_forward_known", test_cfft_forward_known), &
      new_unittest("COST_known_values", test_cost_known_values), &
      new_unittest("SINT_known_values", test_sint_known_values), &
      new_unittest("FFT_normalize", test_fft_normalize), &
      new_unittest("FFT2D_roundtrip", test_fft2d_roundtrip), &
      new_unittest("CFFT2D_roundtrip", test_cfft2d_roundtrip) &
   ]
end subroutine

! FFTPACK must be compiled with 4-byte default REAL, even in a DOUBLE_PRECISION
! build: the NWTC_FFTPACK wrapper hands it REAL(SiKi) wSave/wWork buffers sized in
! 4-byte elements.  A build that forgot to suppress default-real promotion for
! fftpack5.1.f would have FFTPACK write 8-byte elements into them and corrupt memory.
subroutine test_fftpack_real_kind(error)
   type(error_type), allocatable, intent(out) :: error
   integer, external :: FFTPACK_REALKIND

   call check(error, FFTPACK_REALKIND(), SiKi)
end subroutine

! Forward then backward (no normalization) gives x*N
subroutine test_fft_roundtrip(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 16
   real(SiKi) :: x(N), x_orig(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i

   do i = 1, N
      x(i) = sin(2.0_SiKi * Pi_D * real(i-1, SiKi) / real(N, SiKi)) &
            + 0.5_SiKi * cos(4.0_SiKi * Pi_D * real(i-1, SiKi) / real(N, SiKi))
   end do
   x_orig = x

   call InitFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT_f(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Un-normalized roundtrip: forward then backward = x * N
   do i = 1, N
      call check(error, real(x(i), kind=8), real(N, kind=8) * real(x_orig(i), kind=8), thr=real(tol, kind=8))
      if (allocated(error)) return
   end do

   call ExitFFT(fft, ErrStat)
end subroutine

! Forward FFT of known signal should produce expected coefficients
subroutine test_fft_forward_known(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 8
   real(SiKi) :: x(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat

   ! Constant signal: all 3.0 => DC = N*3 = 24, all others zero
   x = 3.0_SiKi

   call InitFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT_f(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! DC component R(1) = sum of all values = N*3 = 24
   call check(error, real(x(1), kind=8), 24.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return

   ! All other coefficients should be zero
   call check(error, real(x(2), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return
   call check(error, real(x(3), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return
   call check(error, real(x(N), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return

   call ExitFFT(fft, ErrStat)
end subroutine

! Forward FFT sign convention: FFTPACK 4.1 stores R(2k-1) = -sum h*sin(...)
subroutine test_fft_forward_sign(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 8
   real(SiKi) :: x(N), expected(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i, k
   real(SiKi) :: pi_val, arg

   pi_val = real(Pi_D, SiKi)

   ! Signal with nonzero imaginary spectral content
   do i = 1, N
      x(i) = sin(2.0_SiKi * pi_val * real(i-1, SiKi) / real(N, SiKi)) &
           + 0.5_SiKi * cos(4.0_SiKi * pi_val * real(i-1, SiKi) / real(N, SiKi)) &
           + 1.25_SiKi
   end do

   ! Reference: FFTPACK 4.1 un-normalized forward transform
   expected = 0.0_SiKi
   do i = 1, N
      expected(1) = expected(1) + x(i)
      expected(N) = expected(N) + ((-1.0_SiKi)**(i-1)) * x(i)
   end do
   do k = 2, N/2
      do i = 1, N
         arg = real(k-1, SiKi) * real(i-1, SiKi) * 2.0_SiKi * pi_val / real(N, SiKi)
         expected(2*k-2) = expected(2*k-2) + x(i) * cos(arg)
         expected(2*k-1) = expected(2*k-1) - x(i) * sin(arg)
      end do
   end do

   call InitFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT_f(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   do i = 1, N
      call check(error, real(x(i), kind=8), real(expected(i), kind=8), thr=1.0d-4)
      if (allocated(error)) return
   end do

   call ExitFFT(fft, ErrStat)
end subroutine

! Backward FFT sign convention: a pure imaginary coefficient must yield -2*sin(...)
subroutine test_fft_backward_sign(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 8
   real(SiKi) :: x(N), expected(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i
   real(SiKi) :: pi_val

   pi_val = real(Pi_D, SiKi)

   ! Spectrum with only the imaginary part of bin 1 set
   x = 0.0_SiKi
   x(3) = 1.0_SiKi

   ! FFTPACK 4.1 backward: h(i) = -2*R(3)*sin(2*pi*(i-1)/N)
   do i = 1, N
      expected(i) = -2.0_SiKi * sin(2.0_SiKi * pi_val * real(i-1, SiKi) / real(N, SiKi))
   end do

   call InitFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   do i = 1, N
      call check(error, real(x(i), kind=8), real(expected(i), kind=8), thr=1.0d-4)
      if (allocated(error)) return
   end do

   call ExitFFT(fft, ErrStat)
end subroutine

! ApplyFFT_cx sign convention: h(i) = sum 2*Re{H(k)*exp(+i*2*pi*(k-1)*(i-1)/N)}
subroutine test_fft_cx_sign(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 8
   real(SiKi) :: x(N), expected(N)
   complex(SiKi) :: H(N/2+1)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i
   real(SiKi) :: pi_val, a, b, arg

   pi_val = real(Pi_D, SiKi)
   a = 0.75_SiKi
   b = 1.5_SiKi

   ! Single complex amplitude at bin 1, endpoints must stay real
   H = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
   H(2) = CMPLX(a, b, SiKi)

   ! Inverse transform uses exp(+i...), so h(i) = 2*(a*cos - b*sin)
   do i = 1, N
      arg = 2.0_SiKi * pi_val * real(i-1, SiKi) / real(N, SiKi)
      expected(i) = 2.0_SiKi * (a * cos(arg) - b * sin(arg))
   end do

   call InitFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT_cx(x, H, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   do i = 1, N
      call check(error, real(x(i), kind=8), real(expected(i), kind=8), thr=1.0d-4)
      if (allocated(error)) return
   end do

   call ExitFFT(fft, ErrStat)
end subroutine

! ApplyFFT_cx should produce same result as ApplyFFT for equivalent input
subroutine test_fft_cx_vs_real(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 8
   real(SiKi) :: x(N), x_real(N), x_cx(N)
   complex(SiKi) :: H(N/2+1)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i

   do i = 1, N
      x(i) = sin(2.0_SiKi * Pi_D * real(i-1, SiKi) / real(N, SiKi)) + 2.0_SiKi
   end do

   call InitFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Get spectral coefficients
   x_real = x
   call ApplyFFT_f(x_real, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Build complex array from real interleaved format
   H(1) = CMPLX(x_real(1), 0.0_SiKi, SiKi)
   do i = 2, N/2
      H(i) = CMPLX(x_real(2*i-2), x_real(2*i-1), SiKi)
   end do
   H(N/2+1) = CMPLX(x_real(N), 0.0_SiKi, SiKi)

   ! Backward via complex interface
   call ApplyFFT_cx(x_cx, H, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Backward via real interface (reuse x_real which has spectral data)
   call ApplyFFT(x_real, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Both should produce identical results
   do i = 1, N
      call check(error, real(x_cx(i), kind=8), real(x_real(i), kind=8), thr=real(tol, kind=8))
      if (allocated(error)) return
   end do

   call ExitFFT(fft, ErrStat)
end subroutine

! Complex FFT forward of a pure cosine should produce energy at correct bin
subroutine test_cfft_forward_known(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 16
   complex(SiKi) :: c(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i

   ! Pure cosine at frequency bin 2: cos(2*pi*2*t/N)
   do i = 1, N
      c(i) = CMPLX(cos(4.0_SiKi * Pi_D * real(i-1, SiKi) / real(N, SiKi)), 0.0_SiKi, SiKi)
   end do

   call InitCFFT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyCFFT_f(c, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! DC should be zero
   call check(error, real(abs(c(1)), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return

   ! Bin 2 (index 3) should have energy; check it's nonzero
   call check(error, real(abs(c(3)), kind=8) > real(tol, kind=8), .true.)
   if (allocated(error)) return

   ! Bins 1 (index 2) and 3 (index 4) should be zero
   call check(error, real(abs(c(2)), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return
   call check(error, real(abs(c(4)), kind=8), 0.0d0, thr=real(tol, kind=8))
   if (allocated(error)) return

   call ExitCFFT(fft, ErrStat)
end subroutine

! Cosine transform of known input vs analytical formula
subroutine test_cost_known_values(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 9  ! must be odd
   real(SiKi) :: x(N), y(N), y_expected(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i, j
   real(SiKi) :: pi_val

   pi_val = real(Pi_D, SiKi)

   ! Input: x(j) = j for j=1..N
   do i = 1, N
      x(i) = real(i, SiKi)
   end do
   y = x

   ! Analytical formula for un-normalized DCT-I:
   ! y(J) = X(1) + (-1)^(J-1)*X(N) + sum_{K=2}^{N-1} 2*X(K)*cos((K-1)*(J-1)*pi/(N-1))
   do j = 1, N
      y_expected(j) = x(1) + ((-1.0_SiKi)**(j-1)) * x(N)
      do i = 2, N-1
         y_expected(j) = y_expected(j) + 2.0_SiKi * x(i) * &
            cos(real(i-1, SiKi) * real(j-1, SiKi) * pi_val / real(N-1, SiKi))
      end do
   end do

   call InitCOST(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyCOST(y, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   do i = 1, N
      call check(error, real(y(i), kind=8), real(y_expected(i), kind=8), thr=0.01d0)
      if (allocated(error)) return
   end do

   call ExitCOST(fft, ErrStat)
end subroutine

! Sine transform of known input vs analytical formula
subroutine test_sint_known_values(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 9  ! must be odd
   real(SiKi) :: x(N), y(N), y_expected(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i, j
   real(SiKi) :: pi_val

   pi_val = real(Pi_D, SiKi)

   ! Input: endpoints zero, interior = sin(k*pi/(N-1))
   x(1) = 0.0_SiKi
   do i = 2, N-1
      x(i) = sin(real(i-1, SiKi) * pi_val / real(N-1, SiKi))
   end do
   x(N) = 0.0_SiKi
   y = x

   ! Analytical formula for un-normalized DST-I:
   ! y(J) = sum_{K=2}^{N-1} 2*X(K)*sin((K-1)*(J-1)*pi/(N-1))
   y_expected(1) = 0.0_SiKi
   y_expected(N) = 0.0_SiKi
   do j = 2, N-1
      y_expected(j) = 0.0_SiKi
      do i = 2, N-1
         y_expected(j) = y_expected(j) + 2.0_SiKi * x(i) * &
            sin(real(i-1, SiKi) * real(j-1, SiKi) * pi_val / real(N-1, SiKi))
      end do
   end do

   call InitSINT(N, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplySINT(y, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   do i = 1, N
      call check(error, real(y(i), kind=8), real(y_expected(i), kind=8), thr=0.01d0)
      if (allocated(error)) return
   end do

   call ExitSINT(fft, ErrStat)
end subroutine

! FFT with normalization on both: forward(/N) then backward(/N) gives x/N
subroutine test_fft_normalize(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: N = 8
   real(SiKi) :: x(N), x_orig(N)
   type(FFT_DataType) :: fft
   integer :: ErrStat, i

   do i = 1, N
      x(i) = real(i, SiKi)
   end do
   x_orig = x

   call InitFFT(N, fft, NormalizeIn=.TRUE., ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT_f(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT(x, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Both forward and backward divide by N, so result = x*N / N^2 = x/N
   do i = 1, N
      call check(error, real(x(i), kind=8), real(x_orig(i) / real(N, SiKi), kind=8), thr=real(tol, kind=8))
      if (allocated(error)) return
   end do

   call ExitFFT(fft, ErrStat)
end subroutine

! 2D real FFT forward then backward (no normalization) recovers original
subroutine test_fft2d_roundtrip(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: L = 8, M = 6
   real(SiKi) :: R(L, M), R_orig(L, M)
   type(FFT2D_DataType) :: fft
   integer :: ErrStat, i, j

   do j = 1, M
      do i = 1, L
         R(i,j) = cos(2.0_SiKi * Pi_D * real(i-1, SiKi) / real(L, SiKi)) &
                + sin(2.0_SiKi * Pi_D * real(j-1, SiKi) / real(M, SiKi))
      end do
   end do
   R_orig = R

   call InitFFT2D(L, M, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT2D_f(R, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyFFT2D(R, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Forward * Backward = Identity (internally normalized)
   do j = 1, M
      do i = 1, L
         call check(error, real(R(i,j), kind=8), &
                    real(R_orig(i,j), kind=8), thr=real(tol, kind=8))
         if (allocated(error)) return
      end do
   end do

   call ExitFFT2D(fft, ErrStat)
end subroutine

! 2D complex FFT forward then backward (no normalization) recovers original
subroutine test_cfft2d_roundtrip(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: L = 8, M = 6
   complex(SiKi) :: C(L, M), C_orig(L, M)
   type(FFT2D_DataType) :: fft
   integer :: ErrStat, i, j

   do j = 1, M
      do i = 1, L
         C(i,j) = CMPLX( cos(2.0_SiKi * Pi_D * real(i-1, SiKi) / real(L, SiKi)), &
                         sin(4.0_SiKi * Pi_D * real(j-1, SiKi) / real(M, SiKi)), SiKi)
      end do
   end do
   C_orig = C

   call InitCFFT2D(L, M, fft, ErrStat=ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyCFFT2D_f(C, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   call ApplyCFFT2D(C, fft, ErrStat)
   call check(error, ErrStat, ErrID_None)
   if (allocated(error)) return

   ! Forward * Backward = Identity (internally normalized)
   do j = 1, M
      do i = 1, L
         call check(error, real(real(C(i,j)), kind=8), &
                    real(real(C_orig(i,j)), kind=8), thr=real(tol, kind=8))
         if (allocated(error)) return
         call check(error, real(aimag(C(i,j)), kind=8), &
                    real(aimag(C_orig(i,j)), kind=8), thr=real(tol, kind=8))
         if (allocated(error)) return
      end do
   end do

   call ExitCFFT2D(fft, ErrStat)
end subroutine

end module
