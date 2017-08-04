! This is a wrapper module for FFTW library
! It is used only by the correlation module

module fftw
  use, intrinsic :: iso_c_binding
  private
  public c_double, c_double_complex, c_ptr, &
         fftw_measure, fftw_plan_dft_r2c_1d, &
         fftw_plan_dft_c2r_1d, fftw_execute_dft_r2c, &
         fftw_execute_dft_c2r
  include 'fftw3.f03'
end module

