program fftw_test
  use fftw3
  implicit none
  real(C_DOUBLE), dimension(16) :: a, b, c
  complex(C_DOUBLE_COMPLEX), dimension(9) :: aa, bb, cc
  type(C_PTR) :: plan1, plan2, plan3
  integer, parameter :: maxlag = 5
  integer :: m

  plan1 = fftw_plan_dft_r2c_1d(16, a, aa, FFTW_MEASURE)
  plan2 = fftw_plan_dft_r2c_1d(16, b, bb, FFTW_MEASURE)
  plan3 = fftw_plan_dft_c2r_1d(16, cc, c, FFTW_MEASURE)

  a = [0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0]
  b = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 0, 0, 0, 0, 0]
  
!  m = 2**nextpow2(size(a))
  
  call fftw_execute_dft_r2c(plan1, a, aa)
  call fftw_execute_dft_r2c(plan2, b, bb)

  cc = aa * dconjg(bb)

  call fftw_execute_dft_c2r(plan3, cc, c)

  write(*,*) c(1:2*maxlag+1)/16


contains
  function nextpow2(n)
    implicit none
    integer, intent(in) :: n
    integer :: nextpow2
    nextpow2 = ceiling(log(dble(n))/log(2d0))
  end function nextpow2

end program fftw_test
