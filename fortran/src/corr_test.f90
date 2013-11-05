program corr_test
  use correlation
  implicit none
  real(C_DOUBLE), dimension(5) :: a, b
  integer, parameter :: maxlag = 3

  a = [1,2,3,4,5]
  b = [1,3,5,7,9]
  write(*,*) "1st"
  write(*,*) corr(a, b, maxlag)
  write(*,*) "2nd"
  write(*,*) corr(b, a, maxlag)
!  write(*,*) corr(a, b, maxlag)
!  write(*,*) corr(d, maxlag)

contains
  function nextpow2(n)
    implicit none
    integer, intent(in) :: n
    integer :: nextpow2
    nextpow2 = ceiling(log(dble(n))/log(2d0))
  end function nextpow2

end program corr_test
