program corr_test
  use correlation
  implicit none
  real(C_DOUBLE), dimension(5) :: a, b
  integer, parameter :: maxlag = 3

  a = [1,2,3,4,5]
  b = [1,3,5,7,9]
!  write(*,*) "Run 1"
!  write(*,*) corr(a, b, maxlag)
!  write(*,*) "Run 2"
!  write(*,*) corr(b, a, maxlag)
  write(*,*) "Run 3"
  write(*,*) corr(a, maxlag)
  write(*,*) "Run 4"
  write(*,*) corr(a, a, maxlag)

contains
  function nextpow2(n)
    implicit none
    integer, intent(in) :: n
    integer :: nextpow2
    nextpow2 = ceiling(log(dble(n))/log(2d0))
  end function nextpow2

end program corr_test
