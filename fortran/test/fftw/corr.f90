module corr_mod
  use fftw
  implicit none
  private

  interface corr
    procedure corr_vec1, corr_vec2, corr_mat
  end interface
contains
  function corr_vec1(vec, maxlag)
    implicit none
    real(8), dimension(maxlag+1) :: corr_vec1
    real(8), dimension(:), intent(in) :: vec
    integer, intent(in) :: maxlag
  end function

  function corr_vec2(vec1, vec2, maxlag)
    implicit none
    real(8), dimension(maxlag+1) :: corr_vec2
    real(8), dimension(:), intent(in) :: vec1, vec2
    integer, intent(in) :: maxlag
  end function

  function corr_mat(mat, maxlag)
    implicit none
!    real(8), dimension(maxlag+1, size(mat, 2)) :: corr_mat
    real(8), dimension(:,:), allocatable :: corr_mat
    real(8), dimension(:, :), intent(in) :: mat
    integer, intent(in) :: maxlag
    integer :: stat

    allocate(corr_mat(maxlag+1, size(mat, 2)**2), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: corr_mat"
      call exit(1)
    end if
  end function

end module corr_mod
