module correlation
  use fftw
  implicit none

  interface corr
    module procedure corr_vec1, corr_vec2, corr_mat
  end interface

contains
  function corr_vec1(vec, maxlag)
  ! strangely this function cannot be run on sigma, memory problem occurs!
    implicit none
    real(C_DOUBLE), dimension(2*maxlag+1) :: corr_vec1
    real(C_DOUBLE), dimension(:), intent(in) :: vec
    integer, intent(in) :: maxlag
    integer, save :: n, m, k, maxlag_save
    real(C_DOUBLE), dimension(:), allocatable, save :: vec_postpad, cor
    complex(C_DOUBLE_COMPLEX), dimension(:), allocatable, save :: vec_postpad_c, cor_c
    type(C_PTR), save :: plan1, plan2
    logical, save :: first_entry = .true.

    if (first_entry) then
      n = 0
      maxlag_save = 0
    end if

    if ((.not. allocated(vec_postpad)) .or. (n /= size(vec)) .or. (maxlag /= maxlag_save)) then
      maxlag_save = maxlag
      n = size(vec)

      ! free previous memory
      if (.not. first_entry) then
        deallocate(vec_postpad)
        deallocate(cor)
        deallocate(vec_postpad_c)
        deallocate(cor_c)
      end if

      ! memory allocation
      m = n + maxlag
      allocate(vec_postpad(m))
      allocate(cor(m))

      k = m / 2 + 1
      allocate(vec_postpad_c(k))
      allocate(cor_c(k))

      ! fftw plan
      plan1 = fftw_plan_dft_r2c_1d(m, vec_postpad, vec_postpad_c, FFTW_MEASURE)
      plan2 = fftw_plan_dft_c2r_1d(m, cor_c, cor, FFTW_MEASURE)
    endif

    ! data initialization
    vec_postpad(1:n) = vec
    vec_postpad(n+1:m) = 0d0

    call fftw_execute_dft_r2c(plan1, vec_postpad, vec_postpad_c)
    cor_c = vec_postpad_c * dconjg(vec_postpad_c)
    call fftw_execute_dft_c2r(plan2, cor_c, cor)
    corr_vec1(1:maxlag) = cor(maxlag+1:2:-1)
    corr_vec1(maxlag+1:) = cor(1:maxlag+1)
    corr_vec1 = corr_vec1 / m
    first_entry = .false.
  end function

  function corr_vec2(vec1, vec2, maxlag)
    implicit none
    real(C_DOUBLE), dimension(2*maxlag+1) :: corr_vec2
    real(C_DOUBLE), dimension(:), intent(in) :: vec1, vec2
    integer, intent(in) :: maxlag
    integer, save :: n, m, k, maxlag_save
    real(C_DOUBLE), dimension(:), allocatable, save :: vec1_prepad, vec2_postpad, cor
    complex(C_DOUBLE_COMPLEX), dimension(:), allocatable, save :: vec1_prepad_c, vec2_postpad_c, cor_c
    type(C_PTR), save :: plan1, plan2, plan3
    logical, save :: first_entry = .true.

    if (first_entry) then
      n = 0
      maxlag_save = 0
    end if

    if ((.not. allocated(vec1_prepad)) .or. (n /= size(vec1)) .or. (maxlag /= maxlag_save)) then
      maxlag_save = maxlag
      n = size(vec1)
      if (n /= size(vec2)) then
        write(*,*) "error: corr_vec2(v1, v2), v1 and v2 must have the same dimension"
        call exit(1)
      end if

      if (.not. first_entry) then
        ! free previous memory
        deallocate(vec1_prepad)
        deallocate(vec2_postpad)
        deallocate(cor)
        deallocate(vec1_prepad_c)
        deallocate(vec2_postpad_c)
        deallocate(cor_c)
      end if

      ! memory allocation
      m = n + maxlag
      allocate(vec1_prepad(m))
      allocate(vec2_postpad(m))
      allocate(cor(m))

      k = m / 2 + 1
      allocate(vec1_prepad_c(m))
      allocate(vec2_postpad_c(m))
      allocate(cor_c(m))

      ! fftw plan
      plan1 = fftw_plan_dft_r2c_1d(m, vec1_prepad, vec1_prepad_c, FFTW_MEASURE)
      plan2 = fftw_plan_dft_r2c_1d(m, vec2_postpad, vec2_postpad_c, FFTW_MEASURE)
      plan3 = fftw_plan_dft_c2r_1d(m, cor_c, cor, FFTW_MEASURE)
    endif

    ! data initialization
    vec1_prepad(1:maxlag) = 0d0
    vec1_prepad(maxlag+1:m) = vec1
    vec2_postpad(1:n) = vec2
    vec2_postpad(n+1:m) = 0d0

    call fftw_execute_dft_r2c(plan1, vec1_prepad, vec1_prepad_c)
    call fftw_execute_dft_r2c(plan2, vec2_postpad, vec2_postpad_c)

    cor_c = vec1_prepad_c * dconjg(vec2_postpad_c)

    call fftw_execute_dft_c2r(plan3, cor_c, cor)
    corr_vec2 = cor(1:2*maxlag+1) / m

    first_entry = .false.
  end function

!TODO: unfinished
  function corr_mat(mat, maxlag)
    implicit none
!    real(C_DOUBLE), dimension(maxlag+1, size(mat, 2)) :: corr_mat
    real(C_DOUBLE), dimension(:,:), allocatable :: corr_mat
    real(C_DOUBLE), dimension(:,:), intent(in) :: mat
    integer, intent(in) :: maxlag
    integer :: stat

    allocate(corr_mat(maxlag+1, size(mat, 2)**2), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: corr_mat"
      call exit(1)
    end if

    corr_mat = 5
  end function

!  function nextpow2(n)
!    implicit none
!    integer, intent(in) :: n
!    integer :: nextpow2
!    nextpow2 = ceiling(log(dble(n))/log(2d0))
!  end function nextpow2
end module correlation
