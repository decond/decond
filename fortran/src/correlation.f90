module correlation
  use fftw
  implicit none

  interface corr
    procedure corr_vec1, corr_vec2, corr_mat
  end interface

contains
  function corr_vec1(vec, maxlag)
  ! strangely this function cannot be run on sigma, memory problem occurs!
    implicit none
    real(C_DOUBLE), dimension(2*maxlag+1) :: corr_vec1
    real(C_DOUBLE), dimension(:), intent(in) :: vec
    integer, intent(in) :: maxlag
    integer, save :: n, m, k, maxlag_save
    real(C_DOUBLE), dimension(:), pointer, save :: vec_postpad, cor
    complex(C_DOUBLE_COMPLEX), dimension(:), pointer, save :: vec_postpad_c, cor_c
    type(C_PTR), save :: p, pcor, p_c, pcor_c
    type(C_PTR), save :: plan1, plan2
    logical, save :: first_entry = .true.

    if (first_entry) then
      nullify(vec_postpad)
      n = 0
      maxlag_save = 0
      first_entry = .false.
    end if

    if ((.not. associated(vec_postpad)) .or. (n /= size(vec)) .or. (maxlag /= maxlag_save)) then
      maxlag_save = maxlag
      n = size(vec)

      ! free previous memory
      call fftw_free(p)
      call fftw_free(pcor)
      call fftw_free(p_c)
      call fftw_free(pcor_c)

      ! memory allocation
      m = n + maxlag
      p = fftw_alloc_real(int(m, C_SIZE_T))
      pcor = fftw_alloc_real(int(m, C_SIZE_T))
      call c_f_pointer(p, vec_postpad, [m])
      call c_f_pointer(pcor, cor, [m])

      k = m / 2 + 1
      p_c = fftw_alloc_complex(int(k, C_SIZE_T))
      pcor_c = fftw_alloc_complex(int(k, C_SIZE_T))
      call c_f_pointer(p_c, vec_postpad_c, [k])
      call c_f_pointer(pcor_c, cor_c, [k])

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
  end function

  function corr_vec2(vec1, vec2, maxlag)
    implicit none
    real(C_DOUBLE), dimension(2*maxlag+1) :: corr_vec2
    real(C_DOUBLE), dimension(:), intent(in) :: vec1, vec2
    integer, intent(in) :: maxlag
    integer, save :: n, m, k, maxlag_save
    real(C_DOUBLE), dimension(:), pointer, save :: vec1_prepad, vec2_postpad, cor
    complex(C_DOUBLE_COMPLEX), dimension(:), pointer, save :: vec1_prepad_c, vec2_postpad_c, cor_c
    type(C_PTR), save :: p1, p2, pcor, p1_c, p2_c, pcor_c
    type(C_PTR), save :: plan1, plan2, plan3
    logical, save :: first_entry = .true.

    if (first_entry) then
      nullify(vec1_prepad)
      n = 0
      maxlag_save = 0
      first_entry = .false.
    end if

    if ((.not. associated(vec1_prepad)) .or. (n /= size(vec1)) .or. (maxlag /= maxlag_save)) then
      maxlag_save = maxlag
      n = size(vec1)
      if (n /= size(vec2)) then
        write(*,*) "error: corr_vec2(v1, v2), v1 and v2 must have the same dimension"
        call exit(1)
      end if

!      if (.not. first_entry) then
      ! free previous memory
      call fftw_free(p1)
      call fftw_free(p2)
      call fftw_free(pcor)
      call fftw_free(p1_c)
      call fftw_free(p2_c)
      call fftw_free(pcor_c)
!      end if

      ! memory allocation
      m = n + maxlag
      p1 = fftw_alloc_real(int(m, C_SIZE_T))
      p2 = fftw_alloc_real(int(m, C_SIZE_T))
      pcor = fftw_alloc_real(int(m, C_SIZE_T))
      call c_f_pointer(p1, vec1_prepad, [m])
      call c_f_pointer(p2, vec2_postpad, [m])
      call c_f_pointer(pcor, cor, [m])

      k = m / 2 + 1
      p1_c = fftw_alloc_complex(int(k, C_SIZE_T))
      p2_c = fftw_alloc_complex(int(k, C_SIZE_T))
      pcor_c = fftw_alloc_complex(int(k, C_SIZE_T))
      call c_f_pointer(p1_c, vec1_prepad_c, [k])
      call c_f_pointer(p2_c, vec2_postpad_c, [k])
      call c_f_pointer(pcor_c, cor_c, [k])

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
  end function

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
