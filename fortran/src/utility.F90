! This module contains small helper subroutines

module utility
  use, intrinsic :: iso_c_binding
  implicit none
  private
  integer, parameter :: line_len = 1024

  public count_record_in_string, &
         get_pairindex_upper_diag, parse_version, swap, &
         get_pairindex_upper_nodiag, f2c_string

#ifndef FORTRAN2008
  public newunit
#endif

  interface swap
    module procedure swap_int, swap_char
  end interface

contains
#ifndef FORTRAN2008
  integer function newunit(unit) result(n)
    ! http://www.fortran90.org/src/best-practices.html
    ! returns lowest i/o unit number not in use
    integer, intent(out), optional :: unit
    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent
    do n = nmin, nmax
      inquire(unit=n, opened=inuse)
      if (.not. inuse) then
        if (present(unit)) unit=n
        return
      end if
    end do
    stop "newunit ERROR: available unit not found."
  end function
#endif /* FORTRAN2008 */

  function f2c_string(string)
    implicit none
    character(len=*) :: string
    character(len=1, kind=C_CHAR), dimension(len_trim(string)+1) :: f2c_string
    integer :: i
    
    do i = 1, len_trim(string)+1
      f2c_string(i) = string(i:i)
    end do
    f2c_string(len_trim(string)+1) = C_NULL_CHAR
  end function

  integer function count_record_in_string(string)
    implicit none
    character(len=*) :: string
    character :: prev_c
    integer :: cur
    logical :: is_prev_record_ended_by_comma

    count_record_in_string = 0
    string = adjustl(string)
    cur = 0
    is_prev_record_ended_by_comma = .false.
    prev_c = ' '

    do cur = 1, len_trim(string)
      select case (string(cur:cur))
        case (',')
          if (prev_c /= ' ' .or. cur == 1 .or. is_prev_record_ended_by_comma) then
            count_record_in_string = count_record_in_string + 1
            is_prev_record_ended_by_comma = .true.
          end if
        case (' ')
          if (prev_c /= ' ' .and. prev_c /= ',') then
            count_record_in_string = count_record_in_string + 1
            is_prev_record_ended_by_comma = .false.
          end if
      end select
      prev_c = string(cur:cur)
    end do
    if (prev_c /= ' ' .and. prev_c /= ',') then
      count_record_in_string = count_record_in_string + 1
    end if
  end function

  integer function get_pairindex_upper_diag(i, j, n)
    implicit none
    integer, intent(in) :: i, j, n
    integer :: r, c
    !          c
    !    | 1  2  3  4
    !  --+------------
    !  1 | 1  2  3  4
    !    |
    !  2 |    5  6  7
    !r   |
    !  3 |       8  9
    !    |
    !  4 |         10
    !
    !  index(r, c) = (r - 1) * n + c - r * (r - 1) / 2
    !  where n = size(c) = size(r), r <= c

    r = min(i, j)
    c = max(i, j)
    if (r < 1 .or. c < 1 .or. n < 1 .or. r > n .or. c > n) then
      get_pairindex_upper_diag = -1
    else
      get_pairindex_upper_diag = (r - 1) * n + c - r * (r - 1) / 2
    end if
  end function get_pairindex_upper_diag

  integer function get_pairindex_upper_nodiag(i, j, n)
    implicit none
    integer, intent(in) :: i, j, n
    integer :: r, c
    !          c
    !    | 1  2  3  4
    !  --+------------
    !  1 |    1  2  3
    !    |
    !  2 |       4  5
    !r   |
    !  3 |          6
    !    |
    !  4 |
    !
    !  index(r, c) = (r - 1) * n + c - ((r + 1) * r) / 2
    !  where n = size(c) = size(r), r < c

    r = min(i, j)
    c = max(i, j)
    if (r < 1 .or. c < 1 .or. n < 1 .or. r == c .or. r > n .or. c > n) then
      get_pairindex_upper_nodiag = -1
    else
      get_pairindex_upper_nodiag = (r - 1) * n + c - ((r + 1) * r) / 2
    end if
  end function get_pairindex_upper_nodiag

  subroutine parse_version(ver, major, minor, patch)
    implicit none
    character(len=11), intent(in) :: ver
    integer, intent(out) :: major
    integer, optional, intent(out) :: minor, patch
    integer :: p1, p2, p_null

    p1 = scan(ver, '.')
    p2 = scan(ver, '.', .true.)
    p_null = scan(ver, char(0))
    read(ver(1:p1-1), *) major
    if (present(minor)) read(ver(p1+1:p2-1), *) minor
    if (present(patch)) read(ver(p2+1:p_null-1), *) patch
  end subroutine parse_version

  subroutine swap_int(list, i, j)
    implicit none
    integer, intent(inout) :: list(:)
    integer, intent(in) :: i, j
    integer :: tmp

    tmp = list(i)
    list(i) = list(j)
    list(j) = tmp
  end subroutine swap_int

  subroutine swap_char(list, i, j)
    implicit none
    character(len=line_len), intent(inout) :: list(:)
    integer, intent(in) :: i, j
    character(len=line_len) :: tmp

    tmp = list(i)
    list(i) = list(j)
    list(j) = tmp
  end subroutine swap_char
end module utility
