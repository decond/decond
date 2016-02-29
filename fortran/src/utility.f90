! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! As a special exception, you may use this file as part of a free
! software
! without restriction.  Specifically, if other files instantiate
! templates or use macros or inline functions from this file, or you
! compile
! this file and link it with other files to produce an executable, this
! file does not by itself cause the resulting executable to be covered
! by
! the GNU General Public License.  
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
! USA.
!

module utility
  use, intrinsic :: iso_c_binding
  implicit none
  private
  integer, parameter :: line_len = 1024

  public count_record_in_string, &
         get_pairindex_upper_diag, parse_version, swap, &
         get_pairindex_upper_nodiag, f2c_string

  interface swap
    module procedure swap_int, swap_char
  end interface

contains
  ! The following function is a snippet from Fortran wiki and in public
  ! domain.
  ! 
  ! This is a simple function to search for an available unit.
  ! LUN_MIN and LUN_MAX define the range of possible LUNs to check.
  ! The UNIT value is returned by the function, and also by the optional
  ! argument. This allows the function to be used directly in an OPEN
  ! statement, and optionally save the result in a local variable.
  ! If no units are available, -1 is returned.
  integer function newunit(unit)
    implicit none
    integer, intent(out), optional :: unit
    ! local
    integer, parameter :: LUN_MIN=100, LUN_MAX=200
    logical :: opened
    integer :: lun
    ! begin
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
       inquire(unit=lun,opened=opened)
       if (.not. opened) then
          newunit=lun
          exit
       end if
    end do
    if (present(unit)) unit=newunit
  end function newunit

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
