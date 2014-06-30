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

  integer, parameter :: LINE_LEN = 128

  type handle
     integer :: iohandle
     character(len=LINE_LEN) :: filename
  end type handle

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
end module utility

