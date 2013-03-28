! -*- F90 -*-
! ERMOD - Eneregy Representation Module
! Copyright (C) 2000-2012 Nobuyuki Matubayasi
! Copyright (C) 2010-2012 Shun Sakuraba
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
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

! module for DCD I/O, as an example of fortran I/O module
module trajectory
  type handle
     integer :: iohandle
     logical :: have_cell_info
  end type handle
  
contains

  ! Open trajectory and returns handle as htraj. 
  ! Should open fails, the program abends.
  subroutine open_trajectory(htraj, fname)
    use utility, only: newunit
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname
    integer(4) :: dcd_header(21)
    integer(4), parameter :: dcd_magic_little = X'44524f43', dcd_magic_big = X'434f5244'

    open(unit=newunit(htraj%iohandle), file=fname, action="READ", form="UNFORMATTED")
    ! Check dcd magic
    read(htraj%iohandle) dcd_header(:)
    if(dcd_header(1) /= dcd_magic_little .and. dcd_header(1) /= dcd_magic_big) then
       stop "incorrect dcd format"
    end if
    htraj%have_cell_info = (dcd_header(12) == 1)
    read(htraj%iohandle)
    read(htraj%iohandle)
  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    implicit none
    type(handle), intent(inout) :: htraj

    close(htraj%iohandle)
  end subroutine close_trajectory

  ! Read trajectory and returns [crd] as a coordinates, and [cell] as a
  ! periodic cell, represented in Angstrom.
  ! [status] is non-zero if any error occurs. In such a case, [crd] and
  ! [cell] can be an arbitrary value.
  ! [cell] may be an arbitrary value if the trajectory does not contain
  ! cell information.
  ! The coordinate is not guaranteed to be within a unit cell.
  subroutine read_trajectory(htraj, natom, is_periodic, crd, cell, status)
    implicit none
    type(handle), intent(in) :: htraj
    integer, intent(in) :: natom
    logical, intent(in) :: is_periodic
    real(8), intent(out) :: crd(3,natom)
    real(8), intent(out) :: cell(3,3)
    integer, intent(out) :: status
    real(4), allocatable :: buffer(:)

    cell(:, :) = 0.
    allocate(buffer(natom))

    if(is_periodic) then
       if(.not. htraj%have_cell_info) stop "Cell info is requested, but does not exist!"
       read(htraj%iohandle, err=999, end=999) cell(1,1), cell(1,2), cell(2,2), cell(1,3), cell(2,3), cell(3,3)
    end if

    read(htraj%iohandle, err=999, end=999) buffer(:)
    crd(1, :) = buffer(:)
    read(htraj%iohandle, err=999, end=999) buffer(:)
    crd(2, :) = buffer(:)
    read(htraj%iohandle, err=999, end=999) buffer(:)
    crd(3, :) = buffer(:)

    status = 0
    deallocate(buffer)
    return

999 status = 1
    deallocate(buffer)
    return
  end subroutine read_trajectory

end module trajectory

