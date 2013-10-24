module g96
  use utility, only : handle, newunit
  character(len=*), parameter :: key_title = "TITLE"
  character(len=*), parameter :: key_timestep = "TIMESTEP"
  character(len=*), parameter :: key_positionred = "POSITIONRED"
  character(len=*), parameter :: key_position = "POSITION"
  character(len=*), parameter :: key_velocityred = "VELOCITYRED"
  character(len=*), parameter :: key_velocity = "VELOCITY"
  character(len=*), parameter :: key_box = "BOX"
  character(len=*), parameter :: key_end = "END"
  integer, parameter :: line_len = 128

contains
  ! Open trajectory and returns handle as htraj. 
  ! Should open fails, the program abends.
  subroutine open_trajectory(htraj, fname)
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname
    character(len=line_len) :: line

    ! use form="FORMATTED" for human-redable file format
    open(unit=newunit(htraj%iohandle), file=fname, action="READ", form="FORMATTED")
    read(htraj%iohandle, *) line
    if (line /= key_title) then
      stop "incorrect g96 format: no " // key_title
    end if
    read(htraj%iohandle, *) line
    do while (line /= key_end)
      read(htraj%iohandle, *) line
    end do
    
    ! Read header info, skip header, etc.
  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    implicit none
    type(handle), intent(inout) :: htraj

    close(htraj%iohandle)

    ! release extra resources if necessary
  end subroutine close_trajectory

  ! Read trajectory and returns [crd] as a coordinates, and [cell] as a
  ! periodic cell, represented in nanometer.
  ! [status] is non-zero if any error occurs. In such a case, [crd] and
  ! [cell] may be an arbitrary value if the trajectory does not contain
  ! cell information.
  ! The coordinate is not guaranteed to be within a unit cell.
  subroutine read_trajectory(htraj, natom, is_periodic, crd, vel, cell, time, status)
    implicit none
    type(handle), intent(in) :: htraj
    integer, intent(in) :: natom
    logical, intent(in) :: is_periodic
    real(8), intent(out) :: crd(3,natom), vel (3,natom)
    real(8), intent(out) :: cell(3)
    real(8), intent(out) :: time
    integer, intent(out) :: status
    character(len=line_len) :: line
    integer :: i, dummy_i
    
    read(htraj%iohandle, *, err = 901, end = 900) line
    if (line /= key_timestep) then
      stop "incorrect g96 format: no " // key_timestep
    end if
    read(htraj%iohandle, *, err = 901, end = 900) dummy_i, time
    read(htraj%iohandle, *, err = 901, end = 900) !END

    read(htraj%iohandle, *, err = 901, end = 900) line
    if (line /= key_positionred .and. line /= key_position) then
      stop "incorrect g96 format: no " // key_positionred // " nor " // key_position
    end if
    do i = 1, natom
      read(htraj%iohandle, "(3F15.9)", err = 901, end = 900) crd(:, i)
    end do
    read(htraj%iohandle, *, err = 901, end = 900) !END

    read(htraj%iohandle, *, err = 901, end = 900) line
    if (line /= key_velocityred .and. line /= key_velocity) then
      stop "Currently velocity information is mandatory: " // key_velocityred // " or " // key_velocity
    end if
    do i = 1, natom
      read(htraj%iohandle, "(3F15.9)", err = 901, end = 900) vel(:, i)
    end do
    read(htraj%iohandle, *, err = 901, end = 900) !END

    if (is_periodic) then
      read(htraj%iohandle, *, err = 901, end = 900) line
      if (line == key_box) then
        read(htraj%iohandle, "(3F15.9)", err = 901, end = 900) cell
        read(htraj%iohandle, *, err = 901, end = 900) !END
      else
        stop "Trajectory with PBC should contain information of cell size"
      end if
    end if

    status = 0
    return

900 status = -1
    return
    
901 status = 1
    return
  end subroutine read_trajectory
  
  function get_natom(fname)
    implicit none
    integer :: get_natom
    character(len=*), intent(in):: fname
    type(handle) :: htraj
    character(len=line_len) :: line

    get_natom = 0
    call open_trajectory(htraj, fname)

    read(htraj%iohandle, *, err = 902, end = 902) line
    do while (line /= key_position .and. line /= key_positionred)
      read(htraj%iohandle, *, err = 902, end = 902) line
    end do

    read(htraj%iohandle, *, err = 902, end = 902) line
    do while (line /= key_end)
      get_natom = get_natom + 1
      read(htraj%iohandle, *, err = 902, end = 902) line
    end do

    call close_trajectory(htraj)
    return

902 get_natom = -1
    call close_trajectory(htraj)
    return
  end function get_natom

  subroutine read_cell(fname, cell, status)
    implicit none
    character(len=*), intent(in):: fname
    real(8), dimension(3), intent(out) :: cell
    integer, intent(out) :: status
    type(handle) :: htraj
    character(len=line_len) :: line

    call open_trajectory(htraj, fname)

    read(htraj%iohandle, *, err = 903, end = 903) line
    do while (line /= key_box)
      read(htraj%iohandle, *, err = 903, end = 903) line
    end do
    read(htraj%iohandle, *, err = 903, end = 903) cell

    status = 0
    call close_trajectory(htraj)
    return

903 status = -1
    call close_trajectory(htraj)
    return
  end subroutine read_cell
end module g96

