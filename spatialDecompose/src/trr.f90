module trr
  use utility, only : handle, newunit
  integer, parameter :: line_len = 128

contains
  ! Open trajectory and returns handle as htraj. 
  ! Should open fails, the program abends.
  subroutine open_trajectory(htraj, fname)
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname

    ! use form="FORMATTED" for human-redable file format
    open(unit=newunit(htraj%iohandle), file=fname, action="READ", &
        &form="UNFORMATTED", access="STREAM", convert="BIG_ENDIAN")
!    rewind(htraj%iohandle) !not sure if it is necessary
    
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
  subroutine read_trajectory(htraj, natom_to_read, is_periodic, crd, vel, cell, time, status)
    implicit none
    type(handle), intent(in) :: htraj
    integer, intent(in) :: natom_to_read
    logical, intent(in) :: is_periodic
    real(8), intent(out) :: crd(3,natom_to_read), vel (3,natom_to_read)
    real(8), intent(out) :: cell(3)
    real(8), intent(out) :: time
    integer, intent(out) :: status
    integer :: i
    character(len=52) :: intro_words
    integer(4) :: data_present(3) !coordinates, velocity, force
    logical :: is_data_present(3) !coordinates, velocity, force
    integer(4) :: natom, frame_step(2)
    real(8) :: box_params(10), dummy_data(3)

    read(htraj%iohandle) intro_words
!    write(*,*) intro_words(13:24)
    read(htraj%iohandle) data_present
    is_data_present = (data_present /= 0)
!    write(*,*) "data_present (x v f) = (", is_data_present, ")"
    read(htraj%iohandle) natom
!    write(*,*) "natom = ", natom
    read(htraj%iohandle) frame_step
!    write(*,*) "frame_step = ", frame_step
    read(htraj%iohandle) time
!    write(*,*) "frame_time = ", time
    read(htraj%iohandle) box_params
    cell(1:3) = box_params(2:10:4)
!    write(*,*) "cell = ", cell

    ! read coordinates
    if (is_data_present(1)) then
      read(htraj%iohandle, err = 901, end = 900) crd
      do i = 1, natom - natom_to_read
        read(htraj%iohandle, err = 901, end = 900) dummy_data
      end do
    end if

    ! read velocities
    if (is_data_present(2)) then
      read(htraj%iohandle, err = 901, end = 900) vel
      do i = 1, natom - natom_to_read
        read(htraj%iohandle, err = 901, end = 900) dummy_data
      end do
    end if

    status = 0
    return

900 status = -1
    return
    
901 status = 1
    return
  end subroutine read_trajectory
  
!  function get_natom(fname)
!    implicit none
!    integer :: get_natom
!    character(len=*), intent(in):: fname
!    type(handle) :: htraj
!    character(len=line_len) :: line
!
!    get_natom = 0
!    call open_trajectory(htraj, fname)
!
!    read(htraj%iohandle, *, err = 902, end = 902) line
!    do while (line /= key_position .and. line /= key_positionred)
!      read(htraj%iohandle, *, err = 902, end = 902) line
!    end do
!
!    read(htraj%iohandle, *, err = 902, end = 902) line
!    do while (line /= key_end)
!      get_natom = get_natom + 1
!      read(htraj%iohandle, *, err = 902, end = 902) line
!    end do
!
!    call close_trajectory(htraj)
!    return
!
!902 get_natom = -1
!    call close_trajectory(htraj)
!    return
!  end function get_natom

!  subroutine read_cell(fname, cell, status)
!    implicit none
!    character(len=*), intent(in):: fname
!    real(8), dimension(3), intent(out) :: cell
!    integer, intent(out) :: status
!    type(handle) :: htraj
!    character(len=line_len) :: line
!
!    call open_trajectory(htraj, fname)
!
!    read(htraj%iohandle, *, err = 903, end = 903) line
!    do while (line /= key_box)
!      read(htraj%iohandle, *, err = 903, end = 903) line
!    end do
!    read(htraj%iohandle, *, err = 903, end = 903) cell
!
!    status = 0
!    call close_trajectory(htraj)
!    return
!
!903 status = -1
!    call close_trajectory(htraj)
!    return
!  end subroutine read_cell
end module trr

