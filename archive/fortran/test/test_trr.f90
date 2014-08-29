program test_trr
  use trr
  implicit none
  type(handle) :: dataFileHandle
  character(len=128) :: dataFilename
  integer :: status
  integer, parameter :: natom_to_read = 36
  logical :: is_periodic
  real(8) :: crd(3, natom_to_read), vel(3, natom_to_read)
  real(8) :: cell(3), time

  is_periodic = .true.

  call get_command_argument(1, dataFilename)
  write(*,*) "data.trr = ", dataFilename

  call open_trajectory(dataFileHandle, dataFilename)
  call read_trajectory(dataFileHandle, natom_to_read, is_periodic, crd, vel, cell, time, status)
  call read_trajectory(dataFileHandle, natom_to_read, is_periodic, crd, vel, cell, time, status)
  call close_trajectory(dataFileHandle)

  write(*,*) "crd"
  write(*,*) crd
  write(*,*) "********************"
  write(*,*) "vel"
  write(*,*) vel
end program test_trr
