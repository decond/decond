program test_xdr
  use utility, only : handle
  use xdr
  implicit none
  integer :: natoms
  character(len=128) :: filename, filename2
  type(handle) :: filehandle, filehandle2
  logical :: is_periodic

  real(8), allocatable :: crd(:, :), vel(:, :)
  real(8) :: cell(3), time
  integer :: status

  is_periodic = .true.
  filename = "../test/test.trr"
  filename2 = "../test/test2.trr"
  
  natoms = get_natom(filename)
  write(*,*) "natoms=", natoms

  allocate(crd(3,natoms))
  allocate(vel(3,natoms))

  call open_trajectory(filehandle, filename)
  write(*,*) "handle=", filehandle%iohandle
  call close_trajectory(filehandle)
  write(*,*) "handle closed"

!  call read_trajectory(filehandle, natoms, is_periodic, crd, vel, cell, time, status)
!  write(*,*) "crd=", crd
!  write(*,*) "vel=", vel
!  write(*,*) "cell=", cell
!  write(*,*) "time=", time
!  write(*,*) "status=", status


  call open_trajectory(filehandle2, filename2)
  write(*,*) "handle2=", filehandle2%iohandle
  call close_trajectory(filehandle2)
  write(*,*) "handle2 closed"

write(*,*) "--------------"

  call open_trajectory(filehandle, filename)
  write(*,*) "handle=", filehandle%iohandle
  call open_trajectory(filehandle2, filename2)
  write(*,*) "handle2=", filehandle2%iohandle
  call close_trajectory(filehandle)
  write(*,*) "handle closed"
  call close_trajectory(filehandle2)
  write(*,*) "handle2 closed"
end program test_xdr
