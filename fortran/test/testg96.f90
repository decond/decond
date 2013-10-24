program testg96
  use g96
  implicit none
  character(len=128) :: data_filename
  integer, parameter :: natom = 18
  logical :: is_periodic
  real(8), dimension(3, natom) :: crd, vel
  real(8), dimension(3) :: cell
  real(8) :: time
  type(handle) :: data_filehandle
  integer :: stat, i
  
  is_periodic = .true.

  call get_command_argument(1, data_filename, status=stat)
  if (stat /= 0) then
    stop "error reading data filename"
  end if
  write(*,*) "data filename: ", trim(data_filename)

  write(*,*) "natom =", get_natom(data_filename)
  call read_cell(data_filename, cell, stat)
  write(*,*) "cell =", cell
  write(*,*) "status =", stat

  call open_trajectory(data_filehandle, data_filename)

  do i = 1, 2
    write(*,*) "frame", i
    call read_trajectory(data_filehandle, natom, is_periodic, crd, vel, cell, time, stat)
    if (stat > 0) then
      stop "error reading trajectory"
    else if (stat < 0) then
      stop "error eof"
    end if
    write(*,*) "crd:"
    write(*,*) crd
    write(*,*) "vel:"
    write(*,*) vel
    write(*,*) "cell:"
    write(*,*) cell
    write(*,*) "time:"
    write(*,*) time
  end do
  
  call close_trajectory(data_filehandle)
end program testg96
