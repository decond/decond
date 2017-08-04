! This is the main decond program.
! It calls subroutines from modules mpiproc and manager.

program decond
  use mpiproc
  use varpars, only: rk
  use manager, only: init_config, read_config, prepare, decompose, output, finish
  implicit none
  real(rk) :: prog_starttime

  call mpi_setup('init')

  ! record starting time for estimating the total run time later
  prog_starttime = mpi_wtime()

  ! initialize config to default or undefined values
  call init_config()

  ! read config from the command line
  call read_config()

  ! do MPI domain decomposition, allocate memories
  ! read MD trajectories, distribute data to each node
  call prepare()

  ! perform spatial decomposition
  call decompose()

  call output()
  call finish()

  if (myrank == root) then
    write(*,*)
    write(*,*) "time for the whole program (sec):", mpi_wtime() - prog_starttime
  end if

  call mpi_setup('stop')
  stop
end program decond
