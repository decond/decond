program decond
  use mpiproc
  use manager, only: is_pa_mode, is_pm_mode, init_config, read_config, &
                     prepare, decompose, output, finish
  implicit none

  call mpi_setup('init')

  prog_starttime = MPI_Wtime()

  call init_config()

  call read_config()

  call prepare()

  call decompose()

  call output()

  if (myrank == root) write(*,*)
  if (myrank == root) write(*,*) "time for the whole program (sec):", MPI_Wtime() - prog_starttime

  call finish()
  call mpi_setup('stop')
  stop
end program decond
