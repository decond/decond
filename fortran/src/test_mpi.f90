program test_mpi
  use mpi
  implicit none
  integer :: totNumAtom
  integer :: numFrame, stat
  character(len=10) :: tmp_str
  !one frame data (dim=3, atom) 
  real(8), allocatable :: pos(:, :, :), vel(:, :, :)

  !MPI variables
  integer :: ierr, nprocs, myrank
  integer, parameter :: root = 0
  real(8) :: starttime, endtime

  !initialize
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)

  !welcome message
!  write(*,*) "hello from rank ", myrank

write(*,*) "rank=",myrank, ", tag0"
  !read parameters for all ranks
  call get_command_argument(1, tmp_str) ! in the unit of frame number
  read(tmp_str, *) numFrame 

  call get_command_argument(2, tmp_str) ! in the unit of frame number
  read(tmp_str, *) totNumAtom

  !rank root output parameters read
  if (myrank == root) then
    write(*,*) "numFrame= ", numFrame
    write(*,*) "totNumAtom= ", totNumAtom
  end if

write(*,*) "rank=",myrank, ", tag1"
  !prepare memory for all ranks
  allocate(pos(3, numFrame, totNumAtom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: pos"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  allocate(vel(3, numFrame, totNumAtom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  if (myrank == root) then
    pos=1d0
    vel=2d0
  else
    pos=0d0
    vel=0d0
  end if

write(*,*) "Before broadcasting trajectory"
write(*,*) "rank=",myrank, ", pos=", pos(1,1,1)
write(*,*) "rank=",myrank, ", vel=", vel(1,1,1)
  !distribute trajectory data collectively
  if (myrank == root) write(*,*) "Start broadcasting trajectory"
  starttime = MPI_Wtime()
  call mpi_bcast(pos, 3*numFrame*totNumAtom, mpi_double_precision, root, MPI_COMM_WORLD, ierr)
  call mpi_bcast(vel, 3*numFrame*totNumAtom, mpi_double_precision, root, MPI_COMM_WORLD, ierr)
write(*,*) "rank=",myrank, ", pos=", pos(1,1,1)
write(*,*) "rank=",myrank, ", vel=", vel(1,1,1)
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished broadcasting trajectory. It took ", endtime - starttime, " seconds"

  call mpi_finalize(ierr)
  stop
end program test_mpi
