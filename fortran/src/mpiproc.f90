module mpiproc
  ! MPI
  use mpi
  implicit none

  ! mpi common variables
  integer :: ierr
  integer :: myrank, nprocs       ! rank number and total number of processes
  integer, parameter :: root = 0

  integer:: numDomain_r, numDomain_c, numMolPerDomain_r, numMolPerDomain_c
  integer :: r_start, r_end, c_start, c_end
  real(8) :: starttime, endtime, starttime2, prog_starttime
  integer :: r_start_offset, c_start_offset
  integer :: residueMol_r, residueMol_c, num_r, num_c
  integer :: row_comm, col_comm, r_group_idx, c_group_idx, offset
  integer, dimension(:), allocatable :: displs_r, displs_c, scounts_r, scounts_c

  real(8) :: dummy_null

contains
  subroutine mpi_setup(type)
    implicit none
    character(len=4) :: type

    if(type == 'init') then
      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
    else if(type == 'stop') then
      call mpi_finalize(ierr)
    endif
  end subroutine mpi_setup

  subroutine mpi_abend()
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end subroutine mpi_abend

  subroutine halt_with_message(message, rank)
    implicit none
    character(len=*), intent(in) :: message
    integer, optional, intent(in) :: rank
    if (present(rank)) then
      if (myrank == rank) write(*, *) message
    else
      write(*, *) message
    end if
    call mpi_abend()
  end subroutine halt_with_message

  subroutine domainDecomposition(totNumMol, numFrame)
    !domain decomposition for atom pairs (numDomain_r * numDomain_c = nprocs)
    !numMolPerDomain_r * numDomain_r ~= totNumMol
    implicit none
    integer, intent(in) :: totNumMol, numFrame
    integer :: i, stat

    if (numDomain_r == 0 .and. numDomain_c == 0) then
      numDomain_c = nint(sqrt(dble(nprocs)))
      do while(mod(nprocs, numDomain_c) /= 0)
        numDomain_c = numDomain_c - 1
      end do
      numDomain_r = nprocs / numDomain_c
    else if (numDomain_r > 0) then
      numDomain_c = nprocs / numDomain_r
    else if (numDomain_c > 0) then
      numDomain_c = nprocs / numDomain_r
    else
      write(*,*) "Invalid domain decomposition: ", numDomain_r, " x ", numDomain_c
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    if (numDomain_r * numDomain_c /= nprocs) then
      write(*,*) "Domain decomposition failed: ", numDomain_r, " x ", numDomain_c, " /= ", nprocs
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    !Determine row and column position for the node
    r_group_idx = mod(myrank, numDomain_r) !column-major mapping
    c_group_idx = myrank / numDomain_r

    !Split comm into row and column comms
    call mpi_comm_split(MPI_COMM_WORLD, c_group_idx, r_group_idx, col_comm, ierr)
    !color by row, rank by column
    call mpi_comm_split(MPI_COMM_WORLD, r_group_idx, c_group_idx, row_comm, ierr)
    !color by column, rank by row

    numMolPerDomain_r = totNumMol / numDomain_r
    numMolPerDomain_c = totNumMol / numDomain_c
    residueMol_r = mod(totNumMol, numDomain_r)
    residueMol_c = mod(totNumMol, numDomain_c)

    allocate(displs_r(numDomain_r), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: displs_r"
      call exit(1)
    end if
    allocate(displs_c(numDomain_c), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: displs_c"
      call exit(1)
    end if
    allocate(scounts_r(numDomain_r), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: scounts_r"
      call exit(1)
    end if
    allocate(scounts_c(numDomain_c), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: scounts_c"
      call exit(1)
    end if

    offset = 0
    do i = 1, numDomain_r
      displs_r(i) = offset
      if (i-1 < residueMol_r) then
        scounts_r(i) = numMolPerDomain_r + 1
      else
        scounts_r(i) = numMolPerDomain_r
      end if
      offset = offset + scounts_r(i)
    end do

    offset = 0
    do i = 1, numDomain_c
      displs_c(i) = offset
      if (i-1 < residueMol_c) then
        scounts_c(i) = numMolPerDomain_c + 1
      else
        scounts_c(i) = numMolPerDomain_c
      end if
      offset = offset + scounts_c(i)
    end do

    num_r = scounts_r(r_group_idx + 1)
    num_c = scounts_c(c_group_idx + 1)
    r_start = displs_r(r_group_idx + 1) + 1
    r_end = r_start + num_r - 1
    c_start = displs_c(c_group_idx + 1) + 1
    c_end = c_start + num_c - 1

    displs_r = displs_r * 3 * numFrame
    displs_c = displs_c * 3 * numFrame
    scounts_r = scounts_r * 3 * numFrame
    scounts_c = scounts_c * 3 * numFrame

    !check if myrank is at the ending boundary and if indexes are coincident
    if (r_group_idx == numDomain_r - 1) then
      if (r_end /= totNumMol) then
        write(*,*) "Error: r_end /= totNumMol, r_end =", r_end
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
    end if
    if (c_group_idx == numDomain_c - 1) then
      if (c_end /= totNumMol) then
        write(*,*) "Error: c_end /= totNumMol"
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
    end if
    if (myrank == root) then
      write(*,*) "numDomain_r x numDomain_c = ", numDomain_r, " x ", numDomain_c
    end if
  !  write(*,*) "my rank =", myrank
  !  write(*,*) "r_start, r_end =", r_start, r_end
  !  write(*,*) "c_start, c_end =", c_start, c_end
  !  write(*,*)
  end subroutine domainDecomposition

end module mpiproc
