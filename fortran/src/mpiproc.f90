module mpiproc
  use mpi
  use varpars, only: qnt_dim
  implicit none
  private
  public :: mpi_setup, domain_dec, mpi_abend
  integer, public :: ierr
  integer, public :: myrank, nprocs       ! rank number and total number of processes
  integer, public, parameter :: root = 0
  public :: mpi_comm_world, mpi_in_place, mpi_sum, mpi_double_precision
  public :: mpi_integer, mpi_character, mpi_min, mpi_max, mpi_info_null
  public :: mpi_wtime
  real(8), public :: starttime, endtime, starttime2
  integer, public :: r_start, r_end, c_start, c_end
  real(8), public :: dummy_null
  integer, public :: num_domain_r, num_domain_c, num_r, num_c
  integer, public :: row_comm, col_comm, r_group_idx, c_group_idx, offset
  integer, public, dimension(:), allocatable :: displs_r, displs_c, scounts_r, scounts_c

  integer :: nummol_per_domain_r, nummol_per_domain_c
  integer :: r_start_offset, c_start_offset
  integer :: residueMol_r, residueMol_c

contains
  subroutine mpi_setup(type)
    implicit none
    character(len=4) :: type

    if(type == 'init') then
      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      call mpi_comm_rank(mpi_comm_world, myrank, ierr)
    else if(type == 'stop') then
      call mpi_finalize(ierr)
    endif
  end subroutine mpi_setup

  subroutine mpi_abend()
    call mpi_abort(mpi_comm_world, 1, ierr);
    stop 1
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

  subroutine domain_dec(totnummol, numframe)
    !domain decomposition for atom pairs (num_domain_r * num_domain_c = nprocs)
    !nummol_per_domain_r * num_domain_r ~= totnummol
    implicit none
    integer, intent(in) :: totnummol, numframe
    integer :: i, stat

    if (num_domain_r == 0 .and. num_domain_c == 0) then
      num_domain_c = nint(sqrt(dble(nprocs)))
      do while(mod(nprocs, num_domain_c) /= 0)
        num_domain_c = num_domain_c - 1
      end do
      num_domain_r = nprocs / num_domain_c
    else if (num_domain_r > 0) then
      num_domain_c = nprocs / num_domain_r
    else if (num_domain_c > 0) then
      num_domain_c = nprocs / num_domain_r
    else
      write(*,*) "Invalid domain decomposition: ", num_domain_r, " x ", num_domain_c
      call mpi_abend()
    end if

    if (num_domain_r * num_domain_c /= nprocs) then
      write(*,*) "Domain decomposition failed: ", num_domain_r, " x ", num_domain_c, " /= ", nprocs
      call mpi_abend()
    end if

    !Determine row and column position for the node
    r_group_idx = mod(myrank, num_domain_r) !column-major mapping
    c_group_idx = myrank / num_domain_r

    !Split comm into row and column comms
    call mpi_comm_split(mpi_comm_world, c_group_idx, r_group_idx, col_comm, ierr)
    !color by row, rank by column
    call mpi_comm_split(mpi_comm_world, r_group_idx, c_group_idx, row_comm, ierr)
    !color by column, rank by row

    nummol_per_domain_r = totnummol / num_domain_r
    nummol_per_domain_c = totnummol / num_domain_c
    residueMol_r = mod(totnummol, num_domain_r)
    residueMol_c = mod(totnummol, num_domain_c)

    allocate(displs_r(num_domain_r), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: displs_r"
      call exit(1)
    end if
    allocate(displs_c(num_domain_c), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: displs_c"
      call exit(1)
    end if
    allocate(scounts_r(num_domain_r), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: scounts_r"
      call exit(1)
    end if
    allocate(scounts_c(num_domain_c), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: scounts_c"
      call exit(1)
    end if

    offset = 0
    do i = 1, num_domain_r
      displs_r(i) = offset
      ! distribute the residue molecules one by one from the first node
      if (i-1 < residueMol_r) then
        scounts_r(i) = nummol_per_domain_r + 1
      else
        scounts_r(i) = nummol_per_domain_r
      end if
      offset = offset + scounts_r(i)
    end do

    offset = 0
    do i = 1, num_domain_c
      displs_c(i) = offset
      ! distribute the residue molecules one by one from the first node
      if (i-1 < residueMol_c) then
        scounts_c(i) = nummol_per_domain_c + 1
      else
        scounts_c(i) = nummol_per_domain_c
      end if
      offset = offset + scounts_c(i)
    end do

    ! index in view of moleulce
    num_r = scounts_r(r_group_idx + 1)  ! number of molecules (in row/col dimension)
    num_c = scounts_c(c_group_idx + 1)  ! distributed to this rank

    r_start = displs_r(r_group_idx + 1) + 1  ! molecule index begins from 1, so plus 1
    c_start = displs_c(c_group_idx + 1) + 1

    r_end = r_start + num_r - 1
    c_end = c_start + num_c - 1

    ! in view of memory, for distributing array in parallel
    displs_r = displs_r * qnt_dim * numframe
    displs_c = displs_c * qnt_dim * numframe
    scounts_r = scounts_r * qnt_dim * numframe
    scounts_c = scounts_c * qnt_dim * numframe

    !check if myrank is at the ending boundary and if indexes are coincident
    if (r_group_idx == num_domain_r - 1) then
      if (r_end /= totnummol) then
        write(*,*) "Error: r_end /= totnummol, r_end =", r_end
        call mpi_abend()
      end if
    end if
    if (c_group_idx == num_domain_c - 1) then
      if (c_end /= totnummol) then
        write(*,*) "Error: c_end /= totnummol"
        call mpi_abend()
      end if
    end if
    if (myrank == root) then
      write(*,*) "num_domain_r x num_domain_c = ", num_domain_r, " x ", num_domain_c
    end if
  !  write(*,*) "my rank =", myrank
  !  write(*,*) "r_start, r_end =", r_start, r_end
  !  write(*,*) "c_start, c_end =", c_start, c_end
  !  write(*,*)
  end subroutine domain_dec

end module mpiproc
