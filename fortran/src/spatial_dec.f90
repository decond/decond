module spatial_dec
  use mpiproc
  use varpars, only: numframe, totnummol, world_dim
  implicit none
  private
  public sd_init, com_pos, sd_prep_corrmemory, sd_getbinindex, &
         sd_cal_num_rbin, sd_broadcastpos, sd_prep, &
         sd_collectcorr, sd_average, sd_make_rbins, sd_finish
  real(8), public :: rbinwidth
  !sdcorr: spatially decomposed correlation (lag, rBin, moltypepair_idx)
  !sdpaircount: (num_rbin, moltypepair_idx)
  real(8), public, allocatable :: sdpaircount(:, :), sdcorr(:, :, :), pos(:, :, :)
  integer, public :: num_rbin
  integer, public, allocatable :: sd_binIndex(:)
  real(8), public, allocatable :: rbins(:)
  !MPI variables
  real(8), public, allocatable :: pos_r(:, :, :), pos_c(:, :, :)

  integer :: stat

contains
  subroutine sd_init()
    rbinwidth = 0.01
  end subroutine sd_init

  subroutine sd_prep()
    allocate(pos_r(world_dim, numframe, num_r), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_r"
      call mpi_abend()
    end if 

    allocate(pos_c(world_dim, numframe, num_c), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_c"
      call mpi_abend()
    end if 

    if (myrank == root) then
      allocate(pos(world_dim, numframe, totnummol), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: pos"
        call mpi_abend()
      end if 
    else
      !not root, allocate dummy pos to inhibit error messages
      allocate(pos(1, 1, 1), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: dummy pos on rank", myrank
        call mpi_abend()
      end if
    end if
  end subroutine sd_prep

  subroutine com_pos(com_p, pos, start_index, sys, cell)
    use top, only : system
    implicit none
    real(8), dimension(:, :), intent(out) :: com_p
    real(8), dimension(:, :), intent(in) :: pos
    real(8), dimension(:, :), allocatable :: pos_gathered
    integer, dimension(:), intent(in) :: start_index
    type(system), intent(in) :: sys
    real(8), dimension(world_dim), intent(in) :: cell
    integer :: d, i, j, k, idx_begin, idx_end, idx_com, num_atom

    idx_com = 0
    do i = 1, size(sys%mol)
      num_atom = size(sys%mol(i)%atom)
      allocate(pos_gathered(world_dim, num_atom), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: pos_gathered"
        call mpi_abend()
      end if
      do j = 1, sys%mol(i)%num
        idx_begin = start_index(i) + (j-1) * num_atom
        idx_end = idx_begin + num_atom - 1
        idx_com = idx_com + 1
        call gatherMolPos(pos_gathered, pos(:, idx_begin:idx_end), cell)
        do d = 1, world_dim
          com_p(d, idx_com) = sum(pos_gathered(d, :) * sys%mol(i)%atom(:)%mass) / sum(sys%mol(i)%atom(:)%mass)
        end do
      end do
      deallocate(pos_gathered)
    end do
  end subroutine com_pos

  ! make sure the gathered positions are from the same moleucle, instead of different images
  subroutine gatherMolPos(pos_gathered, pos, cell)
    implicit none
    real(8), dimension(:, :), intent(out) :: pos_gathered
    real(8), dimension(:, :), intent(in) :: pos
    real(8), dimension(world_dim), intent(in) :: cell
    real(8), dimension(world_dim) :: ref_pos
    integer :: d

    ref_pos = pos(:, 1)
    do d = 1, world_dim
      pos_gathered(d, :) = pos(d, :) - ref_pos(d)
      pos_gathered(d, :) = ref_pos(d) + pos_gathered(d, :) - &
                           nint(pos_gathered(d, :) / cell(d)) * cell(d)
    end do
  end subroutine gatherMolPos

  subroutine sd_broadcastpos()
    if (r_group_idx == 0) then
      call mpi_scatterv(pos, scounts_c, displs_c, mpi_double_precision, pos_c,&
                        scounts_c(c_group_idx + 1), mpi_double_precision, root, row_comm, ierr)
    end if
    call mpi_bcast(pos_c, scounts_c(c_group_idx + 1), mpi_double_precision, root, col_comm, ierr)

    if (c_group_idx == 0) then
      call mpi_scatterv(pos, scounts_r, displs_r, mpi_double_precision, pos_r,&
                        scounts_r(r_group_idx + 1), mpi_double_precision, root, col_comm, ierr)
    end if
    call mpi_bcast(pos_r, scounts_r(r_group_idx + 1), mpi_double_precision, root, row_comm, ierr)
    deallocate(pos)
  end subroutine sd_broadcastpos

  subroutine sd_cal_num_rbin(cell)
    implicit none
    real(8), intent(in) :: cell(world_dim)

    num_rbin = ceiling(cell(1) / 2d0 * sqrt(3d0) / rbinwidth)
    ! *sqrt(3) to accommodate the longest distance inside a cubic (diagonal)
    if (myrank == root) write(*,*) "num_rbin = ", num_rbin
  end subroutine sd_cal_num_rbin

  subroutine sd_prep_corrmemory(maxlag, nummoltype, numframe)
    implicit none
    integer, intent(in) :: maxlag, nummoltype, numframe
    integer :: num_moltypepair

    num_moltypepair = nummoltype * (nummoltype + 1) / 2
    allocate(sdcorr(maxlag+1, num_rbin, num_moltypepair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: sdcorr"
      call mpi_abend()
    end if
    sdcorr = 0d0

    allocate(sdpaircount(num_rbin, num_moltypepair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: sdpaircount"
      call mpi_abend()
    end if
    sdpaircount = 0d0

    allocate(sd_binIndex(numframe), stat=stat)
    if (stat /= 0) then
      write(*,*) "Allocation failed: sd_binIndex"
      call mpi_abend()
    end if
  end subroutine sd_prep_corrmemory

  subroutine sd_getbinindex(r, c, cell, sd_binIndex)
    implicit none
    integer, intent(in) :: r, c
    real(8), intent(in) :: cell(world_dim)
    integer, intent(out) :: sd_binIndex(:)
    real(8) :: pp(world_dim, size(sd_binIndex))
    integer :: d

    pp = pos_r(:,:,r) - pos_c(:,:,c)
    do d = 1, world_dim
      pp(d, :) = pp(d, :) - nint(pp(d, :) / cell(d)) * cell(d)
    end do
    sd_binIndex = ceiling(sqrt(sum(pp*pp, 1)) / rbinwidth)
    where (sd_binIndex == 0)
      sd_binIndex = 1
    end where
!    where (sd_binIndex >= ceiling(cellLength / 2.d0 / rbinwidth))
!    where (sd_binIndex > num_rbin)
!      sd_binIndex = -1
!    end where
  end subroutine sd_getbinindex

  subroutine sd_collectcorr()
    implicit none
    if (myrank == root) then
      write(*,*) "collecting sdcorr"
      call mpi_reduce(MPI_IN_PLACE, sdcorr, size(sdcorr), mpi_double_precision, MPI_SUM, root, mpi_comm_world, ierr)
    else
      call mpi_reduce(sdcorr, dummy_null, size(sdcorr), mpi_double_precision, MPI_SUM, root, mpi_comm_world, ierr)
    end if
    call mpi_barrier(mpi_comm_world, ierr)

    if (myrank == root) then
      write(*,*) "collecting sdpaircount"
      call mpi_reduce(MPI_IN_PLACE, sdpaircount, size(sdpaircount), mpi_double_precision, MPI_SUM, root, mpi_comm_world, ierr)
    else
      call mpi_reduce(sdpaircount, dummy_null, size(sdpaircount), mpi_double_precision, MPI_SUM, root, mpi_comm_world, ierr)
    end if
    call mpi_barrier(mpi_comm_world, ierr)
  end subroutine sd_collectcorr

  subroutine sd_average(numframe, nummoltype, framecount)
    use utility, only: get_pairindex_upper_diag
    implicit none
    integer, intent(in) :: numframe, nummoltype, framecount(:)
    integer :: i, t1, t2, n, moltypepair_idx

    sdpaircount = sdpaircount / numframe
    do n = 1, nummoltype * (nummoltype + 1) / 2
      do i = 1, num_rbin
        sdcorr(:,i,n) = sdcorr(:,i,n) / framecount / sdpaircount(i, n)
      end do
    end do

    do t2 = 1, nummoltype
      do t1 = t2, nummoltype
        if (t1 /= t2) then
          moltypepair_idx = get_pairindex_upper_diag(t1, t2, nummoltype)
          sdpaircount(:, moltypepair_idx) = sdpaircount(:, moltypepair_idx) / 2d0
        end if
      end do
    end do
  end subroutine sd_average

  subroutine sd_make_rbins()
    implicit none
    integer :: i, stat
    allocate(rbins(num_rbin), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: rbins"
      call mpi_abend()
    end if 
    rbins = [ (i - 0.5d0, i = 1, num_rbin) ] * rbinwidth
  end subroutine sd_make_rbins

  subroutine sd_finish()
    implicit none
    deallocate(pos_r, pos_c, sdcorr, sdpaircount)
    ! sd_binIndex has been deallocated in the main program
  end subroutine sd_finish
end module spatial_dec
