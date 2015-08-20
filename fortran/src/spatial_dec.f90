module spatial_dec
  use mpiproc
  implicit none

  integer :: num_rBin
  integer, allocatable :: sd_binIndex(:)
  real(8) :: rBinWidth
  real(8), allocatable :: sdPairCount(:, :), sdCorr(:, :, :), pos(:, :, :), rBins(:)
  !sdCorr: spatially decomposed correlation (lag, rBin, molTypePairIndex)
  !sdPairCount: (num_rBin, molTypePairIndex)
  integer :: stat

  !MPI variables
  real(8), allocatable :: pos_r(:, :, :), pos_c(:, :, :)

contains
  subroutine sd_init()
    implicit none
    rBinWidth = 0.01
  end subroutine sd_init

  subroutine sd_prepPosMemory(numFrame, totNumMol, num_r, num_c)
    implicit none
    integer, intent(in) :: numFrame, totNumMol, num_r, num_c

    allocate(pos_r(3, numFrame, num_r), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_r"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 

    allocate(pos_c(3, numFrame, num_c), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_c"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 

    if (myrank == root) then
      allocate(pos(3, numFrame, totNumMol), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: pos"
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if 
    else
      !not root, allocate dummy pos to inhibit error messages
      allocate(pos(1, 1, 1), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: dummy pos on rank", myrank
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
    end if
  end subroutine sd_prepPosMemory

  subroutine com_pos(com_p, pos, start_index, sys, cell)
    use top, only : system
    implicit none
    real(8), dimension(:, :), intent(out) :: com_p
    real(8), dimension(:, :), intent(in) :: pos
    real(8), dimension(:, :), allocatable :: pos_gathered
    integer, dimension(:), intent(in) :: start_index
    type(system), intent(in) :: sys
    real(8), dimension(3), intent(in) :: cell
    integer :: d, i, j, k, idx_begin, idx_end, idx_com, num_atom

    idx_com = 0
    do i = 1, size(sys%mol)
      num_atom = size(sys%mol(i)%atom)
      allocate(pos_gathered(3, num_atom), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: pos_gathered"
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
      do j = 1, sys%mol(i)%num
        idx_begin = start_index(i) + (j-1) * num_atom
        idx_end = idx_begin + num_atom - 1
        idx_com = idx_com + 1
        call gatherMolPos(pos_gathered, pos(:, idx_begin:idx_end), cell)
        do d = 1, 3
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
    real(8), dimension(3), intent(in) :: cell
    real(8), dimension(3) :: ref_pos
    integer :: d

    ref_pos = pos(:, 1)
    do d = 1, 3
      pos_gathered(d, :) = pos(d, :) - ref_pos(d)
      pos_gathered(d, :) = ref_pos(d) + pos_gathered(d, :) - &
                           nint(pos_gathered(d, :) / cell(d)) * cell(d)
    end do
  end subroutine gatherMolPos

  subroutine sd_broadcastPos()
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
  end subroutine sd_broadcastPos

  subroutine sd_cal_num_rBin(cell)
    implicit none
    real(8), intent(in) :: cell(3)

    num_rBin = ceiling(cell(1) / 2d0 * sqrt(3d0) / rBinWidth)
    ! *sqrt(3) to accommodate the longest distance inside a cubic (diagonal)
    if (myrank == root) write(*,*) "num_rBin = ", num_rBin
  end subroutine sd_cal_num_rBin

  subroutine sd_prepCorrMemory(maxLag, numMolType, numFrame)
    implicit none
    integer, intent(in) :: maxLag, numMolType, numFrame
    integer :: numMolTypePair

    numMolTypePair = numMolType * (numMolType + 1) / 2
    allocate(sdCorr(maxLag+1, num_rBin, numMolTypePair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: sdCorr"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    sdCorr = 0d0

    allocate(sdPairCount(num_rBin, numMolTypePair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: sdPairCount"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    sdPairCount = 0d0

    allocate(sd_binIndex(numFrame), stat=stat)
    if (stat /= 0) then
      write(*,*) "Allocation failed: sd_binIndex"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end subroutine sd_prepCorrMemory

  subroutine sd_getBinIndex(r, c, cell, sd_binIndex)
    implicit none
    integer, intent(in) :: r, c
    real(8), intent(in) :: cell(3)
    integer, intent(out) :: sd_binIndex(:)
    real(8) :: pp(3, size(sd_binIndex))
    integer :: d

    pp = pos_r(:,:,r) - pos_c(:,:,c)
    do d = 1, 3
      pp(d, :) = pp(d, :) - nint(pp(d, :) / cell(d)) * cell(d)
    end do
    sd_binIndex = ceiling(sqrt(sum(pp*pp, 1)) / rBinWidth)
    where (sd_binIndex == 0)
      sd_binIndex = 1
    end where
!    where (sd_binIndex >= ceiling(cellLength / 2.d0 / rBinWidth))
!    where (sd_binIndex > num_rBin)
!      sd_binIndex = -1
!    end where
  end subroutine sd_getBinIndex

  subroutine sd_collectCorr()
    implicit none
    if (myrank == root) then
      call mpi_reduce(MPI_IN_PLACE, sdCorr, size(sdCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
      call mpi_reduce(MPI_IN_PLACE, sdPairCount, size(sdPairCount), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    else
      call mpi_reduce(sdCorr, dummy_null, size(sdCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
      call mpi_reduce(sdPairCount, dummy_null, size(sdPairCount), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    end if
  end subroutine sd_collectCorr

  subroutine sd_average(numFrame, numMolType, frameCount)
    use utility, only: getMolTypePairIndexFromTypes
    implicit none
    integer, intent(in) :: numFrame, numMolType, frameCount(:)
    integer :: i, t1, t2, n, molTypePairIndex

    sdPairCount = sdPairCount / numFrame
    do n = 1, numMolType * (numMolType + 1) / 2
      do i = 1, num_rBin
        sdCorr(:,i,n) = sdCorr(:,i,n) / frameCount / sdPairCount(i, n)
      end do
    end do

    do t2 = 1, numMolType
      do t1 = t2, numMolType
        if (t1 /= t2) then
          molTypePairIndex = getMolTypePairIndexFromTypes(t1, t2, numMolType)
          sdPairCount(:, molTypePairIndex) = sdPairCount(:, molTypePairIndex) / 2d0
        end if
      end do
    end do
  end subroutine sd_average

  subroutine sd_make_rBins()
    implicit none
    integer :: i, stat
    allocate(rBins(num_rBin), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: rBins"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    rBins = [ (i - 0.5d0, i = 1, num_rBin) ] * rBinWidth
  end subroutine sd_make_rBins

  subroutine sd_finish()
    implicit none
    deallocate(pos_r, pos_c, sdCorr, sdPairCount)
    ! sd_binIndex has been deallocated in the main program
  end subroutine sd_finish
end module spatial_dec
