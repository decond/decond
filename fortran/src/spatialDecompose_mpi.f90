program spatialDecompose_mpi
  use mpi
  use utility, only : handle
  use xdr, only : open_trajectory, close_trajectory, read_trajectory, get_natom
  use top, only : open_top, close_top, read_top, system, print_sys
  implicit none
  integer, parameter :: num_parArg = 9
  integer, parameter :: num_argPerData = 2
  integer :: num_dataArg, i, j, k, n, totNumMol, t, sysNumAtom
  character(len=128) :: outFilename 
  character(len=128) :: dataFilename
  character(len=128) :: topFilename
  type(handle) :: dataFileHandle, topFileHandle
  integer :: numFrame, maxLag, num_rBin, stat, numMolType, numFrameRead
  integer :: molTypePairIndex, tmp_i, skip
  integer, allocatable :: charge(:), rBinIndex(:), norm(:), start_index(:)
  character(len=10) :: tmp_str
  real(8) :: cell(3), timestep, rBinWidth, tmp_r
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :), vv(:)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: pos(:, :, :), vel(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), rho(:, :), sdCorr(:, :, :), dummy_null
  !sdCorr: spatially decomposed correlation (lag, rBin, molTypePairIndex)
  !rho: (num_rBin, molTypePairIndex)
  logical :: is_periodic
  type(system) :: sys

  !MPI variables
  integer :: ierr, nprocs, myrank
  integer:: numDomain_r, numDomain_c, numMolPerDomain_r, numMolPerDomain_c
  integer, parameter :: root = 0
  integer :: r_start, r_end, c_start, c_end
  real(8) :: starttime, endtime
  integer :: r_start_offset, c_start_offset
  integer :: residueMol_r, residueMol_c, num_r, num_c
  real(8), allocatable :: pos_r(:, :, :), pos_c(:, :, :), vel_r(:, :, :), vel_c(:, :, :)
  integer :: row_comm, col_comm, r_group_idx, c_group_idx, offset
  integer, dimension(:), allocatable :: displs_r, displs_c, scounts_r, scounts_c

  !initialize
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)

  num_dataArg = command_argument_count() - num_parArg

  !root checks the input arguments
  if (myrank == root) then
    is_periodic = .true.
    if (num_dataArg < num_argPerData .or. mod(num_dataArg, num_argPerData) /= 0) then
      write(*,*) "usage: $spatialdecompose <outfile> <infile.trr> <topfile.top> <numFrameToRead> &
                  <skip> <maxlag> <rbinwidth(nm)> <numDomain_r> <numDomain_c> &
                  <molecule1> <start_index1> [<molecule2> <start_index2>...]"
      write(*,*) "Note: skip=1 means no frames are skipped. skip=2 means reading every 2nd frame."
      write(*,*) "Note: maxlag is counted in terms of the numFrameToRead."
      write(*,*) "Note: numDomain_r and numDomain_c can be 0 and let the program decide (may not be the best)."
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end if

  !read parameters for all ranks
  call get_command_argument(1, outFilename)

  call get_command_argument(2, dataFilename)

  call get_command_argument(3, topFilename)

  call get_command_argument(4, tmp_str) ! in the unit of frame number
  read(tmp_str, *) numFrame 

  call get_command_argument(5, tmp_str) ! in the unit of frame number
  read(tmp_str, *) skip

  call get_command_argument(6, tmp_str) ! in the unit of frame number
  read(tmp_str, *) maxLag

  call get_command_argument(7, tmp_str)
  read(tmp_str, *) rBinWidth

  call get_command_argument(8, tmp_str) 
  read(tmp_str, *) numDomain_r

  call get_command_argument(9, tmp_str) 
  read(tmp_str, *) numDomain_c

  numMolType = num_dataArg / num_argPerData

  !rank root output parameters read
  if (myrank == root) then
    write(*,*) "outFile = ", outFilename
    write(*,*) "inFile.trr = ", dataFilename
    write(*,*) "topFile.trr = ", topFilename
    write(*,*) "numFrame= ", numFrame
    write(*,*) "maxLag = ", maxLag 
    write(*,*) "rBinWidth = ", rBinWidth
    write(*,*) "numMolType = ", numMolType
    write(*,*) "numDomain_r = ", numDomain_r
    write(*,*) "numDomain_c = ", numDomain_c
  end if

  allocate(sys%mol(numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: sys%mol"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  allocate(charge(numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: charge"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  allocate(start_index(numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: start_index"
    call exit(1)
  end if 

  do n = 1, numMolType
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 1, sys%mol(n)%type) 
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 2, tmp_str) 
    read(tmp_str, *) start_index(n)
  end do
  if (myrank == root) then
    write(*,*) "sys%mol%type = ", sys%mol%type
    write(*,*) "start_index = ", start_index
  end if

  !read topFile
  topFileHandle = open_top(topFilename)
  call read_top(topFileHandle, sys)
  call close_top(topFileHandle)
  if (myrank == root) call print_sys(sys)

  do n = 1, numMolType
    charge(n) = sum(sys%mol(n)%atom(:)%charge)
  end do
  totNumMol = sum(sys%mol(:)%num)

  !domain decomposition for atom pairs (numDomain_r * numDomain_c = nprocs)
  !numMolPerDomain_r * numDomain_r ~= totNumMol
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

  !prepare memory for all ranks
  allocate(pos_r(3, numFrame, num_r), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: pos_r"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  allocate(vel_r(3, numFrame, num_r), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel_r"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  allocate(pos_c(3, numFrame, num_c), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: pos_c"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  allocate(vel_c(3, numFrame, num_c), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel_r"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  !read trajectory at root
  if (myrank == root) then
    write(*,*) "start reading trajectory..."
    starttime = MPI_Wtime()
    sysNumAtom = get_natom(dataFilename)
    write(*,*) "sysNumAtom=", sysNumAtom

    allocate(pos(3, numFrame, totNumMol), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    allocate(vel(3, numFrame, totNumMol), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: vel"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    allocate(pos_tmp(3, sysNumAtom), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_tmp"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    allocate(vel_tmp(3, sysNumAtom), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: vel_tmp"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    allocate(time(numFrame), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: time"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 

    numFrameRead = 0
    call open_trajectory(dataFileHandle, dataFilename)
    do i = 1, numFrame
      call read_trajectory(dataFileHandle, sysNumAtom, is_periodic, pos_tmp, vel_tmp, cell, time(i), stat)
      if (stat /= 0) then
        write(*,*) "Reading trajectory error"
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if 
      numFrameRead = numFrameRead + 1
      call com_pos(pos(:, i, :), pos_tmp, start_index, sys, cell)
      call com_vel(vel(:, i, :), vel_tmp, start_index, sys)
      do j = 1, skip-1
        call read_trajectory(dataFileHandle, sysNumAtom, is_periodic, pos_tmp, vel_tmp, cell, tmp_r, stat)
        if (stat > 0) then
          write(*,*) "Reading trajectory error"
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        else if (stat < 0) then
          !end of file
          exit
        end if 
      end do
    end do
    call close_trajectory(dataFileHandle)
    if (myrank == root) write(*,*) "numFrameRead = ", numFrameRead
    if (numFrameRead /= numFrame) then
      write(*,*) "Number of frames expected to read is not the same as actually read!"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    timestep = time(2) - time(1)
    deallocate(pos_tmp)
    deallocate(vel_tmp)
    deallocate(time)
    endtime = MPI_Wtime()
    write(*,*) "finished reading trajectory. It took ", endtime - starttime, "seconds"
    write(*,*) "timestep = ", timestep
    write(*,*) "cell = ", cell
  end if

  !distribute trajectory data collectively
  if (myrank == root) write(*,*) "start broadcasting trajectory"
  starttime = MPI_Wtime()
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

  if (r_group_idx == 0) then
    call mpi_scatterv(vel, scounts_c, displs_c, mpi_double_precision, vel_c,&
                      scounts_c(c_group_idx + 1), mpi_double_precision, root, row_comm, ierr)
  end if
  call mpi_bcast(vel_c, scounts_c(c_group_idx + 1), mpi_double_precision, root, col_comm, ierr)

  if (c_group_idx == 0) then
    call mpi_scatterv(vel, scounts_r, displs_r, mpi_double_precision, vel_r,&
                      scounts_r(r_group_idx + 1), mpi_double_precision, root, col_comm, ierr)
  end if
  call mpi_bcast(vel_r, scounts_r(r_group_idx + 1), mpi_double_precision, root, row_comm, ierr)

  call mpi_bcast(cell, 3, mpi_double_precision, root, MPI_COMM_WORLD, ierr)
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished broadcasting trajectory. It took ", endtime - starttime, " seconds"

  if (myrank == root) then
    deallocate(pos)
    deallocate(vel)
  end if

  ! *sqrt(3) to accommodate the longest distance inside a cubic (diagonal)
  num_rBin = ceiling(cell(1) / 2d0 * sqrt(3d0) / rBinWidth)
  if (myrank == root) write(*,*) "num_rBin = ", num_rBin

  !spatial decomposition correlation
  if (myrank == root) write(*,*) "start spatial decomposition"
  starttime = MPI_Wtime()
  allocate(sdCorr(maxLag+1, num_rBin, numMolType*numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: sdCorr"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  sdCorr = 0d0

  allocate(rho(num_rBin, numMolType*numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: rho"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  rho = 0d0

  allocate(rBinIndex(numFrame), stat=stat)
  if (stat /= 0) then
    write(*,*) "Allocation failed: rBinIndex"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  allocate(vv(numFrame))
  if (stat /=0) then
    write(*,*) "Allocation failed: vv"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  do j = c_start, c_end
    do i = r_start, r_end
      if (i /= j) then
        call getBinIndex(pos_r(:,:,i-r_start+1), pos_c(:,:,j-c_start+1), cell(1), rBinWidth, rBinIndex)
        molTypePairIndex = getMolTypePairIndex(i, j, sys%mol(:)%num)
        do k = 1, maxLag+1      
          vv = sum(vel_r(:, k:numFrame, i-r_start+1) * vel_c(:, 1:numFrame-k+1, j-c_start+1), 1)
          do n = 1, numFrame-k+1
            tmp_i = rBinIndex(n)
            if (tmp_i <= num_rBin) then
              sdCorr(k, tmp_i, molTypePairIndex) = sdCorr(k, tmp_i, molTypePairIndex) + vv(n)
            end if
          end do
        end do
        do t = 1, numFrame
          tmp_i = rBinIndex(t)
          if (tmp_i <= num_rBin) then
            rho(tmp_i, molTypePairIndex) = rho(tmp_i, molTypePairIndex) + 1d0
          end if
        end do
      end if
    end do
  end do
  deallocate(rBinIndex)
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished spatial decomposition. It took ", endtime - starttime, " seconds"

  !collect sdCorr and rho
  if (myrank == root) write(*,*) "start collecting results"
  starttime = MPI_Wtime()
  if (myrank == root) then
    call mpi_reduce(MPI_IN_PLACE, sdCorr, size(sdCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(MPI_IN_PLACE, rho, size(rho), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  else
    call mpi_reduce(sdCorr, dummy_null, size(sdCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(rho, dummy_null, size(rho), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  end if
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished collecting results. It took ", endtime - starttime, " seconds"

  !normalization at root and then output
  if (myrank == root) then
    allocate(norm(maxLag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: norm"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 

    rho = rho / numFrame
    sdCorr = sdCorr / 3d0

    norm = [ (numFrame - (i-1), i = 1, maxLag+1) ]
    forall (i = 1:num_rBin, n = 1:numMolType*numMolType )
      sdCorr(:,i,n) = sdCorr(:,i,n) / norm
    end forall

    deallocate(norm)

!    forall (i = 1:maxLag+1, n = 1:numMolType*numMolType )
!      sdCorr(i,:,n) = sdCorr(i,:,n) / rho(:, n)
!    end forall
!    where (isnan(sdCorr))
!      sdCorr = 0
!    end where

    do n = 1, numMolType*numMolType
      do j = 1, num_rBin
        sdCorr(:,j,n) = sdCorr(:,j,n) / rho(j, n)
      end do
    end do

    !output results
    write(*,*) "start writing outputs"
    starttime = MPI_Wtime()
    call output()
    endtime = MPI_Wtime()
    write(*,*) "finished writing outputs. It took ", endtime - starttime, " seconds"
  end if

  call mpi_finalize(ierr)
  stop

contains
  subroutine getBinIndex(p1, p2, cell, rBinWidth, rBinIndex)
    implicit none
    real(8), intent(in) :: p1(:,:), p2(:,:), cell(3), rBinWidth
    !p1(dim,timeFrame)
    integer, intent(out) :: rBinIndex(:)
    real(8) :: pp(size(p1,1), size(p1,2))
    integer :: d
    
    pp = p1 - p2
    do d = 1, 3
      pp(d, :) = pp(d, :) - nint(pp(d, :) / cell(d)) * cell(d)
    end do
    rBinIndex = ceiling(sqrt(sum(pp*pp, 1)) / rBinWidth)
    where (rBinIndex == 0)
      rBinIndex = 1
    end where
!    where (rBinIndex >= ceiling(cellLength / 2.d0 / rBinWidth))
!    where (rBinIndex > num_rBin)
!      rBinIndex = -1
!    end where
  end subroutine getBinIndex

  integer function getMolTypeIndex(i, numMol)
    implicit none
    integer, intent(in) :: i, numMol(:)
    integer :: n, numMol_acc
!    integer, save :: numMolType = size(numMol)

    getMolTypeIndex = -1
    numMol_acc = 0
    do n = 1, numMolType
      numMol_acc = numMol_acc + numMol(n)
      if (i <= numMol_acc) then
        getMolTypeIndex = n
        return
      end if
    end do
  end function getMolTypeIndex
  
  integer function getMolTypePairIndex(i, j, numMol)
    implicit none
    integer, intent(in) :: i, j, numMol(:)
    integer :: ii, jj
!    integer, save :: numMolType = size(numMol)
  
    ii = getMolTypeIndex(i, numMol)
    jj = getMolTypeIndex(j, numMol)
    getMolTypePairIndex = (ii-1)*numMolType + jj
  end function getMolTypePairIndex

  subroutine com_pos(com_p, pos, start_index, sys, cell)
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

  subroutine com_vel(com_v, vel, start_index, sys)
    implicit none
    real(8), dimension(:, :), intent(out) :: com_v
    real(8), dimension(:, :), intent(in) :: vel 
    integer, dimension(:), intent(in) :: start_index
    type(system), intent(in) :: sys
    integer :: d, i, j, k, idx_begin, idx_end, idx_com

    com_v = 0d0
    idx_com = 0
    do i = 1, size(sys%mol)
      do j = 1, sys%mol(i)%num
        idx_begin = start_index(i) + (j-1) * size(sys%mol(i)%atom)
        idx_end = idx_begin + size(sys%mol(i)%atom) - 1
        idx_com = idx_com + 1
        do d = 1, 3
          com_v(d, idx_com) = com_v(d, idx_com) + &
              sum(vel(d, idx_begin:idx_end) * sys%mol(i)%atom(:)%mass) / sum(sys%mol(i)%atom(:)%mass)
        end do
      end do
    end do
  end subroutine com_vel

  subroutine output()
    use H5DS
    use H5LT
    use HDF5
    implicit none
    real(8), allocatable :: timeLags(:), rBins(:)
    integer :: ierr
    integer(hid_t) :: fid, did1, did2, sid1, sid2

    allocate(timeLags(maxLag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: timeLags"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    timeLags = [ (dble(i), i = 0, maxLag) ] * timestep
    
    allocate(rBins(num_rBin), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: rBins"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    rBins = [ (i - 0.5d0, i = 1, num_rBin) ] * rBinWidth

    call H5open_f(ierr)

    ! create a HDF5 file
    call H5Fcreate_f(outFilename, H5F_ACC_TRUNC_F, fid, ierr)

    ! create and write dataset
    call H5LTset_attribute_double_f(fid, "/", "timestep", [timestep], int(1, kind=size_t), ierr)
    call H5LTset_attribute_int_f(fid, "/", "charge", charge, size(charge, kind=size_t), ierr)
    call H5LTset_attribute_int_f(fid, "/", "numMol", sys%mol(:)%num, size(sys%mol(:)%num, kind=size_t), ierr)
    call H5LTset_attribute_double_f(fid, "/", "cell", cell, size(cell, kind=size_t), ierr)

    call H5LTmake_dataset_double_f(fid, "sdCorr", 3, &
        [size(sdCorr, 1, kind=hsize_t), size(sdCorr, 2, kind=hsize_t), size(sdCorr, 3, kind=hsize_t)], sdCorr, ierr)
    call H5Dopen_f(fid, "sdCorr", did1, ierr)

    call H5LTmake_dataset_double_f(fid, "rho", 2, &
        [size(rho, 1, kind=hsize_t), size(rho, 2, kind=hsize_t)], rho, ierr)
    call H5Dopen_f(fid, "rho", did2, ierr)

    call H5LTmake_dataset_double_f(fid, "timeLags", 1, [size(timeLags, kind=hsize_t)], timeLags, ierr)
    call H5Dopen_f(fid, "timeLags", sid1, ierr)

    call H5LTmake_dataset_double_f(fid, "rBins", 1, [size(rBins, kind=hsize_t)], rBins, ierr)
    call H5Dopen_f(fid, "rBins", sid2, ierr)

    ! attach scale dimension
    call H5DSattach_scale_f(did1, sid1, 1, ierr)
    call H5DSattach_scale_f(did1, sid2, 2, ierr)
    call H5DSattach_scale_f(did2, sid2, 1, ierr)

    call H5Dclose_f(sid1, ierr)
    call H5Dclose_f(sid2, ierr)
    call H5Dclose_f(did1, ierr)
    call H5Dclose_f(did2, ierr)
    call H5Fclose_f(fid, ierr)
    call H5close_f(ierr)
  end subroutine output

end program spatialDecompose_mpi
