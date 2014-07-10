program spatialDecompose_mpi
  use mpi
  use g96
  implicit none
  integer, parameter :: num_parArg = 8, stdout = 6
  integer, parameter :: num_argPerData = 2
  integer :: num_dataArg, i, j, k, n, totNumAtom, t
  character(len=128) :: outFilename 
  character(len=128) :: dataFilename
  type(handle) :: dataFileHandle
  integer :: numFrame, maxLag, num_rBin, stat, numAtomType, numFrameRead
  integer :: atomTypePairIndex, tmp_i, skip
  integer, allocatable :: numAtom(:), charge(:), rBinIndex(:), norm(:)
  character(len=10) :: tmp_str
  real(8) :: cell(3), timestep, rBinWidth, tmp_r
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :), vv(:)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: pos(:, :, :), vel(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), rho(:, :), sdCorr(:, :, :), rho_tmp(:, :), sdCorr_tmp(:, :, :)
  !sdCorr: spatially decomposed correlation (lag, rBin, atomTypePairIndex)
  !rho: (num_rBin, atomTypePairIndex)
  logical :: is_periodic

  !MPI variables
  integer :: ierr, nprocs, myrank
  integer:: numDomain_r, numDomain_c, numAtomPerDomain_r, numAtomPerDomain_c
  integer, parameter :: root = 0
  integer :: r_start, r_end, c_start, c_end
  real(8) :: starttime, endtime

  !initialize
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)

  num_dataArg = command_argument_count() - num_parArg

  !welcome message
!  write(*,*) "hello from rank ", myrank

  !root check the input arguments
  if (myrank == root) then
    is_periodic = .true.
    if (num_dataArg < num_argPerData .or. mod(num_dataArg, num_argPerData) /= 0) then
      write(*,*) "usage: $spatialdecompose <outfile> <infile.g96> <numFrameToRead> <skip> <maxlag> <rbinwidth(nm)> &
                  &<numDomain_r> <numDomain_c> <numatom1> <charge1> [<numatom2> <charge2>...]"
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

  call get_command_argument(3, tmp_str) ! in the unit of frame number
  read(tmp_str, *) numFrame 

  call get_command_argument(4, tmp_str) ! in the unit of frame number
  read(tmp_str, *) skip

  call get_command_argument(5, tmp_str) ! in the unit of frame number
  read(tmp_str, *) maxLag

  call get_command_argument(6, tmp_str)
  read(tmp_str, *) rBinWidth

  call get_command_argument(7, tmp_str) 
  read(tmp_str, *) numDomain_r

  call get_command_argument(8, tmp_str) 
  read(tmp_str, *) numDomain_c

  numAtomType = num_dataArg / num_argPerData

  !rank root output parameters read
  if (myrank == root) then
    write(*,*) "outFile = ", outFilename
    write(*,*) "inFile.g96 = ", dataFilename
    write(*,*) "numFrame= ", numFrame
    write(*,*) "maxLag = ", maxLag 
    write(*,*) "rBinWidth = ", rBinWidth
    write(*,*) "numAtomType = ", numAtomType
    write(*,*) "numDomain_r = ", numDomain_r
    write(*,*) "numDomain_c = ", numDomain_c
  end if

  allocate(numAtom(numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: numAtom"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  allocate(charge(numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: charge"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  do n = 1, numAtomType
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 1, tmp_str) 
    read(tmp_str, *) numAtom(n)
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 2, tmp_str) 
    read(tmp_str, *) charge(n)
  end do
  if (myrank == root) then
    write(*,*) "numAtom = ", numAtom
    write(*,*) "charge = ", charge
  end if
  totNumAtom = sum(numAtom)

  !domain decomposition for atom pairs (numDomain_r * numDomain_c = nprocs)
  !numAtomPerDomain_r * numDomain_r ~= totNumAtom
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

  numAtomPerDomain_r = totNumAtom / numDomain_r
  numAtomPerDomain_c = totNumAtom / numDomain_c
  r_start = mod(myrank, numDomain_r) * numAtomPerDomain_r + 1
  r_end = r_start + numAtomPerDomain_r - 1
  c_start = (myrank / numDomain_r) * numAtomPerDomain_c + 1
  c_end = c_start + numAtomPerDomain_c - 1
  !check if myrank is at the ending boundary
  if (mod(myrank, numDomain_r) == numDomain_r - 1) r_end = totNumAtom
  if (myrank >= (numDomain_c - 1) * numDomain_r) c_end = totNumAtom
  if (myrank == root) then
    write(*,*) "numDomain_r x numDomain_c = ", numDomain_r, " x ", numDomain_c 
  end if
!  write(*,*) "my rank =", myrank
!  write(*,*) "r_start, r_end =", r_start, r_end
!  write(*,*) "c_start, c_end =", c_start, c_end
!  write(*,*)

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

  !read trajectory at root
  if (myrank == root) then
    write(*,*) "start reading trajectory..."
    starttime = MPI_Wtime()
    allocate(pos_tmp(3, totNumAtom), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_tmp"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    allocate(vel_tmp(3, totNumAtom), stat=stat)
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
      call read_trajectory(dataFileHandle, totNumAtom, is_periodic, pos_tmp, vel_tmp, cell, time(i), stat)
      if (stat /= 0) then
        write(*,*) "Reading trajectory error"
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if 
      numFrameRead = numFrameRead + 1
      pos(:,i,:) = pos_tmp
      vel(:,i,:) = vel_tmp
      do j = 1, skip-1
        call read_trajectory(dataFileHandle, totNumAtom, is_periodic, pos_tmp, vel_tmp, cell, tmp_r, stat)
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
  call mpi_bcast(pos, 3*numFrame*totNumAtom, mpi_double_precision, root, MPI_COMM_WORLD, ierr)
  call mpi_bcast(vel, 3*numFrame*totNumAtom, mpi_double_precision, root, MPI_COMM_WORLD, ierr)
  call mpi_bcast(cell, 3, mpi_double_precision, root, MPI_COMM_WORLD, ierr)
  call mpi_bcast(numFrameRead, 1, mpi_double_precision, root, MPI_COMM_WORLD, ierr)
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished broadcasting trajectory. It took ", endtime - starttime, " seconds"

  ! *sqrt(3) to accommodate the longest distance inside a cubic (diagonal)
  num_rBin = ceiling(cell(1) / 2d0 * sqrt(3d0) / rBinWidth)
  if (myrank == root) write(*,*) "num_rBin = ", num_rBin

  !spatial decomposition correlation
  if (myrank == root) write(*,*) "start spatial decomposition"
  starttime = MPI_Wtime()
  allocate(sdCorr(maxLag+1, num_rBin, numAtomType*numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: sdCorr"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  sdCorr = 0d0

  allocate(rho(num_rBin, numAtomType*numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: rho"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  rho = 0d0

  allocate(rBinIndex(numFrameRead), stat=stat)
  if (stat /= 0) then
    write(*,*) "Allocation failed: rBinIndex"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  allocate(vv(numFrameRead))
  if (stat /=0) then
    write(*,*) "Allocation failed: vv"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  do j = c_start, c_end
    do i = r_start, r_end
      if (i /= j) then
        call getBinIndex(pos(:,:,i), pos(:,:,j), cell(1), rBinWidth, rBinIndex)
        atomTypePairIndex = getAtomTypePairIndex(i, j, numAtom)
        do k = 1, maxLag+1      
          vv = sum(vel(:, k:numFrameRead, i) * vel(:, 1:numFrameRead-k+1, j), 1)
          do n = 1, numFrameRead-k+1
            tmp_i = rBinIndex(n)
            if (tmp_i <= num_rBin) then
              sdCorr(k, tmp_i, atomTypePairIndex) = sdCorr(k, tmp_i, atomTypePairIndex) + vv(n)
            end if
          end do
        end do
        do t = 1, numFrameRead
          tmp_i = rBinIndex(t)
          if (tmp_i <= num_rBin) then
            rho(tmp_i, atomTypePairIndex) = rho(tmp_i, atomTypePairIndex) + 1d0
          end if
        end do
      end if
    end do
  end do
  deallocate(rBinIndex)
  deallocate(pos)
  deallocate(vel)
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished spatial decomposition. It took ", endtime - starttime, " seconds"

  !collect sdCorr and rho
  if (myrank == root) write(*,*) "start collecting results"
  starttime = MPI_Wtime()
  if (myrank == root) then
    call mpi_reduce(MPI_IN_PLACE, sdCorr, size(sdCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(MPI_IN_PLACE, rho, size(rho), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  else
    call mpi_reduce(sdCorr, NULL(), size(sdCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    call mpi_reduce(rho, NULL(), size(rho), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
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

    rho = rho / numFrameRead
    sdCorr = sdCorr / 3d0

    norm = [ (numFrameRead - (i-1), i = 1, maxLag+1) ]
    forall (i = 1:num_rBin, n = 1:numAtomType*numAtomType )
      sdCorr(:,i,n) = sdCorr(:,i,n) / norm
    end forall

    deallocate(norm)

!    forall (i = 1:maxLag+1, n = 1:numAtomType*numAtomType )
!      sdCorr(i,:,n) = sdCorr(i,:,n) / rho(:, n)
!    end forall
!    where (isnan(sdCorr))
!      sdCorr = 0
!    end where

    do n = 1, numAtomType*numAtomType
      do j = 1, num_rBin
        if (rho(j, n) > 0d0) then
            sdCorr(:,j,n) = sdCorr(:,j,n) / rho(j, n)
        else if (any(abs(sdCorr(:,j,n)) > 0d0)) then
          !rho is zero, sdCorr should also be zero
          write(*,*) "error: rho(",j,",",n,") = ", rho(j,n)
          write(*,*) "but sdCorr(:,", j, ",", n, ") are not all zero"
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if
      end do
    end do

    !output results
    call output()
  end if

  call mpi_finalize(ierr)
  stop

contains
  elemental real(8) function wrap(r, l)
    implicit none
    real(8), intent(in) :: r, l
    real(8) :: half_l
    half_l = l / 2.d0
    wrap = r
    do while (wrap > half_l)
      wrap = abs(wrap - l)
    end do
  end function wrap
  
  subroutine getBinIndex(p1, p2, cellLength, rBinWidth, rBinIndex)
    implicit none
    real(8), intent(in) :: p1(:,:), p2(:,:), cellLength, rBinWidth
    !p1(dim,timeFrame)
    integer, intent(out) :: rBinIndex(:)
    real(8) :: pp(size(p1,1), size(p1,2))
    
    pp = abs(p1 - p2)
    pp = wrap(pp, cellLength)
    rBinIndex = ceiling(sqrt(sum(pp*pp, 1)) / rBinWidth)
    where (rBinIndex == 0)
      rBinIndex = 1
    end where
!    where (rBinIndex >= ceiling(cellLength / 2.d0 / rBinWidth))
!    where (rBinIndex > num_rBin)
!      rBinIndex = -1
!    end where
  end subroutine getBinIndex

  integer function getAtomTypeIndex(i, numAtom)
    implicit none
    integer, intent(in) :: i, numAtom(:)
    integer :: n, numAtom_acc
!    integer, save :: numAtomType = size(numAtom)

    getAtomTypeIndex = -1
    numAtom_acc = 0
    do n = 1, numAtomType
      numAtom_acc = numAtom_acc + numAtom(n)
      if (i <= numAtom_acc) then
        getAtomTypeIndex = n
        return
      end if
    end do
  end function getAtomTypeIndex
  
  integer function getAtomTypePairIndex(i, j, numAtom)
    implicit none
    integer, intent(in) :: i, j, numAtom(:)
    integer :: ii, jj
!    integer, save :: numAtomType = size(numAtom)
  
    ii = getAtomTypeIndex(i, numAtom)
    jj = getAtomTypeIndex(j, numAtom)
    getAtomTypePairIndex = (ii-1)*numAtomType + jj
  end function getAtomTypePairIndex

  subroutine output()
    use octave_save
    implicit none
    type(handle) :: htraj
    real(8), allocatable :: timeLags(:), rBins(:)

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

    call create_octave(htraj, outFilename)
    call write_octave_scalar(htraj, "timestep", timestep)
    call write_octave_vec(htraj, "charge", dble(charge))
    call write_octave_vec(htraj, "numAtom", dble(numAtom))
    call write_octave_vec(htraj, "cell", cell)
    call write_octave_vec(htraj, "timeLags", timeLags)
    call write_octave_vec(htraj, "rBins", rBins)
    call write_octave_mat3(htraj, "sdCorr", sdCorr)
    call write_octave_mat2(htraj, "rho", rho)
    call close_octave(htraj)
  end subroutine output

end program spatialDecompose_mpi
