program spatialDecompose_mpi
  use mpi
  use g96
  implicit none
  integer, parameter :: num_parArg = 7
  integer, parameter :: num_argPerData = 2
  integer :: num_dataArg, i, j, k, n, totNumAtom
  character(len=128) :: outFilename 
  character(len=128) :: dataFilename
  type(handle) :: dataFileHandle
  integer :: numFrame, maxLag, num_rBin, stat, numAtomType
  integer :: atomTypePairIndex, tmp_i
  integer, allocatable :: numAtom(:), charge(:), rBinIndex(:), norm(:)
  character(len=10) :: numFrame_str, maxLag_str, rBinWidth_str, charge_str, numAtom_str
  real(8) :: cell(3), timestep, rBinWidth
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: pos(:, :, :), vel(:, :, :)
!  real(8), allocatable :: pos_r(:, :, :), pos_c(:, :, :), vel_r(:, :, :), vel_c(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), rho(:, :), sdCorr(:, :, :), timeLags(:), rBins(:)
  !sdCorr: spatially decomposed correlation (lag, rBin, atomTypePairIndex)
  !rho: (num_rBin, atomTypePairIndex)
  logical :: is_periodic

  !MPI variables
  integer :: ierr, nprocs, myrank
  character(len=10) :: numDomain_r_str, numDomain_c_str
  integer:: numDomain_r, numDomain_c, numAtomPerDomain_r, numAtomPerDomain_c
  integer, parameter :: root = 0

  !initialize
  call mpi_init(ierr) !MPI
  call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)

  !root check the input arguments
  if (myrank == root)
    is_periodic = .true.

    num_dataArg = command_argument_count() - num_parArg
    if (num_dataArg < num_argPerData .or. mod(num_dataArg, num_argPerData) /= 0) then
      write(*,*) "Usage: $spatialDecompose <outFile> <inFile.g96> <numFrame> <maxLag> <rBinWidth(nm)> &
                  &<numDomain_r> <numDomain_c> <numAtom1> <charge1> [<numAtom2> <charge2>...]"
      write(*,*) "where the product of <numDomain_r> and <numDomain_c> must match the total number of processors." 
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end if

  !read parameters for all ranks
  call get_command_argument(1, outFilename)
  write(*,*) "outFile = ", outFilename

  call get_command_argument(2, dataFilename)
  write(*,*) "inFile.g96 = ", dataFilename

  call get_command_argument(3, numFrame_str) ! in the unit of frame number
  read(numFrame_str, *) numFrame 

  call get_command_argument(4, maxLag_str) ! in the unit of frame number
  read(maxLag_str, *) maxLag

  call get_command_argument(5, rBinWidth_str)
  read(rBinWidth_str, *) rBinWidth

  call get_command_argument(6, numDomain_r_str)
  read(numDomain_r_str, *) numDomain_r

  call get_command_argument(7, numDomain_c_str)
  read(numDomain_c_str, *) numDomain_c

  if (myrank == root) then
    if (numDomain_r * numDomain_c /= nprocs) then
      write(*,*) "error: the product of <numDomain_r> and <numDomain_c> must match the total number of processors." 
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end if

  numAtomType = num_dataArg / num_argPerData

  !rank root output parameters read
  if (myrank == root) then
    write(*,*) "numFrame= ", numFrame
    write(*,*) "maxLag = ", maxLag 
    write(*,*) "rBinWidth = ", rBinWidth
    write(*,*) "numAtomType = ", numAtomType
  end if

  allocate(numAtom(numAtomType))
  if (stat /=0) then
    write(*,*) "Allocation failed: numAtom"
    call exit(1)
  end if 

  allocate(charge(numAtomType))
  if (stat /=0) then
    write(*,*) "Allocation failed: charge"
    call exit(1)
  end if 

  do n = 1, numAtomType
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 1, numAtom_str) 
    read(numAtom_str, *) numAtom(n)
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 2, charge_str) 
    read(charge_str, *) charge(n)
  end do
  totNumAtom = sum(numAtom)

  !domain decomposition for atom pairs (numDomain_r * numDomain_c = nprocs)
  !numAtomPerDomain_r * numDomain_r ~= totNumAtom
  numDomain_c = nint(sqrt(dble(nprocs)))
  do while(mod(nprocs, numDomain_c) /= 0)
    numDomain_c = numDomain_c - 1
  end do
  numDomain_r = nprocs / numDomain_c
  numAtomPerDomain_r = totNumAtom / numDomain_r
  numAtomPerDomain_c = totNumAtom / numDomain_c
  r_start = mod(myrank, numDomain_r) * numAtomPerDomain_r + 1
  r_end = r_start + numAtomPerDomain_r - 1
  c_start = (myrank / numDomain_r) * numAtomPerDomain_c + 1
  c_end = c_start + numAtomPerDomain_c - 1

  !prepare memory for all ranks
  allocate(pos(3, numFrame, totNumAtom))
  if (stat /=0) then
    write(*,*) "Allocation failed: pos"
    call exit(1)
  end if 
  allocate(vel(3, numFrame, totNumAtom))
  if (stat /=0) then
    write(*,*) "Allocation failed: vel"
    call exit(1)
  end if 

  !read trajectory at root
  if (myrank == root) then
    allocate(pos_tmp(3, totNumAtom))
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_tmp"
      call exit(1)
    end if 
    allocate(vel_tmp(3, totNumAtom))
    if (stat /=0) then
      write(*,*) "Allocation failed: vel_tmp"
      call exit(1)
    end if 
    allocate(time(numFrame))
    if (stat /=0) then
      write(*,*) "Allocation failed: time"
      call exit(1)
    end if 

    call open_trajectory(dataFileHandle, dataFilename)
    do i = 1, numFrame
      call read_trajectory(dataFileHandle, totNumAtom, is_periodic, pos_tmp, vel_tmp, cell, time(i), stat)
      pos(:,i,:) = pos_tmp
      vel(:,i,:) = vel_tmp
    end do
    call close_trajectory(dataFileHandle)

    timestep = time(2) - time(1)
    deallocate(pos_tmp)
    deallocate(vel_tmp)
    deallocate(time)
  end if

  !distribute trajectory data collectively
  call mpi_bcast(pos, 3*numFrame*totNumAtom, mpi_double_precision, root, MPI_COMM_WORLD, ierr)
  call mpi_bcast(vel, 3*numFrame*totNumAtom, mpi_double_precision, root, MPI_COMM_WORLD, ierr)
  call mpi_bcast(cell, 3, mpi_double_precision, root, MPI_COMM_WORLD, ierr)

  num_rBin = ceiling(cell(1) / rBinWidth)
  if (myrank == root) write(*,*) "num_rBin = ", num_rBin

  !spatial decomposition correlation
  allocate(sdCorr(maxLag+1, num_rBin, numAtomType*numAtomType))
  if (stat /=0) then
    write(*,*) "Allocation failed: sdCorr"
    call exit(1)
  end if 
  sdCorr = 0

  allocate(rho(num_rBin, numAtomType*numAtomType))
  if (stat /=0) then
    write(*,*) "Allocation failed: rho"
    call exit(1)
  end if 
  rho = 0

  allocate(rBinIndex(numFrame))
  if (stat /=0) then
    write(*,*) "Allocation failed: rBinIndex"
    call exit(1)
  end if 

  do j = c_start, c_end
    do i = r_start, r_end
      if (i /= j) then
        call getBinIndex(pos(:,:,i), pos(:,:,j), cell(1), rBinWidth, rBinIndex)
        atomTypePairIndex = getAtomTypePairIndex(i, j, numAtom)
        do k = 1, maxLag+1      
          sdCorr(k, rBinIndex, atomTypePairIndex) = sdCorr(k, rBinIndex, atomTypePairIndex) + &
          & sum(vel(:, k:numFrame, i) * vel(:, 1:numFrame-k+1, j), 1)
        end do
      end if
    end do
    do k = 1, numFrame
      tmp_i = rBinIndex(k)
      rho(tmp_i, atomTypePairIndex) = rho(tmp_i, atomTypePairIndex) + 1
    end do
  end do
  deallocate(rBinIndex)

  !normalization
  allocate(norm(maxLag+1))
  if (stat /=0) then
    write(*,*) "Allocation failed: norm"
    call exit(1)
  end if 

  rho = rho / numFrame
  sdCorr = sdCorr / 3d0

  norm = [ (numFrame - (i-1), i = 1, maxLag+1) ]
  forall (i = 1:num_rBin, n = 1:numAtomType*numAtomType )
    sdCorr(:,i,n) = sdCorr(:,i,n) / norm
  end forall
  forall (i = 1:maxLag+1, n = 1:numAtomType*numAtomType )
    sdCorr(i,:,n) = sdCorr(i,:,n) / rho(:, n)
  end forall
  where (isnan(sdCorr))
    sdCorr = 0
  end where

  !output results at root
  if (myrank == root)
    call output()
  end if
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
    real(8) :: pp(size(p1,1), size(p1,2)), tmp_r
    
    pp = abs(p1 - p2)
    pp = wrap(pp, cellLength)
    rBinIndex = ceiling(sqrt(sum(pp*pp, 1)) / rBinWidth)
    where (rBinIndex == 0)
      rBinIndex = 1
    end where
  end subroutine getBinIndex

  integer function getAtomTypeIndex(i, numAtom)
    implicit none
    integer, intent(in) :: i, numAtom(:)
    integer :: n, numAtom_acc

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
    integer :: numAtomType
    integer :: ii, jj
  
    numAtomType = size(numAtom)
    ii = getAtomTypeIndex(i, numAtom)
    jj = getAtomTypeIndex(j, numAtom)
    getAtomTypePairIndex = (ii-1)*numAtomType + jj
  end function getAtomTypePairIndex

  subroutine output()
    use octave_save
    implicit none
    type(handle) :: htraj
    real(8), allocatable :: timeLags(:), rBins(:)

    allocate(timeLags(maxLag+1))
    if (stat /=0) then
      write(*,*) "Allocation failed: timeLags"
      call exit(1)
    end if 
    timeLags = [ (dble(i), i = 0, maxLag) ] * timestep
    
    allocate(rBins(num_rBin))
    if (stat /=0) then
      write(*,*) "Allocation failed: rBins"
      call exit(1)
    end if 
    rBins = [ (i - 0.5d0, i = 1, num_rBin) ] * rBinWidth

    call create_octave(htraj, outFilename)
    call write_octave_scalar(htraj, "timestep", timestep)
    call write_octave_vec(htraj, "charge", dble(charge))
    call write_octave_vec(htraj, "timeLags", timeLags)
    call write_octave_vec(htraj, "rBins", rBins)
    call write_octave_mat3(htraj, "sdCorr", sdCorr)
    call write_octave_mat2(htraj, "rho", rho)
    call close_octave(htraj)
  end subroutine output

end program spatialDecompose_mpi
