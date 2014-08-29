program spatialDecompose
  use omp_lib
  use g96
  implicit none
  integer, parameter :: num_parArg = 5
  integer, parameter :: num_argPerData = 2
  integer :: num_dataArg, i, j, k, n, totNumAtom
  character(len=128) :: outFilename 
  character(len=128) :: dataFilename
  type(handle) :: dataFileHandle
  integer :: numFrame, maxLag, num_rBin, stat, numAtomType
  integer :: atomTypePairIndex, tmp_i, atomPairIndex, typei, typej, atomi, atomj
  integer, allocatable :: numAtom(:), charge(:), rBinIndex(:), norm(:)
  character(len=10) :: numFrame_str, maxLag_str, rBinWidth_str, charge_str, numAtom_str
  real(8) :: cell(3), timestep, rBinWidth
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :), vv(:)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: pos(:, :, :), vel(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), rho(:, :), rho_tmp(:, :), sdCorr(:, :, :), sdCorr_tmp(:, :, :)
  !sdCorr: spatially decomposed correlation (lag, rBin, atomTypePairIndex)
  !rho: (num_rBin, atomTypePairIndex)
  logical :: is_periodic

  !openmp
  integer :: myID, numThreads 
  real(8), allocatable :: allocationTime(:)
  real(8) :: runtime, starttime, allocationStartTime

  is_periodic = .true.

  num_dataarg = command_argument_count() - num_pararg
  if (num_dataarg < num_argperdata .or. mod(num_dataarg, num_argperdata) /= 0) then
    write(*,*) "usage: $spatialDecompose <outfile> <infile.g96> <numframe> <maxlag> <rbinwidth(nm)> &
                &<numatom1> <charge1> [<numatom2> <charge2>...]"
    call exit(1)
  end if

  call get_command_argument(1, outfilename)
  write(*,*) "outfile = ", outfilename

  call get_command_argument(2, dataFilename)
  write(*,*) "inFile.g96 = ", dataFilename

  call get_command_argument(3, numFrame_str) ! in the unit of frame number
  read(numFrame_str, *) numFrame 
  write(*,*) "numFrame= ", numFrame

  call get_command_argument(4, maxLag_str) ! in the unit of frame number
  read(maxLag_str, *) maxLag
  write(*,*) "maxLag = ", maxLag 

  call get_command_argument(5, rBinWidth_str)
  read(rBinWidth_str, *) rBinWidth
  write(*,*) "rBinWidth = ", rBinWidth
  
  numAtomType = num_dataArg / num_argPerData
  write(*,*) "numAtomType = ", numAtomType

  allocate(numAtom(numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: numAtom"
    call exit(1)
  end if 

  allocate(charge(numAtomType), stat=stat)
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

  allocate(pos_tmp(3, totNumAtom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: pos_tmp"
    call exit(1)
  end if 
  allocate(vel_tmp(3, totNumAtom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel_tmp"
    call exit(1)
  end if 
  allocate(time(numFrame), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: time"
    call exit(1)
  end if 

  allocate(pos(3, numFrame, totNumAtom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: pos"
    call exit(1)
  end if 
  allocate(vel(3, numFrame, totNumAtom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel"
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

  num_rBin = ceiling(cell(1) / rBinWidth)
  write(*,*) "num_rBin = ", num_rBin

  allocate(sdCorr(maxLag+1, num_rBin, numAtomType*numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: sdCorr"
    call exit(1)
  end if 
  sdCorr = 0d0

  allocate(rho(num_rBin, numAtomType*numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: rho"
    call exit(1)
  end if 
  rho = 0d0

  !spatial decomposition correlation
  numThreads = omp_get_max_threads()
  allocate(sdCorr_tmp(maxLag+1, num_rBin, numThreads), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: sdCorr_tmp"
    call exit(1)
  end if 

  allocate(rho_tmp(num_rBin, numThreads), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: rho_tmp"
    call exit(1)
  end if 

  allocate(allocationTime(numThreads), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: allocationTime"
    call exit(1)
  end if 
  allocationTime = 0d0
  starttime = omp_get_wtime()
  do atomTypePairIndex = 1, numAtomType*numAtomType
    sdCorr_tmp = 0d0
    rho_tmp = 0d0
    typei = (atomTypePairIndex - 1) / numAtomType + 1
    typej = mod(atomTypePairIndex - 1, numAtomType) + 1
    allocationStartTime = omp_get_wtime()
!$OMP parallel private(rBinIndex, vv, myID)
    myID = omp_get_thread_num() + 1
    allocate(rBinIndex(numFrame), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: rBinIndex"
      call exit(1)
    end if 
    allocate(vv(numFrame), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: vv"
      call exit(1)
    end if 
    allocationTime(myID) = allocationTime(myID) + omp_get_wtime() - allocationStartTime
!$OMP do private(k, n, tmp_i, atomi, atomj)
    do atomPairIndex = 1, numAtom(typei)*numAtom(typej)
      call atomPairIndex2atomIndexes(atomPairIndex, numAtom, typei, typej, atomi, atomj)
      if (atomi /= atomj) then
        call getBinIndex(pos(:,:,atomi), pos(:,:,atomj), cell(1), rBinWidth, rBinIndex)
        do k = 1, maxLag+1      
          vv = sum(vel(:, k:numFrame, atomi) * vel(:, 1:numFrame-k+1, atomj), 1)
          do n = 1, numFrame-k+1
            tmp_i = rBinIndex(n)
            sdCorr_tmp(k, tmp_i, myID) = sdCorr_tmp(k, tmp_i, myID) + vv(n)
          end do
        end do

        do k = 1, numFrame
          tmp_i = rBinIndex(k)
          rho_tmp(tmp_i, myID) = rho_tmp(tmp_i, myID) + 1d0
        end do
      end if 
    end do
!OMP end do
    deallocate(rBinIndex)
    deallocate(vv)
!$OMP end parallel
    sdCorr(:, :, atomTypePairIndex) = sum(sdCorr_tmp, 3)
    rho(:, atomTypePairIndex) = sum(rho_tmp, 2)
  end do
  runtime = omp_get_wtime() - starttime
  write(*,*) "Total time for spatialDecompose = ", runtime, " seconds"
  write(*,*) "Time spent on repeating allocation = ", sum(allocationTime), " seconds"
  
  deallocate(allocationTime)
  deallocate(sdCorr_tmp)

  deallocate(pos)
  deallocate(vel)

  !normalization
  allocate(norm(maxLag+1), stat=stat)
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
!  forall (i = 1:maxLag+1, n = 1:numAtomType*numAtomType )
!    where (rho > 0d0)
!      sdCorr(i,:,n) = sdCorr(i,:,n) / rho(:, n)
!    end where
!  end forall

  deallocate(norm)

  do n = 1, numAtomType*numAtomType
    do j = 1, num_rBin
      if (rho(j, n) > 0d0) then
          sdCorr(:,j,n) = sdCorr(:,j,n) / rho(j, n)
      else if (any(abs(sdCorr(:,j,n)) > 0d0)) then
        !rho is zero, sdCorr should also be zero
        write(*,*) "error: rho(",j,",",n,") = ", rho(j,n)
        write(*,*) "but sdCorr(:,", j, ",", n, ") are not all zero"
        call exit(1)
      end if
    end do
  end do

!  where (isnan(sdCorr))
!    sdCorr = 0
!  end where

  !output results
  call output()
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

    allocate(timeLags(maxLag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: timeLags"
      call exit(1)
    end if 
    timeLags = [ (dble(i), i = 0, maxLag) ] * timestep
    
    allocate(rBins(num_rBin), stat=stat)
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

    deallocate(timeLags)
    deallocate(rBins)
  end subroutine output
  
  pure subroutine atomPairIndex2atomIndexes(atomPairIndex, numAtom, typei, typej, atomi, atomj)
    implicit none
    integer, intent(in) :: atomPairIndex, numAtom(:), typei, typej
    integer, intent(out) :: atomi, atomj

    atomi = sum(numAtom(1:typei-1)) + (atomPairIndex - 1) / numAtom(typej) + 1
    atomj = sum(numAtom(1:typej-1)) + mod(atomPairIndex-1, numAtom(typej)) + 1
  end subroutine atomPairIndex2atomIndexes

end program spatialDecompose
