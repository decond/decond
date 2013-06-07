program spatialDecompose
  use g96
  implicit none
  integer, parameter :: num_parArg = 6
  integer, parameter :: num_argPerData = 2
  integer :: num_dataArg, i, j, k, n, totNumAtom
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
  real(8), allocatable :: time(:), rho(:, :), sdCorr(:, :, :)
  integer(kind=4), allocatable :: sdCorrNorm(:, :, :)
  !sdCorr: spatially decomposed correlation (lag, rBin, atomTypePairIndex)
  !rho: (num_rBin, atomTypePairIndex)
  logical :: is_periodic

  is_periodic = .true.

  !Reading arguments and initialization
  num_dataarg = command_argument_count() - num_pararg
  if (num_dataarg < num_argperdata .or. mod(num_dataarg, num_argperdata) /= 0) then
    write(*,*) "Usage: $spatialdecompose <outfile> <infile.g96> <numFrameToRead> <skip> <maxlag> <rbinwidth(nm)> &
                &<numatom1> <charge1> [<numatom2> <charge2>...]"
    write(*,*) "Note: skip=1 means no frames are skipped. skip=2 means reading every 2nd frame."
    write(*,*) "Note: maxlag is counted in terms of the numFrameToRead."
    call exit(1)
  end if

  call get_command_argument(1, outFilename)
  write(*,*) "outfile = ", outFilename

  call get_command_argument(2, dataFilename)
  write(*,*) "inFile.g96 = ", dataFilename

  call get_command_argument(3, tmp_str) ! in the unit of frame number
  read(tmp_str, *) numFrame 
  write(*,*) "numFrame = ", numFrame

  call get_command_argument(4, tmp_str) ! in the unit of frame number
  read(tmp_str, *) skip
  write(*,*) "skip = ", skip

  call get_command_argument(5, tmp_str) ! in the unit of frame number
  read(tmp_str, *) maxLag
  write(*,*) "maxLag = ", maxLag 

  call get_command_argument(6, tmp_str)
  read(tmp_str, *) rBinWidth
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
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 1, tmp_str) 
    read(tmp_str, *) numAtom(n)
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 2, tmp_str) 
    read(tmp_str, *) charge(n)
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

  numFrameRead = 0
  call open_trajectory(dataFileHandle, dataFilename)
  do i = 1, numFrame
    call read_trajectory(dataFileHandle, totNumAtom, is_periodic, pos_tmp, vel_tmp, cell, time(i), stat)
    if (stat /=0) then
      write(*,*) "Reading trajectory error"
      call exit(1)
    end if 
    numFrameRead = numFrameRead + 1
    pos(:,i,:) = pos_tmp
    vel(:,i,:) = vel_tmp
    do j = 1, skip-1
      call read_trajectory(dataFileHandle, totNumAtom, is_periodic, pos_tmp, vel_tmp, cell, tmp_r, stat)
      if (stat > 0) then
        write(*,*) "Reading trajectory error"
        call exit(1)
      else if (stat < 0) then
        !end of file
        exit
      end if 
    end do
  end do
  call close_trajectory(dataFileHandle)
  write(*,*) "numFrameRead = ", numFrameRead

  timestep = time(2) - time(1)
  deallocate(pos_tmp)
  deallocate(vel_tmp)
  deallocate(time)

  num_rBin = ceiling(cell(1) / 2d0 / rBinWidth)
  write(*,*) "num_rBin = ", num_rBin

  allocate(sdCorr(maxLag+1, num_rBin, numAtomType*numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: sdCorr"
    call exit(1)
  end if 
  sdCorr = 0d0

  allocate(sdCorrNorm(maxLag+1, num_rBin, numAtomType*numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: sdCorrNorm"
    call exit(1)
  end if 
  sdCorrNorm = 0

  allocate(rho(num_rBin, numAtomType*numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: rho"
    call exit(1)
  end if 
  rho = 0d0

  allocate(rBinIndex(numFrameRead), stat=stat)
  if (stat /= 0) then
    write(*,*) "Allocation failed: rBinIndex"
    call exit(1)
  end if 

  allocate(vv(numFrameRead))
  if (stat /=0) then
    write(*,*) "Allocation failed: vv"
    call exit(1)
  end if 

  !spatial decomposition correlation
  do i = 1, totNumAtom
    do j = 1, totNumAtom
      if (i /= j) then
        !get rBinIndex
        call getBinIndex(pos(:,:,i), pos(:,:,j), cell(1), rBinWidth, rBinIndex)
        !get the index for different atomType pair (ex. Na-Na, Na-Cl, Cl-Na, Cl-Cl)
        atomTypePairIndex = getAtomTypePairIndex(i, j, numAtom)
        do k = 1, maxLag+1      
          vv = sum(vel(:, k:numFrameRead, i) * vel(:, 1:numFrameRead-k+1, j), 1)
          aveNum = count(rBinIndex(1:numFrameRead-k+1) > 0)
          sdCorr_tmp = 0d0
          do n = 1, numFrameRead-k+1
            tmp_i = rBinIndex(n)
            if (tmp_i > 0)
!              sdCorr(k, tmp_i, atomTypePairIndex) = sdCorr(k, tmp_i, atomTypePairIndex) + vv(n)
              sdCorr_tmp(tmp_i) = sdCorr_tmp(tmp_i) + vv(n)
!              sdCorrNorm(k, tmp_i, atomTypePairIndex) = sdCorrNorm(k, tmp_i, atomTypePairIndex) + 1
            end if
          end do
          sdCorr(k, :, atomTypePairIndex) = sdCorr(k, :, atomTypePairIndex) + (sdCorr_tmp / dble(aveNum))
        end do
        rhoCount_tmp = 0
        aveNum = count(rBinIndex > 0)
        do t = 1, numFrameRead
          tmp_i = rBinIndex(t)
          if (tmp_i > 0)
            rhoCount_tmp(tmp_i) = rhoCount_tmp(tmp_i) + 1
          end if
        end do
        rho(:, atomTypePairIndex) = rho(:, atomTypePairIndex) + (dble(rhoCount_tmp) / dble(aveNum))
      end if
    end do
  end do
  deallocate(rBinIndex)
  deallocate(pos)
  deallocate(vel)

  !normalization
!  allocate(norm(maxLag+1), stat=stat)
!  if (stat /=0) then
!    write(*,*) "Allocation failed: norm"
!    call exit(1)
!  end if 

  rho = rho / numFrameRead
  sdCorr = sdCorr / 3d0

  norm = [ (numFrameRead - (i-1), i = 1, maxLag+1) ]
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
    where (rBinIndex > ceiling(cellLength / 2.d0 / rBinWidth))
      rBinIndex = -1
    end where
  end subroutine getBinIndex

  integer function getAtomTypeIndex(i, numAtom)
    implicit none
    integer, intent(in) :: i, numAtom(:)
    integer :: n, numAtom_acc
!    integer, save :: numAtomType = size(numAtom)

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
  end subroutine output

end program spatialDecompose
