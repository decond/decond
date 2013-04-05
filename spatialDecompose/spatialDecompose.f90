program spatialDecompose
  use trajectory
  integer, parameter :: num_parArg = 4
  integer, parameter :: num_argPerData = 3
  integer :: num_dataArg, i, j, n, totNumAtom
  character(len=128) :: outFilename 
  character(len=128), allocatable :: dataFilename(:)
  type(handle) :: dataFileHandle
  integer :: numFrame, maxLag, rBinWidth, num_rBin, stat, numAtomType
  integer :: atomTypePairIndex
  integer, allocatable :: numAtom(:),, charge(:), rBinIndex(:)
  character(len=10) :: numFrame_str, maxLag_str, rBinWidth_str
  real(8) :: cell(3), timestep
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: pos(:, :, :), vel(:, :, :)
  !pos(dim=3, timeFrame, atomType), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), rho(:, :), sdCorr(:, :, :)
  !sdCorr: spatially decomposed correlation (lag, rBin, atomTypePairIndex)
  !rho: (num_rBin, atomTypePairIndex)

  is_periodic = .true.

  num_dataArg = command_argument_count() - num_parArg
  if (num_dataArg < num_argPerData || mod(num_dataArg, num_argPerData) /= 0) then
    write(*,*) "Usage: $spatialDecompose <outFilename> <numFrame> <maxLag> <rBinWidth(nm)> &
                &<data1> <charge1> [<data2> <charge2>...]"
    call exit(1)
  end if

  call get_command_argument(1, outFilename)
  write(*,*) "outFilename = ", outFilename

  call get_command_argument(2, numFrame_str) ! in the unit of frame number
  read(numFrame_str, "(I)") numFrame 
  write(*,*) "numFrame= ", numFrame

  call get_command_argument(3, maxLag_str) ! in the unit of frame number
  read(maxLag_str, "(I)") maxLag
  write(*,*) "maxLag = ", maxLag 

  call get_command_argument(4, rBinWidth_str)
  read(rBinWidth_str, "(I)") rBinWidth
  write(*,*) "rBinWidth = ", rBinWidth
  
  numAtomType = num_dataArg / num_argPerData
  write(*,*) "numAtomType = ", numAtomType

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
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 1, dataFilename) 
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 2, charge(i)) 
    numAtom(n) = get_natom(dataFilename)
  end do
  totNumAtom = sum(numAtom)

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

  allocate(rBinIndex(numFrame))
  if (stat /=0) then
    write(*,*) "Allocation failed: rBinIndex"
    call exit(1)
  end if 

  allocate(pos(3, numFrame, totNumAtom))
  if (stat /=0) then
    write(*,*) "Allocation failed: pos"
    call exit(1)
  end if 
  allocate(vel(3, numFrame, numAtomType))
  if (stat /=0) then
    write(*,*) "Allocation failed: vel"
    call exit(1)
  end if 

  do n = 1, numAtomType
    call open_trajectory(dataFileHandle, dataFilename)
    do i = 1, numFrame
      call read_trajectory(dataFileHandle, numAtom(n), is_periodic, pos_tmp, vel_tmp, cell, time(n), stat)
      pos(:,i,:) = pos_tmp
      vel(:,i,:) = vel_tmp
    end do
    call close_trajectory(dataFileHandle)
  end do

  timestep = time(2) - time(1)
  deallocate(pos_tmp)
  deallocate(vel_tmp)
  deallocate(time)

  call read_cell(dataFilename, cell, stat)
  if (stat /= 0) then
    stop "failed reading cell info from file: ", trim(dataFilename)
  end if

  num_rBin = ceiling(cell(1) / rBinWidth)
  write(*,*) "num_rBin = ", num_rBin

  allocate(sdCorr(maxLag+1, num_rBin, numAtomType*numAtomType))
  if (stat /=0) then
    write(*,*) "Allocation failed: sdCorr"
    call exit(1)
  end if 
  allocate(rho(num_rBin, numAtomType*numAtomType))
  if (stat /=0) then
    write(*,*) "Allocation failed: rho"
    call exit(1)
  end if 

  do i = totNumAtom
    do j = totNumAtom
      call getBinIndex(pos(:,:,i), pos(:,:,j), cell(1), rBinWidth, rBinIndex)
      atomTypePairIndex = getAtomTypePairIndex(i, j, numAtom)
      do k = 1, maxLag+1      
        sdCorr(k, rBinIndex, atomTypePairIndex) = sdCorr(k, rBinIndex, atomTypePairIndex) + &
        & sum(vel(:, k:numFrame, i) * vel(:, 1:numFrame-k+1, j), 1)
        rho(rBinIndex, atomTypePairIndex) = rho(rBinIndex, atomTypePairIndex) + 1
      end do
    end do
  end do

  do n = 1, numAtomType*numAtomType
    sdCorr = sdCorr / (3*
  end do

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
  
  pure subroutine getBinIndex(p1, p2, cellLength, rBinWidth, rBinIndex)
    implicit none
    real(8), intent(in) :: p1(3,:), p2(3,:), cellLength, rBinWidth
    !p1(dim,timeFrame)
    integer, intent(out) :: rBinIndex(size(p1, 2))
    real(8) :: pp(3, size(p1,2))
    
    pp = abs(p1 - p2)
    pp = wrap(pp, cellLength)
    rBinIndex = sqrt(sum(pp*pp, 1))
  end subroutine getBinIndex

  integer function getAtomTypeIndex(i, numAtom)
    implicit none
    integer, intent(in) :: i
    real(8), intent(in) :: numAtom(:)
    integer :: numAtomType = size(numAtom)
    integer :: n, numAtom_acc
    
    numAtom_acc = 0
    do n = 1, numAtomType
      numAtom_acc = numAtom_acc + numAtom(n)
      if (i <= numAtom_acc) then
        getAtomTypeIndex = n
      end if
    end do
  end function getAtomTypeIndex
  
  integer function getAtomTypePairIndex(i, j, numAtom)
    implicit none
    integer, intent(in) :: i, j
    real(8), intent(in) :: numAtom(:)
    integer :: numAtomType = size(numAtom)
    integer :: ii, jj
  
    ii = getAtomTypeIndex(i, numAtom)
    jj = getAtomTypeIndex(j, numAtom)
    getAtomTypePairIndex = (ii-1)*numAtomType + jj
  end function getAtomTypePairIndex

dataIndex{1} = [1:numAtom(1)];
for n = [2:numAtomType]
    dataIndex{n} = [dataIndex{n-1}(end) + 1: dataIndex{n-1}(end) + numAtom(n)];
endfor

function serialIndex = atomPair2SerialIndex(idx1, idx2)
    global totalNumAtoms;
    serialIndex = (idx1 - 1) * totalNumAtoms + idx2;
endfunction

function groupIndex = atomIndex2GroupIndex(idx)
    global dataIndex;
    for i = [1:length(dataIndex)]
        if (any(dataIndex{i} == idx))
            groupIndex = i;
            return;
        endif
    endfor
    error(strcat("Unable to convert atomIndex:", num2str(idx), " to groupIndex"));
endfunction

function wrappedR = wrap(r)
    global boxLength;
    wrappedR = r;
    isOutsideHalfBox = (r > (boxLength / 2));
    wrappedR(isOutsideHalfBox) = abs(r(isOutsideHalfBox) - boxLength);
endfunction

function rBinIndex = getBinIndex(pos1, pos2)
    global rBinWidth;
    r = abs(pos1 - pos2);
    r = wrap(r);
    r = sqrt(sum(r.*r, 2));
    rBinIndex = ceil(r ./ rBinWidth);
    rBinIndex(rBinIndex == 0) = 1;
endfunction

#average 3 dimensions
for i = [1:numAtomType]
%    cAutocorr{i} ./= (3*rhoAutocorr{i});
    for j = [1:numAtomType]
        cCorr{i,j} ./= 3*repmat([num_frames:-1:num_frames-maxLag]', [1, num_rBin]).*repmat(rhoCorr{i,j}', [maxLag+1, 1]);
    endfor
endfor

puts("Tag5\n");
whos


# output time vector for convenience of plotting
timeLags = [0:maxLag]' * timestep;
rBins = ([1:num_rBin] - 0.5)* rBinWidth;

save(strcat(outFilename, ".cCorr"), "timestep", "charge", "numAtom", "timeLags", "rBins", "cCorr", "rhoCorr");
#save(strcat(outFilename, ".cCorrUnaverage"), "timestep", "charge", "numAtom", "timeLags", "rBins", "cCorr", "rhoCorr");

stop
end program spatialDecompose
