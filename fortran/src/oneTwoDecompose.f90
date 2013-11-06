program oneTwoDecompose
  use trr
  use correlation
  implicit none
  integer, parameter :: num_parArg = 5
  integer, parameter :: num_argPerData = 2
  integer :: num_dataArg, i, j, k, n, totNumAtom, t
  character(len=128) :: outFilename 
  character(len=128) :: dataFilename
  type(handle) :: dataFileHandle
  integer :: numFrame, maxLag, stat, numAtomType, numFrameRead
  integer :: idx1, idx2, tmp_i, skip
  integer, allocatable :: numAtom(:), charge(:), norm(:)
  character(len=10) :: tmp_str
  real(8) :: cell(3), timestep, tmp_r
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: vel(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), autoCorr(:,:), crossCorr(:,:), corr_tmp(:)
  logical :: is_periodic

  is_periodic = .true.

  !Reading arguments and initialization
  num_dataarg = command_argument_count() - num_pararg
  if (num_dataarg < num_argperdata .or. mod(num_dataarg, num_argperdata) /= 0) then
    write(*,*) "Usage: $oneTwoDecompose <outfile> <infile.trr> <numFrameToRead> <skip> <maxLag (-1=max)> &
                &<numatom1> <charge1> [<numatom2> <charge2>...]"
    write(*,*) "Note: skip=1 means no frames are skipped. skip=2 means reading every 2nd frame."
    write(*,*) "Note: maxLag is counted in terms of the numFrameToRead."
    call exit(1)
  end if

  call get_command_argument(1, outFilename)
  write(*,*) "outfile = ", outFilename

  call get_command_argument(2, dataFilename)
  write(*,*) "inFile.trr= ", dataFilename

  call get_command_argument(3, tmp_str) ! in the unit of frame number
  read(tmp_str, *) numFrame 
  write(*,*) "numFrame = ", numFrame

  call get_command_argument(4, tmp_str) ! in the unit of frame number
  read(tmp_str, *) skip
  write(*,*) "skip = ", skip

  call get_command_argument(5, tmp_str) ! in the unit of frame number
  read(tmp_str, *) maxLag
  if (maxLag < 0) then
      maxLag = numFrame - 1;
  endif
  write(*,*) "maxLag = ", maxLag 

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

  allocate(autoCorr(maxLag+1, numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: autoCorr"
    call exit(1)
  end if 
  autoCorr = 0d0

  allocate(crossCorr(maxLag+1, numAtomType*numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: crossCorr"
    call exit(1)
  end if 
  crossCorr = 0d0

  allocate(corr_tmp(2*maxLag+1), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: corr_tmp"
    call exit(1)
  end if 
  corr_tmp = 0d0

  !autoCorr part
  do i = 1, totNumAtom
write(*,*) "autoCorr: ", i
    idx1 = getAtomTypeIndex(i, numAtom)
    do k = 1, 3
      corr_tmp = corr(vel(k, :, i), maxLag)
      autoCorr(:, idx1) = autoCorr(:, idx1) + corr_tmp(maxLag+1:)
    end do
  end do

  !crossCorr part
  do i = 1, totNumAtom
    do j = i+1, totNumAtom
write(*,*) "crossCorr: ", i, j
      !get the index for different atomType pair (ex. Na-Na, Na-Cl, Cl-Na, Cl-Cl)
      idx1 = getAtomTypePairIndex(i, j, numAtom)
      idx2 = getAtomTypePairIndex(j, i, numAtom)
      do k = 1, 3
        corr_tmp = corr(vel(k, :, i), vel(k, :, j), maxLag)
        crossCorr(:, idx1) = crossCorr(:, idx1) + corr_tmp(maxLag+1:)
        crossCorr(:, idx2) = crossCorr(:, idx2) + corr_tmp(maxLag+1:1:-1)
      end do
    end do
  end do
  deallocate(vel)

  !normalization
  autoCorr = autoCorr / 3d0
  crossCorr = crossCorr / 3d0

  allocate(norm(maxLag+1), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: norm"
    call exit(1)
  end if 
  norm = [ (numFrameRead - (i-1), i = 1, maxLag+1) ]

  do i = 1, numAtomType
      autoCorr(:, i) = autoCorr(:, i) / norm
  end do
  do i = 1, numAtomType*numAtomType
    crossCorr(:, i) = crossCorr(:, i) / norm
  end do
  deallocate(norm)

  !output results
  call output()
  stop

contains
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
      call exit(1)
    end if 
    timeLags = [ (dble(i), i = 0, maxLag) ] * timestep
    
    call create_octave(htraj, outFilename)
    call write_octave_scalar(htraj, "timestep", timestep)
    call write_octave_vec(htraj, "charge", dble(charge))
    call write_octave_vec(htraj, "cell", cell)
    call write_octave_vec(htraj, "timeLags", timeLags)
    call write_octave_mat2(htraj, "autoCorr", autoCorr)
    call write_octave_mat2(htraj, "crossCorr", crossCorr)
    call close_octave(htraj)
  end subroutine output

end program oneTwoDecompose
