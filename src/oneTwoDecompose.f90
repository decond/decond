program oneTwoDecompose
  use utility, only : handle
  use xdr, only : open_trajectory, close_trajectory, read_trajectory, get_natom
  use top, only : open_top, close_top, read_top, system
  use correlation
  implicit none
  integer, parameter :: num_parArg = 6
  integer, parameter :: num_argPerData = 2
  integer :: num_dataArg, i, j, k, n, totNumMol, t, sysNumAtom
  character(len=128) :: outFilename 
  character(len=128) :: dataFilename
  character(len=128) :: topFilename
  type(handle) :: dataFileHandle, topFileHandle
  integer :: numFrame, maxLag, stat, numMolType, numFrameRead
  integer :: idx1, idx2, tmp_i, skip, one_percent
  integer, allocatable :: charge(:), norm(:), start_index(:)
  character(len=10) :: tmp_str
  real(8) :: cell(3), timestep, tmp_r
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: vel(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), autoCorr(:,:), crossCorr(:,:), corr_tmp(:)
  logical :: is_periodic
  type(system) :: sys

  is_periodic = .true.

  !Reading arguments and initialization
  num_dataarg = command_argument_count() - num_pararg
  if (num_dataarg < num_argperdata .or. mod(num_dataarg, num_argperdata) /= 0) then
    write(*,*) "Usage: $oneTwoDecompose <outfile> <infile.trr> <topfile.top> <numFrameToRead>&
                <skip> <maxLag (-1=max)> <molecule1> <start_index> [<molecule2> <start_index>...]"
    write(*,*) "Note: start_index is the 1st index of a moleculetype in the trajectory file"
    write(*,*) "Note: skip=1 means no frames are skipped. skip=2 means reading every 2nd frame."
    write(*,*) "Note: maxLag is counted in terms of the numFrameToRead."
    call exit(1)
  end if

  call get_command_argument(1, outFilename)
  write(*,*) "outfile = ", outFilename

  call get_command_argument(2, dataFilename)
  write(*,*) "inFile.trr= ", dataFilename

  call get_command_argument(3, topFilename)
  write(*,*) "topFilename= ", topFilename

  call get_command_argument(4, tmp_str) ! in the unit of frame number
  read(tmp_str, *) numFrame 
  write(*,*) "numFrame = ", numFrame

  call get_command_argument(5, tmp_str) ! in the unit of frame number
  read(tmp_str, *) skip
  write(*,*) "skip = ", skip

  call get_command_argument(6, tmp_str) ! in the unit of frame number
  read(tmp_str, *) maxLag
  if (maxLag < 0) then
      maxLag = numFrame - 1;
  endif
  write(*,*) "maxLag = ", maxLag 

  numMolType = num_dataArg / num_argPerData
  write(*,*) "numMolType = ", numMolType

  allocate(charge(numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: charge"
    call exit(1)
  end if 

  allocate(sys%mol(numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: sys%mol"
    call exit(1)
  end if 

  allocate(start_index(numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: start_index"
    call exit(1)
  end if 

  do n = 1, numMolType
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 1, sys%mol(n)%type)
    write(*,*) "molName= ", sys%mol(n)%type
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 2, tmp_str) 
    read(tmp_str, *) start_index(n)
  end do

  topFileHandle = open_top(topFilename)
  call read_top(topFileHandle, sys)
  call close_top(topFileHandle)

  do n = 1, numMolType
    charge(n) = sum(sys%mol(n)%atom(:)%charge)
  end do
  totNumMol = sum(sys%mol(:)%num)

  sysNumAtom = get_natom(dataFilename)
  write(*,*) "sysNumAtom=", sysNumAtom

  allocate(pos_tmp(3, sysNumAtom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: pos_tmp"
    call exit(1)
  end if 
  allocate(vel_tmp(3, sysNumAtom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel_tmp"
    call exit(1)
  end if 
  allocate(time(numFrame), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: time"
    call exit(1)
  end if 
  allocate(vel(3, numFrame, totNumMol), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel"
    call exit(1)
  end if 

  numFrameRead = 0
  call open_trajectory(dataFileHandle, dataFilename)
  one_percent = nint(numFrame / 100d0)
  do i = 1, numFrame
    if (mod(i, one_percent) == 0) write(*,*) "Reading... ", floor(dble(i) / one_percent), "%"
    call read_trajectory(dataFileHandle, sysNumAtom, is_periodic, pos_tmp, vel_tmp, cell, time(i), stat)
    if (stat /=0) then
      write(*,*) "Reading trajectory error"
      call exit(1)
    end if 
    numFrameRead = numFrameRead + 1
    call com_vel(vel(:, i, :), vel_tmp, start_index, sys)
    do j = 1, skip-1
      call read_trajectory(dataFileHandle, sysNumAtom, is_periodic, pos_tmp, vel_tmp, cell, tmp_r, stat)
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

  allocate(autoCorr(maxLag+1, numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: autoCorr"
    call exit(1)
  end if 
  autoCorr = 0d0

  allocate(crossCorr(maxLag+1, numMolType*numMolType), stat=stat)
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
  do i = 1, totNumMol
write(*,*) "autoCorr: ", i
    idx1 = getMolTypeIndex(i, sys%mol(:)%num)
    do k = 1, 3
      corr_tmp = corr(vel(k, :, i), maxLag)
      autoCorr(:, idx1) = autoCorr(:, idx1) + corr_tmp(maxLag+1:)
    end do
  end do

  !crossCorr part
  do i = 1, totNumMol
    do j = i+1, totNumMol
write(*,*) "crossCorr: ", i, j
      !get the index for different atomType pair (ex. Na-Na, Na-Cl, Cl-Na, Cl-Cl)
      idx1 = getMolTypePairIndex(i, j, sys%mol(:)%num)
      idx2 = getMolTypePairIndex(j, i, sys%mol(:)%num)
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

  do i = 1, numMolType
      autoCorr(:, i) = autoCorr(:, i) / norm
  end do
  do i = 1, numMolType*numMolType
    crossCorr(:, i) = crossCorr(:, i) / norm
  end do
  deallocate(norm)

  !output results
  call output()
  stop

contains
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

  subroutine output()
    use H5DS
    use H5LT
    use HDF5
    implicit none
    real(8), allocatable :: timeLags(:)
    integer :: ierr
    integer(hid_t) :: fid, did1, did2, sid

    allocate(timeLags(maxLag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: timeLags"
      call exit(1)
    end if 
    timeLags = [ (dble(i), i = 0, maxLag) ] * timestep
    
    call H5open_f(ierr)

    ! create a HDF5 file 
    call H5Fcreate_f(outFilename, H5F_ACC_TRUNC_F, fid, ierr)

    ! create and write dataset
    call H5LTset_attribute_double_f(fid, "/", "timestep", [timestep], int(1, kind=size_t), ierr)
    call H5LTset_attribute_int_f(fid, "/", "charge", charge, size(charge, kind=size_t), ierr)
    call H5LTset_attribute_int_f(fid, "/", "numMol", sys%mol(:)%num, size(sys%mol(:)%num, kind=size_t), ierr)
    call H5LTset_attribute_double_f(fid, "/", "cell", cell, size(cell, kind=size_t), ierr)


    call H5LTmake_dataset_double_f(fid, "autoCorr", 2, &
        [size(autoCorr, 1, kind=hsize_t), size(autoCorr, 2, kind=hsize_t)], autoCorr, ierr)
    call H5Dopen_f(fid, "autoCorr", did1, ierr)

    call H5LTmake_dataset_double_f(fid, "crossCorr", 2, &
        [size(crossCorr, 1, kind=hsize_t), size(crossCorr, 2, kind=hsize_t)], crossCorr, ierr)
    call H5Dopen_f(fid, "crossCorr", did2, ierr)


    call H5LTmake_dataset_double_f(fid, "timeLags", 1, [size(timeLags, kind=hsize_t)], timeLags, ierr)
    call H5Dopen_f(fid, "timeLags", sid, ierr)

    ! attach to 1st dimension
    call H5DSattach_scale_f(did1, sid, 1, ierr)
    call H5DSattach_scale_f(did2, sid, 1, ierr)

    call H5Dclose_f(sid, ierr)
    call H5Dclose_f(did1, ierr)
    call H5Dclose_f(did2, ierr)
    call H5Fclose_f(fid, ierr)
    call H5close_f(ierr)
  end subroutine output

end program oneTwoDecompose
