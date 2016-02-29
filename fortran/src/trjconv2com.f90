program trjconv2com
  use xdr, only : open_xdr, close_xdr, read_xdr, write_xdr, read_natom_xdr
  use top, only : open_top, close_top, read_top, system, print_sys
  implicit none
  integer, parameter :: num_parArg = 5
  integer, parameter :: num_argPerData = 2
  integer :: num_dataArg, i, j, k, n, totnummol, t, sysnumatom
  character(len=128) :: outFilename 
  character(len=128) :: dataFilename
  character(len=128) :: topFilename
  integer :: outFileHandle, dataFileHandle, topFileHandle
  integer :: numFrame, stat, numMolType, numFrameRead, percent
  integer :: moltypepair_idx, moltypepair_allidx, tmp_i, skip
  integer, allocatable :: charge(:), start_index(:)
  character(len=10) :: tmp_str
  real(8) :: cell(3), timestep, tmp_r
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :), pos_com(:, :), vel_com(:, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), dummy_null

  type(system) :: sys

  num_dataArg = command_argument_count() - num_parArg

  if (num_dataArg < num_argPerData .or. mod(num_dataArg, num_argPerData) /= 0) then
    write(*,*) "usage: $trjconv2com <outfile> <infile.trr> <topfile.top> <numFrameToRead> &
               &<skip> <molecule1> <start_index1> [<molecule2> <start_index2>...]"
    write(*,*) "Note: skip=1 means no frames are skipped. skip=2 means reading every 2nd frame."
    call exit(1)
  end if

  !read parameters for all ranks
  call get_command_argument(1, outFilename)

  call get_command_argument(2, dataFilename)

  call get_command_argument(3, topFilename)

  call get_command_argument(4, tmp_str) ! in the unit of frame number
  read(tmp_str, *) numFrame 

  call get_command_argument(5, tmp_str) ! in the unit of frame number
  read(tmp_str, *) skip

  numMolType = num_dataArg / num_argPerData

  !output parameters read
  write(*,*) "outFile = ", outFilename
  write(*,*) "inFile.trr = ", dataFilename
  write(*,*) "topFile.trr = ", topFilename
  write(*,*) "numFrame= ", numFrame
  write(*,*) "numMolType = ", numMolType

  allocate(sys%mol(numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: sys%mol"
    call exit(1)
  end if 

  allocate(charge(numMolType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: charge"
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
  write(*,*) "sys%mol%type = ", sys%mol%type
  write(*,*) "start_index = ", start_index

  !read topFile
  topFileHandle = open_top(topFilename)
  call read_top(topFileHandle, sys)
  call close_top(topFileHandle)
  call print_sys(sys)

  do n = 1, numMolType
    charge(n) = sum(sys%mol(n)%atom(:)%charge)
  end do
  totnummol = sum(sys%mol(:)%num)

  !read trajectory
  write(*,*) "start converting trajectory..."
  sysnumatom = read_natom_xdr(dataFilename)
  write(*,*) "sysnumatom=", sysnumatom

  allocate(pos_com(3, totnummol), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: pos"
    call exit(1)
  end if 
  allocate(vel_com(3, totnummol), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel"
    call exit(1)
  end if 
  allocate(pos_tmp(3, sysnumatom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: pos_tmp"
    call exit(1)
  end if 
  allocate(vel_tmp(3, sysnumatom), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel_tmp"
    call exit(1)
  end if 
  allocate(time(numFrame), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: time"
    call exit(1)
  end if 

  numFrameRead = 0
  call open_xdr(dataFileHandle, dataFilename)
  call open_xdr(outFileHandle, outFilename, 'w')
  percent = numFrame / 100
  do i = 1, numFrame
    if (mod(i, percent) == 0) write(*,*) "progress: ", i / percent, "%"
    call read_xdr(dataFileHandle, sysnumatom, pos_tmp, vel_tmp, cell, time(i), stat)
    if (stat /= 0) then
      write(*,*) "Reading trajectory error"
      call exit(1)
    end if 
    numFrameRead = numFrameRead + 1
    call com_pos(pos_com, pos_tmp, start_index, sys, cell)
    call com_vel(vel_com, vel_tmp, start_index, sys)
    call write_xdr(outFileHandle, totnummol, pos_com, vel_com, cell, i-1, time(i), stat)
    do j = 1, skip-1
      call read_xdr(dataFileHandle, sysnumatom, pos_tmp, vel_tmp, cell, tmp_r, stat)
      if (stat > 0) then
        write(*,*) "Reading trajectory error"
        call exit(1)
      else if (stat < 0) then
        !end of file
        exit
      end if 
    end do
  end do
  call close_xdr(dataFileHandle)
  call close_xdr(outFileHandle)
  write(*,*) "numFrameRead = ", numFrameRead
  if (numFrameRead /= numFrame) then
    write(*,*) "Number of frames expected to read is not the same as actually read!"
    call exit(1)
  end if

  stop

contains
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

end program trjconv2com
