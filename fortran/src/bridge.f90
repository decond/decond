program bridge
  use trr
  implicit none
  integer, parameter :: num_parArg = 4
  integer, parameter :: num_argPerData = 1
  integer :: num_dataArg, i, j, k, n, totNumAtom, t
  character(len=128) :: outFilename 
  character(len=128) :: dataFilename
  type(handle) :: dataFileHandle
  integer :: numFrame, stat, numAtomType, numFrameRead
  integer :: skip
  integer, allocatable :: numAtom(:)
  character(len=10) :: tmp_str
  real(8) :: cell(3), timestep, numer(2), denom(2), probability(2), tmp_r
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: pos(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:)

  !first coordination shell radius
  real(8), parameter :: lambda_na = 0.575
  real(8), parameter :: lambda_nacl = 0.375
  real(8), parameter :: lambda_cl = 0.615

  !sdCorr: spatially decomposed correlation (lag, rBin, atomTypePairIndex)
  !rho: (num_rBin, atomTypePairIndex)
  logical :: is_periodic
  is_periodic = .true.

  !Reading arguments and initialization
  num_dataarg = command_argument_count() - num_pararg
  if (num_dataarg < num_argperdata .or. mod(num_dataarg, num_argperdata) /= 0) then
    write(*,*) "Usage: $bridge <outfile> <infile.trr> <numFrameToRead> <skip> &
                &<numatom1> [<numatom2> ...]"
    write(*,*) "Note: skip=1 means no frames are skipped. skip=2 means reading every 2nd frame."
    call exit(1)
  end if

  call get_command_argument(1, outFilename)
  write(*,*) "outfile = ", outFilename

  call get_command_argument(2, dataFilename)
  write(*,*) "inFile.trr = ", dataFilename

  call get_command_argument(3, tmp_str) ! in the unit of frame number
  read(tmp_str, *) numFrame 
  write(*,*) "numFrame = ", numFrame

  call get_command_argument(4, tmp_str) ! in the unit of frame number
  read(tmp_str, *) skip
  write(*,*) "skip = ", skip

  
  numAtomType = num_dataArg / num_argPerData
  write(*,*) "numAtomType = ", numAtomType

  allocate(numAtom(numAtomType), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: numAtom"
    call exit(1)
  end if 

  do n = 1, numAtomType
    call get_command_argument(num_parArg + num_argPerData*(n-1) + 1, tmp_str) 
    read(tmp_str, *) numAtom(n)
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

  denom = 0
  numer = 0

  ! Na-Na bridge by Cl
  do t = 1, numFrameRead
    do i = 1, numAtom(1)
      do j = i + 1, numAtom(1)
        if (distance(pos(:,t,i), pos(:,t,j)) < lambda_na) then
          denom(1) = denom(1) + 1
          do k = numAtom(1) + 1, numAtom(1) + numAtom(2)
            if( distance(pos(:,t,i), pos(:,t,k)) < lambda_nacl .and. &
             &distance(pos(:,t,j), pos(:,t,k)) < lambda_nacl ) then
              numer(1) = numer(1) + 1
              exit
            endif
          enddo
        endif
      enddo
    enddo
  enddo

  ! Cl-Cl bridge by Na
  do t = 1, numFrameRead
    do i = numAtom(1) + 1, numAtom(1) + numAtom(2)
      do j = i + 1, numAtom(1) + numAtom(2)
        if (distance(pos(:,t,i), pos(:,t,j)) < lambda_cl ) then
          denom(2) = denom(2) + 1
          do k = 1, numAtom(1)
            if( distance(pos(:,t,i), pos(:,t,k)) < lambda_nacl .and. &
             &distance(pos(:,t,j), pos(:,t,k)) < lambda_nacl ) then
              numer(2) = numer(2) + 1
              exit
            endif
          enddo
        endif
      enddo
    enddo
  enddo

  probability = numer / denom

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
  
  real(8) function distance(p1, p2)
    implicit none
    real(8), intent(in) :: p1(:), p2(:)
    !p1(dim)
    real(8) :: pp(size(p1)), cellLength
    
    cellLength = cell(1)
    pp = abs(p1 - p2)
    pp = wrap(pp, cellLength)
    distance = sqrt(sum(pp*pp))
  end function distance

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
    real(8) :: lambda(3)

    lambda(1) = lambda_na
    lambda(2) = lambda_nacl
    lambda(3) = lambda_cl

    call create_octave(htraj, outFilename)
!    call write_octave_scalar(htraj, "timestep", timestep)
    call write_octave_vec(htraj, "lambda", lambda)
    call write_octave_vec(htraj, "bridge_probability", probability)
!    call write_octave_mat3(htraj, "sdCorr", sdCorr)
!    call write_octave_mat2(htraj, "rho", rho)
    call close_octave(htraj)
  end subroutine output

end program bridge
