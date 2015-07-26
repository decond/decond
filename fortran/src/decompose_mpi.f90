program decompose_mpi
  use mpiproc
  use utility, only : handle
  use xdr, only : open_trajectory, close_trajectory, read_trajectory, get_natom
  use top, only : open_top, close_top, read_top, system, print_sys
  use correlation
  use spatial_dec, only : sd_init, rBinWidth, pos_r, pos_c, rho, sdCorr, pos, num_rBin, &
                          com_pos, sd_binIndex, sd_prepCorrMemory, sd_getBinIndex, &
                          sd_cal_num_rBin, sd_broadcastPos, sd_prepPosMemory, &
                          sd_collectCorr, sd_normalize
  use engdec, only : engtrjFilename

  implicit none
  integer, parameter :: NUM_POSITIONAL_ARG = 2, LEAST_REQUIRED_NUM_ARG = 6
  integer :: num_arg, num_subArg, num_argPerMolType
  integer :: i, j, k, n, totNumMol, t, sysNumAtom
  character(len=128) :: outCorrFilename, dataFilename, topFilename, arg
  type(handle) :: dataFileHandle, topFileHandle
  integer :: numFrame, maxLag, stat, numMolType, numFrameRead, numFrame_k
  integer :: molTypePairIndex, molTypePairAllIndex, tmp_i, skip
  integer, allocatable :: charge(:), norm(:), start_index(:)
  real(8) :: cell(3), timestep, tmp_r
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :), vv(:)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: vel(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), nCorr(:, :), corr_tmp(:)
  !sdCorr: spatially decomposed correlation (lag, rBin, molTypePairIndex)
  !rho: (num_rBin, molTypePairIndex)
  logical :: is_periodic, is_pa_mode, is_pm_mode, is_sd, is_ed
  type(system) :: sys
  integer(hid_t) :: outCorrFileid

  !MPI variables
  real(8), allocatable :: vel_r(:, :, :), vel_c(:, :, :)

  !initialize
  call mpi_setup('init')

  prog_starttime = MPI_Wtime()

  num_arg = command_argument_count()
  is_pa_mode = .false.
  is_pm_mode = .false.

  !default values
  outCorrFilename = 'corr.h5'
  skip = 1
  maxLag = -1
  is_sd = .false.
  is_ed = .false.
  numDomain_r = 0
  numDomain_c = 0
  call sd_init()
  call ed_init()

  !root checks the number of the input arguments
  is_periodic = .true.
  if (num_arg < LEAST_REQUIRED_NUM_ARG) then
    if (myrank == root) call print_usage()
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if

  !read parameters for all ranks
  i = 1
  do while (i <= num_arg)
    call get_command_argument(number=i, value=arg, status=stat)
    if (i <= NUM_POSITIONAL_ARG) then
      select case (i)
      case (1)
        dataFilename = arg
        i = i + 1
      case (2)
        read(arg, *) numFrame ! in the unit of frame number
        i = i + 1
      case default
        if (myrank == root) then
          write(*,*) "Something is wrong in the codes; maybe NUM_POSITIONAL_ARG is not set correctly."
        end if
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end select
    else
      i = i + 1
      select case (arg)
      case ('-pa', '--pa')
        if (is_pm_mode) then
          if (myrank == root) then
            write(*,*) "-pa and -pm cannot be given at the same time!"
            call print_usage()
          end if
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if
        is_pa_mode = .true.
        call get_command_argument(i, topFilename)
        i = i + 1
        num_subArg = count_arg(i, num_arg)
        num_argPerMolType = 2
        if (mod(num_subArg, num_argPerMolType) > 0 .or. num_subArg < num_argPerMolType) then
          if (myrank == root) then
            write(*,*) "Wrong number of arguments for -pm: ", num_subArg + 1
            call print_usage()
          end if
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if

        numMolType = num_subArg / num_argPerMolType

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
          call get_command_argument(i, sys%mol(n)%type)
          i = i + 1
          call get_command_argument(i, arg)
          read(arg, *) start_index(n)
          i = i + 1
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

      case ('-pm', '--pm')
        if (is_pa_mode) then
          if (myrank == root) then
            write(*,*) "-pa and -pm cannot be given at the same time!"
            call print_usage()
          end if
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if
        is_pm_mode = .true.
        num_subArg = count_arg(i, num_arg)
        num_argPerMolType = 3
        if (mod(num_subArg, num_argPerMolType) > 0 .or. num_subArg < num_argPerMolType) then
          if (myrank == root) then
            write(*,*) "Wrong number of arguments for -pm: ", num_subArg
            call print_usage()
          end if
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if

        numMolType = num_subArg / num_argPerMolType

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

        do n = 1, numMolType
          call get_command_argument(i, sys%mol(n)%type) 
          i = i + 1
          call get_command_argument(i, arg) 
          i = i + 1
          read(arg, *) charge(n)
          call get_command_argument(i, arg) 
          i = i + 1
          read(arg, *) sys%mol(n)%num
        end do

        totNumMol = sum(sys%mol(:)%num)
        if (myrank == root) then
          write(*,*) "sys%mol%type = ", sys%mol%type
          write(*,*) "charge = ", charge
          write(*,*) "sys%mol%num = ", sys%mol%num
        end if

      case ('-o')
        call get_command_argument(i, outCorrFilename)
        i = i + 1

      case ('-s')
        call get_command_argument(i, arg) 
        i = i + 1
        read(arg, *) skip

      case ('-l')
        call get_command_argument(i, arg) ! in the unit of frame number
        i = i + 1
        read(arg, *) maxLag

      case ('-sd')
        is_sd = .true.

      case ('-ed')
        is_ed = .true.
        call get_command_argument(i, engtrjFilename)
        i = i + 1

      case ('-r')
        call get_command_argument(i, arg)
        i = i + 1
        read(arg, *) rBinWidth

      case ('-d')
        num_subArg = 2
        call get_command_argument(i, arg) 
        i = i + 1
        read(arg, *) numDomain_r

        call get_command_argument(i, arg) 
        i = i + 1
        read(arg, *) numDomain_c

      case default
        if (myrank == root) write(*,*) "Unknown argument: ", trim(adjustl(arg))
        call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        call exit(1)
      end select
    end if
  end do

  if (maxLag == -1) then
    maxLag = numFrame - 1
  end if

  !rank root output parameters read
  if (myrank == root) then
    write(*,*) "outFile = ", outCorrFilename
    write(*,*) "inFile.trr = ", dataFilename
    if (is_pa_mode) write(*,*) "topFile.top = ", topFilename
    write(*,*) "numFrame= ", numFrame
    write(*,*) "maxLag = ", maxLag 
    if (is_sd) write(*,*) "rBinWidth = ", rBinWidth
    write(*,*) "numMolType = ", numMolType
    write(*,*) "numDomain_r = ", numDomain_r
    write(*,*) "numDomain_c = ", numDomain_c
  end if

  if (myrank == root) then
    ! create an HDF5 file
    call H5open_f(ierr)
    call H5Fcreate_f(outCorrFilename, H5F_ACC_EXCL, outCorrFileid, ierr)
    if (ierr /= 0) then
      write(*,*) "Failed to create HDF5 file: ", outCorrFilename
      write(*,*) "Probably the file already exists?"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end if


  ! prepare eBinIndex for each rank
  if (is_ed) then
    call ed_readEng(engtrjFilename, numFrame, skip)
  end if

  call domainDecomposition(totNumMol, numFrame)

  !prepare memory for all ranks
  allocate(vel_r(3, numFrame, num_r), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel_r"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  allocate(vel_c(3, numFrame, num_c), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel_r"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  if (is_sd) then
    call sd_prepPosMemory(numFrame, totNumMol, num_r, num_c)
  end if

  !read trajectory at root
  if (myrank == root) then
    write(*,*) "start reading trajectory..."
    starttime = MPI_Wtime()
    sysNumAtom = get_natom(dataFilename)
    if (is_pm_mode .and. sysNumAtom /= totNumMol) then
      write(*,*) "sysNumAtom = ", sysNumAtom, ", totNumMol = ", totNumMol
      write(*,*) "In COM mode, sysNumAtom should equal to totNumMol!"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    write(*,*) "sysNumAtom=", sysNumAtom

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
      if (is_pm_mode) then
        vel(:, i, :) = vel_tmp
        if (is_sd) then
          pos(:, i, :) = pos_tmp
        end if
      else
        call com_vel(vel(:, i, :), vel_tmp, start_index, sys)
        if (is_sd) then
          call com_pos(pos(:, i, :), pos_tmp, start_index, sys, cell)
        end if
      end if
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
  else
    !not root, allocate dummy vel to inhibit error messages
    allocate(vel(1, 1, 1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: dummy vel on rank", myrank
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end if

  !distribute trajectory data collectively
  if (myrank == root) write(*,*) "start broadcasting trajectory"
  starttime = MPI_Wtime()
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

  deallocate(vel)

  call mpi_bcast(cell, 3, mpi_double_precision, root, MPI_COMM_WORLD, ierr)

  if (is_sd) call sd_broadcastPos()
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished broadcasting trajectory. It took ", endtime - starttime, " seconds"


  !decomposition
  if (is_sd) then
    if (myrank == root) write(*,*) "start spatial decomposition"
  else
    if (myrank == root) write(*,*) "start one-two decomposition"
  end if
  starttime = MPI_Wtime()

  allocate(vv(numFrame))
  if (stat /=0) then
    write(*,*) "Allocation failed: vv"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  allocate(nCorr(maxLag+1, numMolType*(numMolType+1)), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: nCorr"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  nCorr = 0d0

  if (is_sd) then
    call sd_cal_num_rBin(cell)
    call sd_prepCorrMemory(maxLag, numMolType, numFrame)
  else
    allocate(corr_tmp(2*maxLag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: corr_tmp"
      call exit(1)
    end if 
    corr_tmp = 0d0
  end if

  if (myrank == root) write(*,*) "time for allocation (sec):", MPI_Wtime() - starttime

  do j = c_start, c_end
    do i = r_start, r_end
      if (i == j) then
        if (myrank == root) write(*,*) "loop r =",i-r_start+1, " of ", num_r,&
                                          ", c =", j-c_start+1, " of ", num_c
        starttime2 = MPI_Wtime()
        molTypePairAllIndex = getMolTypeIndex(i, sys%mol(:)%num)
        if (is_sd) then
          do k = 1, maxLag+1
            numFrame_k = numFrame - k + 1
            vv(1:numFrame_k) = sum(vel_r(:, k:numFrame, i-r_start+1) * vel_c(:, 1:numFrame_k, j-c_start+1), 1)
            nCorr(k, molTypePairAllIndex) = nCorr(k, molTypePairAllIndex) + sum(vv(1:numFrame_k))
          end do
        else
          do k = 1, 3
            corr_tmp = corr(vel_r(k, :, i-r_start+1), maxLag)
            nCorr(:, molTypePairAllIndex) = nCorr(:, molTypePairAllIndex) + corr_tmp(maxLag+1:)
          end do
        end if
      else
        if (myrank == root) write(*,*) "loop r =",i-r_start+1, " of ", num_r,&
                                          ", c =", j-c_start+1, " of ", num_c
        starttime2 = MPI_Wtime()
        molTypePairIndex = getMolTypePairIndex(i, j, sys%mol(:)%num)
        molTypePairAllIndex = molTypePairIndex + numMolType
        if (is_sd) then
          call sd_getBinIndex(pos_r(:,:,i-r_start+1), pos_c(:,:,j-c_start+1), cell(1), sd_binIndex)
          do k = 1, maxLag+1
            numFrame_k = numFrame - k + 1
            vv(1:numFrame_k) = sum(vel_r(:, k:numFrame, i-r_start+1) * vel_c(:, 1:numFrame_k, j-c_start+1), 1)
            !TODO: test if this sum should be put here or inside the following loop for better performance
            nCorr(k, molTypePairAllIndex) = nCorr(k, molTypePairAllIndex) + sum(vv(1:numFrame_k))
            do n = 1, numFrame_k
              tmp_i = sd_binIndex(n)
              if (tmp_i <= num_rBin) then
                sdCorr(k, tmp_i, molTypePairIndex) = sdCorr(k, tmp_i, molTypePairIndex) + vv(n)
              end if
              !TODO: need test
              !nCorr(k, molTypePairIndex) = corr(k, molTypePairIndex) + vv(n)
            end do
          end do

          do t = 1, numFrame
            tmp_i = sd_binIndex(t)
            if (tmp_i <= num_rBin) then
              rho(tmp_i, molTypePairIndex) = rho(tmp_i, molTypePairIndex) + 1d0
            end if
          end do

        else ! one-two only
          do k = 1, 3
            corr_tmp = corr(vel_r(k, :, i-r_start+1), vel_c(k, :, j-c_start+1), maxLag)
            nCorr(:, molTypePairAllIndex) = nCorr(:, molTypePairAllIndex) + corr_tmp(maxLag+1:)
          end do
        end if
      end if
      if (myrank == root) write(*,*) "time for this loop (sec):", MPI_Wtime() - starttime2
      if (myrank == root) write(*,*) 
    end do
  end do
  if (is_sd) then
    deallocate(sd_binIndex)
  end if
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished decomposition. It took ", endtime - starttime, " seconds"

  !collect nCorr, sdCorr and rho
  if (myrank == root) write(*,*) "start collecting results"
  starttime = MPI_Wtime()
  if (myrank == root) then
    call mpi_reduce(MPI_IN_PLACE, nCorr, size(nCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  else
    call mpi_reduce(nCorr, dummy_null, size(nCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  end if
  if (is_sd) call sd_collectCorr()
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

    norm = [ (numFrame - (i-1), i = 1, maxLag+1) ] * 3d0
    do n = 1, numMolType*(numMolType+1)
      nCorr(:,n) = nCorr(:,n) / norm
    end do

    ! normalize rho and sdCorr
    if (is_sd) call sd_normalize(numFrame, numMolType, norm)

    deallocate(norm)

    !output results
    write(*,*) "start writing outputs"
    starttime = MPI_Wtime()
    call output()
    endtime = MPI_Wtime()
    write(*,*) "finished writing outputs. It took ", endtime - starttime, " seconds"
  end if

  if (myrank == root) write(*,*)
  if (myrank == root) write(*,*) "time for the whole program (sec):", MPI_Wtime() - prog_starttime
  call mpi_setup('stop')
  stop

contains
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
    integer(hid_t) :: did1, did2, did3, sid1, sid2

    allocate(timeLags(maxLag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: timeLags"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    timeLags = [ (dble(i), i = 0, maxLag) ] * timestep
    
    if (is_sd) then
      allocate(rBins(num_rBin), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: rBins"
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if 
      rBins = [ (i - 0.5d0, i = 1, num_rBin) ] * rBinWidth
    end if

    ! create and write dataset
    call H5LTset_attribute_double_f(outCorrFileid, "/", "timestep", [timestep], int(1, kind=size_t), ierr)
    call H5LTset_attribute_int_f(outCorrFileid, "/", "charge", charge, size(charge, kind=size_t), ierr)
    call H5LTset_attribute_int_f(outCorrFileid, "/", "numMol", sys%mol(:)%num, size(sys%mol(:)%num, kind=size_t), ierr)
    call H5LTset_attribute_double_f(outCorrFileid, "/", "cell", cell, size(cell, kind=size_t), ierr)

    call H5LTmake_dataset_double_f(outCorrFileid, "nCorr", 2, &
        [size(nCorr, 1, kind=hsize_t), size(nCorr, 2, kind=hsize_t)], nCorr, ierr)
    call H5Dopen_f(outCorrFileid, "nCorr", did1, ierr)

    if (is_sd) then
      call H5LTmake_dataset_double_f(outCorrFileid, "sdCorr", 3, &
          [size(sdCorr, 1, kind=hsize_t), size(sdCorr, 2, kind=hsize_t), size(sdCorr, 3, kind=hsize_t)], sdCorr, ierr)
      call H5Dopen_f(outCorrFileid, "sdCorr", did2, ierr)

      call H5LTmake_dataset_double_f(outCorrFileid, "rho", 2, &
          [size(rho, 1, kind=hsize_t), size(rho, 2, kind=hsize_t)], rho, ierr)
      call H5Dopen_f(outCorrFileid, "rho", did3, ierr)
    end if

    call H5LTmake_dataset_double_f(outCorrFileid, "timeLags", 1, [size(timeLags, kind=hsize_t)], timeLags, ierr)
    call H5Dopen_f(outCorrFileid, "timeLags", sid1, ierr)

    if (is_sd) then
      call H5LTmake_dataset_double_f(outCorrFileid, "rBins", 1, [size(rBins, kind=hsize_t)], rBins, ierr)
      call H5Dopen_f(outCorrFileid, "rBins", sid2, ierr)
    end if

    ! attach scale dimension
    call H5DSattach_scale_f(did1, sid1, 1, ierr)
    if (is_sd) then
      call H5DSattach_scale_f(did2, sid1, 1, ierr)
      call H5DSattach_scale_f(did2, sid2, 2, ierr)
      call H5DSattach_scale_f(did3, sid2, 1, ierr)
    end if

    call H5Dclose_f(sid1, ierr)
    call H5Dclose_f(did1, ierr)
    if (is_sd) then
      call H5Dclose_f(sid2, ierr)
      call H5Dclose_f(did2, ierr)
      call H5Dclose_f(did3, ierr)
    end if
    call H5Fclose_f(outCorrFileid, ierr)
    call H5close_f(ierr)
  end subroutine output

  integer function count_arg(i, num_arg)
    implicit none
    integer, intent(in) :: i, num_arg
    character(len=128) :: arg
    integer :: j, stat
    logical :: is_numeric
    !count number of arguments for the (i-1)-th option
    count_arg = 0
    j = i
    do while (.true.)
      if (j > num_arg) then
        return
      end if
      call get_command_argument(number=j, value=arg, status=stat)
      if (stat /= 0) then
        if (myrank == root) then
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          write(*,*) "Error: unable to count the number of arguments for the ", i-1, "-th option"
          call print_usage()
          call exit(1)
        end if
      else if (arg(1:1) == '-' ) then
        is_numeric = verify(arg(2:2), '0123456789') .eq. 0 
        if (.not. is_numeric) return !end of data file arguments
      end if
      j = j + 1
      count_arg = count_arg + 1
    end do
  end function count_arg

  subroutine print_usage()
    implicit none
    write(*, *) "usage: $ decompose_mpi <infile.trr> <numFrameToRead> <-pa | -pm ...> [options]"
    write(*, *) "options: "
    write(*, *) "  -pa <topfile.top> <molecule1> <start_index1> [<molecule2> <start_index2>...]:"
    write(*, *) "   read parameters from topology file. ignored when -pm is given"
    write(*, *) 
    write(*, *) "  -pm <molecule1> <charge1> <number1> [<molecule2> <charge2> <number2>...]:"
    write(*, *) "   manually assign parameters for single-atom-molecule system"
    write(*, *) 
    write(*, *) "  -o <outfile>: output filename. default = corr.h5"
    write(*, *) 
    write(*, *) "  -s <skip>: skip=1 means no frames are skipped, which is default."
    write(*, *) "             skip=2 means reading every 2nd frame."
    write(*, *) 
    write(*, *) "  -l <maxlag>: maximum time lag. default = <numFrameToRead - 1>"
    write(*, *) 
    write(*, *) "  -sd: do spatial decomposition. default no sd."
    write(*, *)
    write(*, *) "  -ed <engtrj.dat>: do energy decomposition. default no ed."
    write(*, *)
    write(*, *) "  -sbwidth <sBinWidth(nm)>: spatial-decomposition bin width. default = 0.01."
    write(*, *) "                            only meaningful when -sd is given."
    write(*, *)
    write(*, *) "  -ebnum <num_eBin>: number of energy-decomposition bins. default = 500"
    write(*, *) "                     only meaningful when -ed is given."
    write(*, *) 
    write(*, *) "  -d <numDomain_r> <numDomain_c>:" 
    write(*, *) "   manually assign the MPI decomposition pattern"
  end subroutine print_usage
end program decompose_mpi
