program decond
  use mpiproc
  use HDF5
  use utility, only: handle, getMolTypePairIndexFromTypes, newunit
  use xdr, only: open_trajectory, close_trajectory, read_trajectory, get_natom
  use top, only: open_top, close_top, read_top, system, print_sys
  use correlation
  use spatial_dec, only: sd_init, rBinWidth, pos_r, pos_c, sdPairCount, sdCorr, pos, num_rBin, &
                         com_pos, sd_binIndex, sd_prepCorrMemory, sd_getBinIndex, &
                         sd_cal_num_rBin, sd_broadcastPos, sd_prepPosMemory, &
                         sd_collectCorr, sd_average, sd_make_rBins, sd_finish
  use energy_dec, only: engtrajFilename, ed_readEng, ed_getBinIndex, ed_binIndex, num_eBin, &
                        ed_binIndex, ed_prepCorrMemory, ed_getBinIndex, skipEng, &
                        ed_collectCorr, ed_average, edPairCount, edCorr, eBinWidth, &
                        ed_init, ed_make_eBins, ed_finish


  implicit none
  character(len=*), parameter :: DECOND_VERSION = "0.4.0"
  integer, parameter :: NUM_POSITIONAL_ARG = 3, LEAST_REQUIRED_NUM_ARG = 7
  integer :: num_arg, num_subArg, num_argPerMolType
  integer :: i, j, k, n, totNumMol, t, sysNumAtom
  character(len=128) :: outCorrFilename, dataFilename, logFilename, topFilename, arg
  type(handle) :: dataFileHandle, topFileHandle
  integer :: numFrame, maxLag, stat, numMolType, numFrameRead, numFrame_k
  integer :: molTypePairIndex, molTypePairAllIndex, tmp_i, skipTrr, numMolTypePairAll
  integer, allocatable :: charge(:), frameCount(:), start_index(:)
  real(8) :: cell(3), timestep, tmp_r, temperature
  real(8), allocatable :: pos_tmp(:, :), vel_tmp(:, :), vv(:)
  !one frame data (dim=3, atom) 
  real(8), allocatable :: vel(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(8), allocatable :: time(:), nCorr(:, :), corr_tmp(:)
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
  outCorrFilename = 'corr.c5'
  skipTrr = 1
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
        logFilename = arg
        i = i + 1
        temperature = getTfromLog(logFilename)
      case (3)
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
          charge(n) = nint(sum(sys%mol(n)%atom(:)%charge))
        end do
        totNumMol = sum(sys%mol(:)%num)
        if (myrank == root) write(*,*) "charge = ", charge

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

      case ('-skiptrr')
        call get_command_argument(i, arg) 
        i = i + 1
        read(arg, *) skipTrr

      case ('-skipeng')
        call get_command_argument(i, arg)
        i = i + 1
        read(arg, *) skipEng


      case ('-l')
        call get_command_argument(i, arg) ! in the unit of frame number
        i = i + 1
        read(arg, *) maxLag

      case ('-sd')
        is_sd = .true.

      case ('-ed')
        is_ed = .true.
        call get_command_argument(i, engtrajFilename)
        i = i + 1

      case ('-sbwidth')
        call get_command_argument(i, arg)
        i = i + 1
        read(arg, *) rBinWidth

      case ('-ebwidth')
        call get_command_argument(i, arg)
        i = i + 1
        read(arg, *) eBinWidth

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

  !auto1, auto2, ...autoN, cross11, cross12, ..., cross1N, cross22, ...cross2N, cross33,..., crossNN
  != [auto part] + [cross part]
  != [numMolType] + [numMolType * (numMolType + 1) / 2]
  numMolTypePairAll = numMolType * (numMolType + 3) / 2

  if (maxLag == -1) then
    maxLag = numFrame - 1
  end if

  !rank root output parameters read
  if (myrank == root) then
    write(*,*) "outFile = ", outCorrFilename
    write(*,*) "inFile.trr = ", dataFilename
    write(*,*) "logFile = ", logFilename
    write(*,*) "temperature = ", temperature
    if (is_pa_mode) write(*,*) "topFile.top = ", topFilename
    write(*,*) "numFrame= ", numFrame
    write(*,*) "maxLag = ", maxLag 
    if (is_sd) write(*,*) "rBinWidth = ", rBinWidth
    if (is_ed) write(*,*) "eBinWidth = ", eBinWidth
    write(*,*) "numMolType = ", numMolType
    write(*,*) "numDomain_r = ", numDomain_r
    write(*,*) "numDomain_c = ", numDomain_c
  end if

  if (myrank == root) then
    ! create an HDF5 file
    call H5open_f(ierr)
    call H5Fcreate_f(outCorrFilename, H5F_ACC_EXCL_F, outCorrFileid, ierr)
    if (ierr /= 0) then
      write(*,*) "Failed to create HDF5 file: ", trim(adjustl(outCorrFilename))
      write(*,*) "Probably the file already exists?"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end if

  call domainDecomposition(totNumMol, numFrame) ! determine r_start, c_start ...etc.

  ! prepare eBinIndex for each rank
  if (is_ed) then
    if (myrank == root) then
      write(*,*) "start preparing eBinIndex..."
      starttime = MPI_Wtime()
    end if
    call ed_readEng(engtrajFilename, numFrame)
    if (myrank == root) then
      endtime = MPI_Wtime()
      write(*,*) "finished preparing eBinIndex. It took ", endtime - starttime, "seconds"
    end if
  end if

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
      do j = 1, skipTrr-1
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
  if (myrank == root) then
    write(*,*) "start one-two decomposition"
    if (is_sd) write(*,*) "start spatial decomposition"
    if (is_ed) write(*,*) "start energy decomposition"
  end if
  starttime = MPI_Wtime()

  allocate(vv(numFrame))
  if (stat /=0) then
    write(*,*) "Allocation failed: vv"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 

  allocate(nCorr(maxLag+1, numMolTypePairAll), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: nCorr"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if 
  nCorr = 0d0

  if (is_sd .or. is_ed) then
    if (is_sd) then
      call sd_cal_num_rBin(cell)
      call sd_prepCorrMemory(maxLag, numMolType, numFrame)
    end if
    if (is_ed) call ed_prepCorrMemory(maxLag, numMolType, numFrame)
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
        molTypePairAllIndex = getMolTypeIndex(i, sys%mol(:)%num, numMolType)
        if (is_sd .or. is_ed) then
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
        molTypePairIndex = getMolTypePairIndex(i, j, sys%mol(:)%num, numMolType)
        molTypePairAllIndex = molTypePairIndex + numMolType
        if (is_sd .or. is_ed) then
          if (is_sd) call sd_getBinIndex(i-r_start+1, j-c_start+1, cell, sd_binIndex)
          if (is_ed) call ed_getBinIndex(i, j, ed_binIndex)
          do k = 1, maxLag+1
            numFrame_k = numFrame - k + 1
            vv(1:numFrame_k) = sum(vel_r(:, k:numFrame, i-r_start+1) * vel_c(:, 1:numFrame_k, j-c_start+1), 1)
            !TODO: test if this sum should be put here or inside the following loop for better performance
            nCorr(k, molTypePairAllIndex) = nCorr(k, molTypePairAllIndex) + sum(vv(1:numFrame_k))
            do n = 1, numFrame_k
              if (is_sd) then
                tmp_i = sd_binIndex(n)
                if (tmp_i <= num_rBin) then
                  sdCorr(k, tmp_i, molTypePairIndex) = sdCorr(k, tmp_i, molTypePairIndex) + vv(n)
                end if
              end if
              if (is_ed) then
                tmp_i = ed_binIndex(n)
                if (tmp_i <= num_rBin) then
                  edCorr(k, tmp_i, molTypePairIndex) = edCorr(k, tmp_i, molTypePairIndex) + vv(n)
                end if
              end if
              !TODO: need test
              !nCorr(k, molTypePairIndex) = corr(k, molTypePairIndex) + vv(n)
            end do
          end do

          do t = 1, numFrame
            if (is_sd) then
              tmp_i = sd_binIndex(t)
              if (tmp_i <= num_rBin) then
                sdPairCount(tmp_i, molTypePairIndex) = sdPairCount(tmp_i, molTypePairIndex) + 1d0
              end if
            end if
            if (is_ed) then
              tmp_i = ed_binIndex(t)
              if (tmp_i <= num_eBin) then
                edPairCount(tmp_i, molTypePairIndex) = edPairCount(tmp_i, molTypePairIndex) + 1d0
              end if
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
  if (is_sd) deallocate(sd_binIndex)
  if (is_ed) deallocate(ed_binIndex)

  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished decomposition. It took ", endtime - starttime, " seconds"

  !collect nCorr
  if (myrank == root) write(*,*) "start collecting results"
  starttime = MPI_Wtime()
  if (myrank == root) then
    call mpi_reduce(MPI_IN_PLACE, nCorr, size(nCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  else
    call mpi_reduce(nCorr, dummy_null, size(nCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  end if
  if (is_sd) call sd_collectCorr()
  if (is_ed) call ed_collectCorr()
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished collecting results. It took ", endtime - starttime, " seconds"

  !average at root and then output
  if (myrank == root) then
    allocate(frameCount(maxLag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: frameCount"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 

    do j = 1, numMolType
      do i = j, numMolType
        if (i /= j) then
          molTypePairIndex = getMolTypePairIndexFromTypes(i, j, numMolType)
          molTypePairAllIndex = molTypePairIndex + numMolType
          nCorr(:, molTypePairAllIndex) = nCorr(:, molTypePairAllIndex) / 2d0
        end if
      end do
    end do

    frameCount = [ (numFrame - (i-1), i = 1, maxLag+1) ] * 3d0
    do n = 1, numMolTypePairAll
      nCorr(:,n) = nCorr(:,n) / frameCount
    end do

    if (is_sd) call sd_average(numFrame, numMolType, frameCount)
    if (is_ed) call ed_average(numFrame, numMolType, frameCount)

    deallocate(frameCount)

    !output results
    write(*,*) "start writing outputs"
    starttime = MPI_Wtime()
    call output()
    endtime = MPI_Wtime()
    write(*,*) "finished writing outputs. It took ", endtime - starttime, " seconds"
  end if

  if (myrank == root) write(*,*)
  if (myrank == root) write(*,*) "time for the whole program (sec):", MPI_Wtime() - prog_starttime

  deallocate(nCorr, charge, vel_r, vel_c, start_index)
  if (is_sd) call sd_finish()
  if (is_ed) call ed_finish()
  call mpi_setup('stop')
  stop

contains
  integer function getMolTypeIndex(i, numMol, numMolType)
    implicit none
    integer, intent(in) :: i, numMol(:), numMolType
    integer :: n, numMol_acc

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

  integer function getMolTypePairIndex(i, j, numMol, numMolType)
    implicit none
    integer, intent(in) :: i, j, numMol(:), numMolType
    integer :: r, c, ii, jj
    !          c
    !    | 1  2  3  4
    !  --+------------
    !  1 | 1  2  3  4
    !    |
    !  2 |    5  6  7
    !r   |
    !  3 |       8  9
    !    |
    !  4 |         10
    !
    !  index(r, c) = (r - 1) * n + c - r * (r - 1) / 2
    !  where n = size(c) = size(r), r <= c


    ii = getMolTypeIndex(i, numMol, numMolType)
    jj = getMolTypeIndex(j, numMol, numMolType)
    r = min(ii, jj)
    c = max(ii, jj)
    getMolTypePairIndex = (r - 1) * numMolType + c - r * (r - 1) / 2
  end function getMolTypePairIndex

  subroutine com_vel(com_v, vel, start_index, sys)
    implicit none
    real(8), dimension(:, :), intent(out) :: com_v
    real(8), dimension(:, :), intent(in) :: vel 
    integer, dimension(:), intent(in) :: start_index
    type(system), intent(in) :: sys
    integer :: d, i, j, idx_begin, idx_end, idx_com

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
    use spatial_dec, only: rBins
    use energy_dec, only: eBins, engMin_global
    implicit none
    real(8), allocatable :: timeLags(:)
    integer :: ierr
    integer(hid_t) :: dset_nCorr, dset_timeLags, &
                      dset_sdCorr, dset_sdPairCount, dset_rBins, &
                      dset_edCorr, dset_edPairCount, dset_eBins, &
                      grp_sd_id, grp_ed_id, space_id, dset_id

    !HDF5:
    character(len=*), parameter :: GROUP_ROOT = "/", &
                                   GROUP_SPATIAL = "spatialDec", &
                                   GROUP_ENERGY = "energyDec"
    !/Attributes
    character(len=*), parameter :: ATTR_VERSION = "version", &
                                   ATTR_TYPE = "type", &
                                   ATTR_UNIT = "unit", &
                                   OUT_TYPE = "CorrFile"
    !/Dataset
    character(len=*), parameter :: DSETNAME_VOLUME = "volume", &
                                   DSETNAME_TEMP = "temperature", &
                                   DSETNAME_CHARGE = "charge", &
                                   DSETNAME_NUMMOL = "numMol", &
                                   DSETNAME_TIMELAGS = "timeLags", &
                                   DSETNAME_NCORR = "nCorr"

    !/GROUP_SPATIAL or GROUP_ENERGY/Dataset
    character(len=*), parameter :: DSETNAME_DECBINS = "decBins", &
                                   DSETNAME_DECCORR = "decCorr", &
                                   DSETNAME_DECPAIRCOUNT = "decPairCount"


    allocate(timeLags(maxLag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: timeLags"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    timeLags = [ (dble(i), i = 0, maxLag) ] * timestep
    
    if (is_sd) then
      call sd_make_rBins()
    end if

    if (is_ed) then
      call ed_make_eBins()
    end if

    !create and write attributes
    call H5LTset_attribute_string_f(outCorrFileid, GROUP_ROOT, ATTR_VERSION, DECOND_VERSION, ierr)
    call H5LTset_attribute_string_f(outCorrFileid, GROUP_ROOT, ATTR_TYPE, OUT_TYPE, ierr)

    !create and write datasets
    !volume
    call H5Screate_f(H5S_SCALAR_F, space_id, ierr)
    call H5Dcreate_f(outCorrFileid, DSETNAME_VOLUME, H5T_NATIVE_DOUBLE, space_id, dset_id, ierr)
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, product(cell), [0_hsize_t], ierr)
    call H5Dclose_f(dset_id, ierr)
    call H5Sclose_f(space_id, ierr)
    call H5LTset_attribute_string_f(outCorrFileid, DSETNAME_VOLUME, ATTR_UNIT, "nm$^3$", ierr)

    !temperature
    call H5Screate_f(H5S_SCALAR_F, space_id, ierr)
    call H5Dcreate_f(outCorrFileid, DSETNAME_TEMP, H5T_NATIVE_DOUBLE, space_id, dset_id, ierr)
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, temperature, [0_hsize_t], ierr)
    call H5Dclose_f(dset_id, ierr)
    call H5Sclose_f(space_id, ierr)
    call H5LTset_attribute_string_f(outCorrFileid, DSETNAME_TEMP, ATTR_UNIT, "K", ierr)

    !charge
    call H5LTmake_dataset_int_f(outCorrFileid, DSETNAME_CHARGE, 1, [size(charge, kind=hsize_t)], charge, ierr)
    call H5LTset_attribute_string_f(outCorrFileid, DSETNAME_CHARGE, ATTR_UNIT, "e", ierr)

    !numMol
    call H5LTmake_dataset_int_f(outCorrFileid, DSETNAME_NUMMOL, 1, &
                                   [size(sys%mol(:)%num, kind=size_t)], sys%mol(:)%num, ierr)

    !timeLags
    call H5LTmake_dataset_double_f(outCorrFileid, DSETNAME_TIMELAGS, 1, [size(timeLags, kind=hsize_t)], timeLags, ierr)
    call H5Dopen_f(outCorrFileid, DSETNAME_TIMELAGS, dset_timeLags, ierr)
    call H5LTset_attribute_string_f(outCorrFileid, DSETNAME_TIMELAGS, ATTR_UNIT, "ps", ierr)

    !nCorr
    call H5LTmake_dataset_double_f(outCorrFileid, DSETNAME_NCORR, 2, &
        [size(nCorr, 1, kind=hsize_t), size(nCorr, 2, kind=hsize_t)], nCorr, ierr)
    call H5Dopen_f(outCorrFileid, DSETNAME_NCORR, dset_nCorr, ierr)
    call H5LTset_attribute_string_f(outCorrFileid, DSETNAME_NCORR, ATTR_UNIT, "nm$^2$ ps$^{-2}$", ierr)

    if (is_sd) then
      !create a group for storing spatial-decomposition data
      call H5Gcreate_f(outCorrFileid, GROUP_SPATIAL, grp_sd_id, ierr)

      !decCorr
      call H5LTmake_dataset_double_f(grp_sd_id, DSETNAME_DECCORR, 3, &
          [size(sdCorr, 1, kind=hsize_t), size(sdCorr, 2, kind=hsize_t), size(sdCorr, 3, kind=hsize_t)], sdCorr, ierr)
      call H5Dopen_f(grp_sd_id, DSETNAME_DECCORR, dset_sdCorr, ierr)
      call H5LTset_attribute_string_f(grp_sd_id, DSETNAME_DECCORR, ATTR_UNIT, "nm$^2$ ps$^{-2}$", ierr)

      !decPairCount
      call H5LTmake_dataset_double_f(grp_sd_id, DSETNAME_DECPAIRCOUNT, 2, &
          [size(sdPairCount, 1, kind=hsize_t), size(sdPairCount, 2, kind=hsize_t)], sdPairCount, ierr)
      call H5Dopen_f(grp_sd_id, DSETNAME_DECPAIRCOUNT, dset_sdPairCount, ierr)

      !decBins
      call H5LTmake_dataset_double_f(grp_sd_id, DSETNAME_DECBINS, 1, [size(rBins, kind=hsize_t)], rBins, ierr)
      call H5Dopen_f(grp_sd_id, DSETNAME_DECBINS, dset_rBins, ierr)
      call H5LTset_attribute_string_f(grp_sd_id, DSETNAME_DECBINS, ATTR_UNIT, "nm", ierr)
    end if

    if (is_ed) then
      !create a group for storing energy-decomposition data
      call H5Gcreate_f(outCorrFileid, GROUP_ENERGY, grp_ed_id, ierr)

      !decCorr
      call H5LTmake_dataset_double_f(grp_ed_id, DSETNAME_DECCORR, 3, &
          [size(edCorr, 1, kind=hsize_t), size(edCorr, 2, kind=hsize_t), size(edCorr, 3, kind=hsize_t)], edCorr, ierr)
      call H5Dopen_f(grp_ed_id, DSETNAME_DECCORR, dset_edCorr, ierr)
      call H5LTset_attribute_string_f(grp_ed_id, DSETNAME_DECCORR, ATTR_UNIT, "nm$^2$ ps$^{-2}$", ierr)

      !decPairCount
      call H5LTmake_dataset_double_f(grp_ed_id, DSETNAME_DECPAIRCOUNT, 2, &
          [size(edPairCount, 1, kind=hsize_t), size(edPairCount, 2, kind=hsize_t)], edPairCount, ierr)
      call H5Dopen_f(grp_ed_id, DSETNAME_DECPAIRCOUNT, dset_edPairCount, ierr)

      !decBins
      call H5LTmake_dataset_double_f(grp_ed_id, DSETNAME_DECBINS, 1, [size(eBins, kind=hsize_t)], eBins, ierr)
      call H5Dopen_f(grp_ed_id, DSETNAME_DECBINS, dset_eBins, ierr)
      call H5LTset_attribute_string_f(grp_ed_id, DSETNAME_DECBINS, ATTR_UNIT, "kcal mol$^{-1}$", ierr)
    end if

    !attach scale dimension
    !dimension index is ordered in row major as (dimN, dimN-1, ..., dim2, dim1)
    call H5DSattach_scale_f(dset_nCorr, dset_timeLags, 2, ierr)
    if (is_sd) then
      call H5DSattach_scale_f(dset_sdCorr, dset_timeLags, 3, ierr)
      call H5DSattach_scale_f(dset_sdCorr, dset_rBins, 2, ierr)
      call H5DSattach_scale_f(dset_sdPairCount, dset_rBins, 2, ierr)
    end if
    if (is_ed) then
      call H5DSattach_scale_f(dset_edCorr, dset_timeLags, 3, ierr)
      call H5DSattach_scale_f(dset_edCorr, dset_eBins, 2, ierr)
      call H5DSattach_scale_f(dset_edPairCount, dset_eBins, 2, ierr)
    end if

    call H5Dclose_f(dset_timeLags, ierr)
    call H5Dclose_f(dset_nCorr, ierr)
    if (is_sd) then
      call H5Dclose_f(dset_rBins, ierr)
      call H5Dclose_f(dset_sdCorr, ierr)
      call H5Dclose_f(dset_sdPairCount, ierr)
    end if
    if (is_ed) then
      call H5Dclose_f(dset_eBins, ierr)
      call H5Dclose_f(dset_edCorr, ierr)
      call H5Dclose_f(dset_edPairCount, ierr)
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

  real(8) function getTfromLog(logFilename)
    implicit none
    character(len=*), intent(in) :: logFilename
    integer :: logio, idx, temp_idx, stat
    character(len=*), parameter :: LINE_LEN_STR = '128'
    integer, parameter :: RECORD_LEN = 15
    character(len=128) :: line
    logical :: found_average

    found_average = .false.
    open(unit=newunit(logio), file=logFilename, status='old', action="READ", form="FORMATTED")
    do while(.true.)
      read(logio, "(A"//LINE_LEN_STR//")", iostat=stat) line
      if (stat > 0) then
        write(*,*) "Error reading line"
        call exit(1)
      else if (stat < 0) then
        write(*,*) "Unable to find 'A V E R A G E S' in logfile ", trim(adjustl(logFilename))
        call exit(1)
      end if

      if (found_average) then
        idx = index(line, "Temperature")
        if (idx > 0) then
          temp_idx = ceiling(real(idx, 8) / RECORD_LEN)
          read(logio, "(A"//LINE_LEN_STR//")", iostat=stat) line
          if (stat > 0) then
            write(*,*) "Error reading temperature record"
            call exit(1)
          end if
          read(line(RECORD_LEN * (temp_idx - 1) + 1: RECORD_LEN * temp_idx), *) getTfromLog
          return
        end if
      else
        idx = index(line, "A V E R A G E S")
        if (idx > 0) then
          found_average = .true.
        end if
      end if
    end do
  end function getTfromLog

  subroutine print_usage()
    implicit none
    write(*, *) "usage: $ decond <trrfile> <logfile> <numFrameToRead> <-pa | -pm ...> [options]"
    write(*, *) "options: "
    write(*, *) "  -pa <topfile.top> <molecule1> <start_index1> [<molecule2> <start_index2>...]:"
    write(*, *) "   read parameters from topology file. ignored when -pm is given"
    write(*, *) 
    write(*, *) "  -pm <molecule1> <charge1> <number1> [<molecule2> <charge2> <number2>...]:"
    write(*, *) "   manually assign parameters for single-atom-molecule system"
    write(*, *) 
    write(*, *) "  -o <outfile>: output filename. default = corr.h5"
    write(*, *) 
    write(*, *) "  -skiptrr <skip>: skip=1 means no frames are skipped, which is default."
    write(*, *) "             skip=2 means reading every 2nd frame."
    write(*, *)
    write(*, *) "  -skipeng <skip>: skip=1 means no frames are skipped, which is default."
    write(*, *) "             skip=2 means reading every 2nd frame."
    write(*, *) 
    write(*, *) "  -l <maxlag>: maximum time lag in frames. default = <numFrameToRead - 1>"
    write(*, *) 
    write(*, *) "  -sd: do spatial decomposition. default no sd."
    write(*, *)
    write(*, *) "  -ed <engtraj.dat>: do energy decomposition. default no ed."
    write(*, *)
    write(*, *) "  -sbwidth <sBinWidth(nm)>: spatial-decomposition bin width. default = 0.01."
    write(*, *) "                            only meaningful when -sd is given."
    write(*, *)
    write(*, *) "  -ebwidth <eBinWidth(kcal/mol)>: energy-decomposition bin width. default = 0.1"
    write(*, *) "                                  only meaningful when -ed is given."
    write(*, *) 
    write(*, *) "  -d <numDomain_r> <numDomain_c>:" 
    write(*, *) "   manually assign the MPI decomposition pattern"
  end subroutine print_usage
end program decond
