program decond
  use mpiproc
  use hdf5
  use utility, only: handle, get_pairindex_upper_diag, newunit
  use xdr, only: open_trajectory, close_trajectory, read_trajectory, get_natom
  use top, only: open_top, close_top, read_top, system, print_sys
  use correlation, only: corr
  use spatial_dec, only: sd_init, rBinWidth, pos_r, pos_c, sdPairCount, &
                         sdCorr, pos, num_rBin, com_pos, sd_binIndex, &
                         sd_prepCorrMemory, sd_getBinIndex, &
                         sd_cal_num_rBin, sd_broadcastPos, sd_prepPosMemory, &
                         sd_collectCorr, sd_average, sd_make_rBins, sd_finish
  use energy_dec, only: engfiles, ed_readEng, ed_getBinIndex, &
                        num_eBin, ed_binIndex, ed_prepCorrMemory, &
                        skipEng, ed_collectCorr, ed_average, &
                        edPairCount, edCorr, eBinWidth, ed_init, &
                        ed_make_eBins, ed_finish, num_engfiles, &
                        ed_prep_engfiles
  use varpars, only: line_len, dp, decond_version, num_positional_arg, &
                     least_required_num_arg, num_arg, num_subarg, &
                     num_arg_per_moltype, totnummol, sysnumatom, &
                     corrfile, datafile, logfile, topfile, arg, &
                     line_len_str, numframe, maxlag, nummoltype, &
                     skiptrr, num_moltypepair_all, charge, framecount, &
                     start_index, cell, timestep, temperature, &
                     is_periodic, is_pa_mode, is_pm_mode, is_sd, is_ed
  implicit none
  integer :: i, j, k, n, t, stat, numframe_k, numframe_read, tmp_i
  integer :: moltypepair_idx, moltypepair_allidx
  real(dp) :: tmp_r
  type(handle) :: datafileio, topfileio
  real(dp), allocatable :: pos_tmp(:, :), vel_tmp(:, :), vv(:)
  !one frame data (dim=3, atom)
  real(dp), allocatable :: vel(:, :, :)
  !pos(dim=3, timeFrame, atom), vel(dim=3, timeFrame, atom)
  real(dp), allocatable :: time(:), ncorr(:, :), corr_tmp(:)
  type(system) :: sys
  !MPI variables
  real(dp), allocatable :: vel_r(:, :, :), vel_c(:, :, :)
  integer(hid_t) :: corrfileio

  !initialize
  call mpi_setup('init')

  prog_starttime = MPI_Wtime()

  num_arg = command_argument_count()
  is_pa_mode = .false.
  is_pm_mode = .false.

  !default values
  corrfile = 'corr.c5'
  skiptrr = 1
  maxlag = -1
  is_sd = .false.
  is_ed = .false.
  numDomain_r = 0
  numDomain_c = 0
  call sd_init()
  call ed_init()

  !root checks the number of the input arguments
  is_periodic = .true.
  if (num_arg < least_required_num_arg) then
    if (myrank == root) call print_usage()
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if

  !read parameters for all ranks
  i = 1
  do while (i <= num_arg)
    call get_command_argument(number=i, value=arg, status=stat)
    if (i <= num_positional_arg) then
      select case (i)
      case (1)
        datafile = arg
        i = i + 1
      case (2)
        logfile = arg
        i = i + 1
        temperature = getTfromLog(logfile)
      case (3)
        read(arg, *) numframe ! in the unit of frame number
        i = i + 1
      case default
        if (myrank == root) then
          write(*,*) "Something is wrong in the codes; maybe num_positional_arg is not set correctly."
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
        call get_command_argument(i, topfile)
        i = i + 1
        num_subarg = count_arg(i, num_arg)
        num_arg_per_moltype = 2
        if (mod(num_subarg, num_arg_per_moltype) > 0 .or. num_subarg < num_arg_per_moltype) then
          if (myrank == root) then
            write(*,*) "Wrong number of arguments for -pm: ", num_subarg + 1
            call print_usage()
          end if
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if

        nummoltype = num_subarg / num_arg_per_moltype

        allocate(sys%mol(nummoltype), stat=stat)
        if (stat /=0) then
          write(*,*) "Allocation failed: sys%mol"
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if

        allocate(charge(nummoltype), stat=stat)
        if (stat /=0) then
          write(*,*) "Allocation failed: charge"
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if

        allocate(start_index(nummoltype), stat=stat)
        if (stat /=0) then
          write(*,*) "Allocation failed: start_index"
          call exit(1)
        end if

        do n = 1, nummoltype
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
        topfileio = open_top(topfile)
        call read_top(topfileio, sys)
        call close_top(topfileio)
        if (myrank == root) call print_sys(sys)

        do n = 1, nummoltype
          charge(n) = nint(sum(sys%mol(n)%atom(:)%charge))
        end do
        totnummol = sum(sys%mol(:)%num)
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
        num_subarg = count_arg(i, num_arg)
        num_arg_per_moltype = 3
        if (mod(num_subarg, num_arg_per_moltype) > 0 .or. num_subarg < num_arg_per_moltype) then
          if (myrank == root) then
            write(*,*) "Wrong number of arguments for -pm: ", num_subarg
            call print_usage()
          end if
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if

        nummoltype = num_subarg / num_arg_per_moltype

        allocate(sys%mol(nummoltype), stat=stat)
        if (stat /=0) then
          write(*,*) "Allocation failed: sys%mol"
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if

        allocate(charge(nummoltype), stat=stat)
        if (stat /=0) then
          write(*,*) "Allocation failed: charge"
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        end if

        do n = 1, nummoltype
          call get_command_argument(i, sys%mol(n)%type)
          i = i + 1
          call get_command_argument(i, arg)
          i = i + 1
          read(arg, *) charge(n)
          call get_command_argument(i, arg)
          i = i + 1
          read(arg, *) sys%mol(n)%num
        end do

        totnummol = sum(sys%mol(:)%num)
        if (myrank == root) then
          write(*,*) "sys%mol%type = ", sys%mol%type
          write(*,*) "charge = ", charge
          write(*,*) "sys%mol%num = ", sys%mol%num
        end if

      case ('-o')
        call get_command_argument(i, corrfile)
        i = i + 1

      case ('-skiptrr')
        call get_command_argument(i, arg)
        i = i + 1
        read(arg, *) skiptrr

      case ('-skipeng')
        call get_command_argument(i, arg)
        i = i + 1
        read(arg, *) skipEng


      case ('-l')
        call get_command_argument(i, arg) ! in the unit of frame number
        i = i + 1
        read(arg, *) maxlag

      case ('-sd')
        is_sd = .true.

      case ('-ed')
        is_ed = .true.
        num_engfiles = count_arg(i, num_arg)
        call ed_prep_engfiles()
        do n = 1, num_engfiles
          call get_command_argument(i, engfiles(n))
          i = i + 1
        end do

      case ('-sbwidth')
        call get_command_argument(i, arg)
        i = i + 1
        read(arg, *) rBinWidth

      case ('-ebwidth')
        call get_command_argument(i, arg)
        i = i + 1
        read(arg, *) eBinWidth

      case ('-d')
        num_subarg = 2
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
  != [nummoltype] + [nummoltype * (nummoltype + 1) / 2]
  num_moltypepair_all = nummoltype * (nummoltype + 3) / 2

  if (maxlag == -1) then
    maxlag = numframe - 1
  end if

  !rank root output parameters read
  if (myrank == root) then
    write(*,*) "outFile = ", corrfile
    write(*,*) "inFile.trr = ", datafile
    write(*,*) "logFile = ", logfile
    write(*,*) "temperature = ", temperature
    if (is_pa_mode) write(*,*) "topFile.top = ", topfile
    write(*,*) "numframe= ", numframe
    write(*,*) "maxlag = ", maxlag
    if (is_sd) write(*,*) "rBinWidth = ", rBinWidth
    if (is_ed) write(*,*) "eBinWidth = ", eBinWidth
    write(*,*) "nummoltype = ", nummoltype
    write(*,*) "numDomain_r = ", numDomain_r
    write(*,*) "numDomain_c = ", numDomain_c
  end if

  if (myrank == root) then
    ! create an HDF5 file
    call H5open_f(ierr)
    call H5Fcreate_f(corrfile, H5F_ACC_EXCL_F, corrfileio, ierr)
    if (ierr /= 0) then
      write(*,*) "Failed to create HDF5 file: ", trim(adjustl(corrfile))
      write(*,*) "Probably the file already exists?"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end if

  call domain_dec(totnummol, numframe) ! determine r_start, c_start ...etc.

  ! prepare eBinIndex for each rank
  if (is_ed) then
    if (myrank == root) then
      write(*,*) "start preparing eBinIndex..."
      starttime = MPI_Wtime()
    end if
    call ed_readEng(numframe, totnummol)
    if (myrank == root) then
      endtime = MPI_Wtime()
      write(*,*) "finished preparing eBinIndex. It took ", endtime - starttime, "seconds"
    end if
  end if

  !prepare memory for all ranks
  allocate(vel_r(3, numframe, num_r), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel_r"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if
  allocate(vel_c(3, numframe, num_c), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: vel_r"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if

  if (is_sd) then
    call sd_prepPosMemory(numframe, totnummol, num_r, num_c)
  end if

  !read trajectory at root
  if (myrank == root) then
    write(*,*) "start reading trajectory..."
    starttime = MPI_Wtime()
    sysnumatom = get_natom(datafile)
    if (is_pm_mode .and. sysnumatom /= totnummol) then
      write(*,*) "sysnumatom = ", sysnumatom, ", totnummol = ", totnummol
      write(*,*) "In COM mode, sysnumatom should equal to totnummol!"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    write(*,*) "sysnumatom=", sysnumatom

    allocate(vel(3, numframe, totnummol), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: vel"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    allocate(pos_tmp(3, sysnumatom), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_tmp"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    allocate(vel_tmp(3, sysnumatom), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: vel_tmp"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    allocate(time(numframe), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: time"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    numframe_read = 0
    call open_trajectory(datafileio, datafile)
    do i = 1, numframe
      do j = 1, skiptrr-1
        call read_trajectory(datafileio, sysnumatom, is_periodic, pos_tmp, vel_tmp, cell, tmp_r, stat)
        if (stat > 0) then
          write(*,*) "Reading trajectory error"
          call mpi_abort(MPI_COMM_WORLD, 1, ierr);
          call exit(1)
        else if (stat < 0) then
          !end of file
          exit
        end if
      end do
      call read_trajectory(datafileio, sysnumatom, is_periodic, pos_tmp, vel_tmp, cell, time(i), stat)
      if (stat /= 0) then
        write(*,*) "Reading trajectory error"
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
      numframe_read = numframe_read + 1
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
    call close_trajectory(datafileio)
    if (myrank == root) write(*,*) "numframe_read = ", numframe_read
    if (numframe_read /= numframe) then
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

  allocate(vv(numframe))
  if (stat /=0) then
    write(*,*) "Allocation failed: vv"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if

  allocate(ncorr(maxlag+1, num_moltypepair_all), stat=stat)
  if (stat /=0) then
    write(*,*) "Allocation failed: ncorr"
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end if
  ncorr = 0d0

  if (is_sd .or. is_ed) then
    if (is_sd) then
      call sd_cal_num_rBin(cell)
      call sd_prepCorrMemory(maxlag, nummoltype, numframe)
    end if
    if (is_ed) call ed_prepCorrMemory(maxlag, nummoltype, numframe)
  else
    allocate(corr_tmp(2*maxlag+1), stat=stat)
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
        moltypepair_allidx = getMolTypeIndex(i, sys%mol(:)%num, nummoltype)
        if (is_sd .or. is_ed) then
          do k = 1, maxlag+1
            numframe_k = numframe - k + 1
            vv(1:numframe_k) = sum(vel_r(:, k:numframe, i-r_start+1) * vel_c(:, 1:numframe_k, j-c_start+1), 1)
            ncorr(k, moltypepair_allidx) = ncorr(k, moltypepair_allidx) + sum(vv(1:numframe_k))
          end do
        else
          do k = 1, 3
            corr_tmp = corr(vel_r(k, :, i-r_start+1), maxlag)
            ncorr(:, moltypepair_allidx) = ncorr(:, moltypepair_allidx) + corr_tmp(maxlag+1:)
          end do
        end if
      else
        if (myrank == root) write(*,*) "loop r =",i-r_start+1, " of ", num_r,&
                                          ", c =", j-c_start+1, " of ", num_c
        starttime2 = MPI_Wtime()
        moltypepair_idx = getMolTypePairIndex(i, j, sys%mol(:)%num, nummoltype)
        moltypepair_allidx = moltypepair_idx + nummoltype
        if (is_sd .or. is_ed) then
          if (is_sd) call sd_getBinIndex(i-r_start+1, j-c_start+1, cell, sd_binIndex)
          if (is_ed) call ed_getBinIndex(i, j, ed_binIndex)
          do k = 1, maxlag+1
            numframe_k = numframe - k + 1
            vv(1:numframe_k) = sum(vel_r(:, k:numframe, i-r_start+1) * vel_c(:, 1:numframe_k, j-c_start+1), 1)
            !TODO: test if this sum should be put here or inside the following loop for better performance
            ncorr(k, moltypepair_allidx) = ncorr(k, moltypepair_allidx) + sum(vv(1:numframe_k))
            do n = 1, numframe_k
              if (is_sd) then
                tmp_i = sd_binIndex(n)
                if (tmp_i <= num_rBin) then
                  sdCorr(k, tmp_i, moltypepair_idx) = sdCorr(k, tmp_i, moltypepair_idx) + vv(n)
                end if
              end if
              if (is_ed) then
                tmp_i = ed_binIndex(n)
                if (tmp_i <= num_rBin) then
                  edCorr(k, tmp_i, moltypepair_idx) = edCorr(k, tmp_i, moltypepair_idx) + vv(n)
                end if
              end if
              !TODO: need test
              !ncorr(k, moltypepair_idx) = corr(k, moltypepair_idx) + vv(n)
            end do
          end do

          do t = 1, numframe
            if (is_sd) then
              tmp_i = sd_binIndex(t)
              if (tmp_i <= num_rBin) then
                sdPairCount(tmp_i, moltypepair_idx) = sdPairCount(tmp_i, moltypepair_idx) + 1d0
              end if
            end if
            if (is_ed) then
              tmp_i = ed_binIndex(t)
              if (tmp_i <= num_eBin) then
                edPairCount(tmp_i, moltypepair_idx) = edPairCount(tmp_i, moltypepair_idx) + 1d0
              end if
            end if
          end do

        else ! one-two only
          do k = 1, 3
            corr_tmp = corr(vel_r(k, :, i-r_start+1), vel_c(k, :, j-c_start+1), maxlag)
            ncorr(:, moltypepair_allidx) = ncorr(:, moltypepair_allidx) + corr_tmp(maxlag+1:)
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

  !collect ncorr
  if (myrank == root) write(*,*) "start collecting results"
  starttime = MPI_Wtime()
  if (myrank == root) then
    write(*,*) "collecting ncorr"
    call mpi_reduce(MPI_IN_PLACE, ncorr, size(ncorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  else
    call mpi_reduce(ncorr, dummy_null, size(ncorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
  end if
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (is_sd) call sd_collectCorr()
  if (is_ed) call ed_collectCorr()
  endtime = MPI_Wtime()
  if (myrank == root) write(*,*) "finished collecting results. It took ", endtime - starttime, " seconds"

  !average at root and then output
  if (myrank == root) then
    allocate(framecount(maxlag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: framecount"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    do j = 1, nummoltype
      do i = j, nummoltype
        if (i /= j) then
          moltypepair_idx = get_pairindex_upper_diag(i, j, nummoltype)
          moltypepair_allidx = moltypepair_idx + nummoltype
          ncorr(:, moltypepair_allidx) = ncorr(:, moltypepair_allidx) / 2d0
        end if
      end do
    end do

    framecount = [ (numframe - (i-1), i = 1, maxlag+1) ] * 3d0
    do n = 1, num_moltypepair_all
      ncorr(:,n) = ncorr(:,n) / framecount
    end do

    if (is_sd) call sd_average(numframe, nummoltype, framecount)
    if (is_ed) call ed_average(numframe, nummoltype, framecount)

    deallocate(framecount)

    !output results
    write(*,*) "start writing outputs"
    starttime = MPI_Wtime()
    call output()
    endtime = MPI_Wtime()
    write(*,*) "finished writing outputs. It took ", endtime - starttime, " seconds"
  end if

  if (myrank == root) write(*,*)
  if (myrank == root) write(*,*) "time for the whole program (sec):", MPI_Wtime() - prog_starttime

  deallocate(ncorr, charge, vel_r, vel_c, start_index)
  if (is_sd) call sd_finish()
  if (is_ed) call ed_finish()
  call mpi_setup('stop')
  stop

contains
  integer function getMolTypeIndex(i, numMol, nummoltype)
    implicit none
    integer, intent(in) :: i, numMol(:), nummoltype
    integer :: n, numMol_acc

    getMolTypeIndex = -1
    numMol_acc = 0
    do n = 1, nummoltype
      numMol_acc = numMol_acc + numMol(n)
      if (i <= numMol_acc) then
        getMolTypeIndex = n
        return
      end if
    end do
  end function getMolTypeIndex

  integer function getMolTypePairIndex(i, j, numMol, nummoltype)
    implicit none
    integer, intent(in) :: i, j, numMol(:), nummoltype
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


    ii = getMolTypeIndex(i, numMol, nummoltype)
    jj = getMolTypeIndex(j, numMol, nummoltype)
    r = min(ii, jj)
    c = max(ii, jj)
    getMolTypePairIndex = (r - 1) * nummoltype + c - r * (r - 1) / 2
  end function getMolTypePairIndex

  subroutine com_vel(com_v, vel, start_index, sys)
    implicit none
    real(dp), dimension(:, :), intent(out) :: com_v
    real(dp), dimension(:, :), intent(in) :: vel
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
    real(dp), allocatable :: timeLags(:)
    integer :: ierr
    integer(hid_t) :: dset_ncorr, dset_timeLags, &
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
                                   DSETNAME_NCORR = "ncorr"

    !/GROUP_SPATIAL or GROUP_ENERGY/Dataset
    character(len=*), parameter :: DSETNAME_DECBINS = "decBins", &
                                   DSETNAME_DECCORR = "decCorr", &
                                   DSETNAME_DECPAIRCOUNT = "decPairCount"


    allocate(timeLags(maxlag+1), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: timeLags"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    timeLags = [ (dble(i), i = 0, maxlag) ] * timestep

    if (is_sd) then
      call sd_make_rBins()
    end if

    if (is_ed) then
      call ed_make_eBins()
    end if

    !create and write attributes
    call H5LTset_attribute_string_f(corrfileio, GROUP_ROOT, ATTR_VERSION, decond_version, ierr)
    call H5LTset_attribute_string_f(corrfileio, GROUP_ROOT, ATTR_TYPE, OUT_TYPE, ierr)

    !create and write datasets
    !volume
    call H5Screate_f(H5S_SCALAR_F, space_id, ierr)
    call H5Dcreate_f(corrfileio, DSETNAME_VOLUME, H5T_NATIVE_DOUBLE, space_id, dset_id, ierr)
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, product(cell), [0_hsize_t], ierr)
    call H5Dclose_f(dset_id, ierr)
    call H5Sclose_f(space_id, ierr)
    call H5LTset_attribute_string_f(corrfileio, DSETNAME_VOLUME, ATTR_UNIT, "nm$^3$", ierr)

    !temperature
    call H5Screate_f(H5S_SCALAR_F, space_id, ierr)
    call H5Dcreate_f(corrfileio, DSETNAME_TEMP, H5T_NATIVE_DOUBLE, space_id, dset_id, ierr)
    call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, temperature, [0_hsize_t], ierr)
    call H5Dclose_f(dset_id, ierr)
    call H5Sclose_f(space_id, ierr)
    call H5LTset_attribute_string_f(corrfileio, DSETNAME_TEMP, ATTR_UNIT, "K", ierr)

    !charge
    call H5LTmake_dataset_int_f(corrfileio, DSETNAME_CHARGE, 1, [size(charge, kind=hsize_t)], charge, ierr)
    call H5LTset_attribute_string_f(corrfileio, DSETNAME_CHARGE, ATTR_UNIT, "e", ierr)

    !numMol
    call H5LTmake_dataset_int_f(corrfileio, DSETNAME_NUMMOL, 1, &
                                   [size(sys%mol(:)%num, kind=size_t)], sys%mol(:)%num, ierr)

    !timeLags
    call H5LTmake_dataset_double_f(corrfileio, DSETNAME_TIMELAGS, 1, [size(timeLags, kind=hsize_t)], timeLags, ierr)
    call H5Dopen_f(corrfileio, DSETNAME_TIMELAGS, dset_timeLags, ierr)
    call H5LTset_attribute_string_f(corrfileio, DSETNAME_TIMELAGS, ATTR_UNIT, "ps", ierr)

    !ncorr
    call H5LTmake_dataset_double_f(corrfileio, DSETNAME_NCORR, 2, &
        [size(ncorr, 1, kind=hsize_t), size(ncorr, 2, kind=hsize_t)], ncorr, ierr)
    call H5Dopen_f(corrfileio, DSETNAME_NCORR, dset_ncorr, ierr)
    call H5LTset_attribute_string_f(corrfileio, DSETNAME_NCORR, ATTR_UNIT, "nm$^2$ ps$^{-2}$", ierr)

    if (is_sd) then
      !create a group for storing spatial-decomposition data
      call H5Gcreate_f(corrfileio, GROUP_SPATIAL, grp_sd_id, ierr)

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
      call H5Gcreate_f(corrfileio, GROUP_ENERGY, grp_ed_id, ierr)

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
    call H5DSattach_scale_f(dset_ncorr, dset_timeLags, 2, ierr)
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
    call H5Dclose_f(dset_ncorr, ierr)
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
    call H5Fclose_f(corrfileio, ierr)
    call H5close_f(ierr)
  end subroutine output

  integer function count_arg(i, num_arg)
    implicit none
    integer, intent(in) :: i, num_arg
    character(len=line_len) :: arg
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

  real(dp) function getTfromLog(logfile)
    implicit none
    character(len=*), intent(in) :: logfile
    integer :: logio, idx, temp_idx, stat
    integer, parameter :: RECORD_LEN = 15
    character(len=line_len) :: line
    logical :: found_average

    found_average = .false.
    open(unit=newunit(logio), file=logfile, status='old', action="READ", form="FORMATTED")
    do while(.true.)
      read(logio, "(A"//line_len_str//")", iostat=stat) line
      if (stat > 0) then
        write(*,*) "Error reading line"
        call exit(1)
      else if (stat < 0) then
        write(*,*) "Unable to find 'A V E R A G E S' in logfile ", trim(adjustl(logfile))
        call exit(1)
      end if

      if (found_average) then
        idx = index(line, "Temperature")
        if (idx > 0) then
          temp_idx = ceiling(real(idx, dp) / RECORD_LEN)
          read(logio, "(A"//line_len_str//")", iostat=stat) line
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
    write(*, *) "  -ed <engtraj> <engtraj> ...: do energy decomposition. default no ed."
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
