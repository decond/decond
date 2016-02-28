module energy_dec
  use mpiproc
  use hdf5
  use varpars, only: numframe, totnummol
  implicit none
  private
  integer, parameter :: line_len = 1024

  public ed_getbinindex, ed_prep_corrmemory, ed_collectcorr, &
         ed_average, ed_init, ed_make_ebins, ed_finish, ed_prep
  character(len=line_len), public, allocatable :: engfiles(:)
  real(8), public, allocatable :: ed_binIndex(:)  !ed_binIndex(numframe)
  integer, public :: num_ebin, skipEng
  real(8), public, allocatable :: edpaircount(:, :)
  !edpaircount(num_ebin, nummoltype*nummoltype)
  real(8), public, allocatable :: edcorr(:, :, :)
  real(8), public :: ebinwidth
  integer, public :: num_engfiles
  real(8), public, allocatable :: eBins(:)  !eBins(num_ebin)
  real(8), public :: engMin_global

  integer, parameter :: MIN_ENGTRJ_VER_MAJOR = 0
  character(len=*), parameter :: ENGDSET_NAME = "energy"
  integer(hid_t), allocatable :: engfileids(:)
  real(8), allocatable :: eBinIndexAll(:, :)  !eBinIndexAll(numframe,
                                              !             uniqueNumMolPair)
  real(8) :: engMax_global
  integer, allocatable :: engLocLookupTable(:, :)
  integer :: eBinIndex_absolute_max, eBinIndex_absolute_min
  integer, allocatable :: sltspec_list(:), sltfirsttag_list(:), &
                         &nummol_list(:), numpair_list(:), &
                         &numframe_list(:), numslt_list(:)
contains
  subroutine ed_init()
    implicit none
    ebinwidth = 0.5
    skipEng = 1

    ! initialize HDF5 Fortran predefined datatypes
    call h5open_f(ierr)
  end subroutine ed_init

  subroutine ed_prep()
    call prep_engfiles()

    ! prepare eBinIndex for each rank
    if (myrank == root) then
      write(*,*) "start preparing eBinIndex..."
      starttime = MPI_Wtime()
    end if
    call readEng()
    if (myrank == root) then
      endtime = MPI_Wtime()
      write(*,*) "finished preparing eBinIndex. It took ", endtime - starttime, "seconds"
    end if
  end subroutine ed_prep

  subroutine prep_engfiles()
    allocate(engfiles(num_engfiles))
    allocate(engfileids(num_engfiles))
    allocate(sltfirsttag_list(num_engfiles))

    if (myrank == root) then
      allocate(sltspec_list(num_engfiles))
      allocate(numslt_list(num_engfiles))
      allocate(nummol_list(num_engfiles))
      allocate(numpair_list(num_engfiles))
      allocate(numframe_list(num_engfiles))
    end if
  end subroutine prep_engfiles

  subroutine check_version()
    use utility, only: parse_version
    implicit none
    integer(hsize_t), parameter :: VERSION_LEN = 11
    integer(hid_t) :: file_id, attr_id, type_id
    integer(hsize_t) :: dims(1)
    character(len=VERSION_LEN) :: buf
    integer :: i, major, minor, patch

    dims = -1  ! value of dims does not seem to have any effect
    do i = 1, num_engfiles
      call h5fopen_f(engfiles(i), H5F_ACC_RDONLY_F, file_id, ierr)

      call h5aopen_f(file_id, 'version', attr_id, ierr)
      call h5aget_type_f(attr_id, type_id, ierr)
      call h5aread_f(attr_id, type_id, buf, dims, ierr)
      call h5aclose_f(attr_id, ierr)
      call h5fclose_f(file_id, ierr)
      call parse_version(buf, major, minor, patch)
      if (major < MIN_ENGTRJ_VER_MAJOR) then
        write(*,*) "Only engtraj with a major version number greater than &
                   &or equal to ", MIN_ENGTRJ_VER_MAJOR, " is supported."
        write(*,*) "The major version of the file '", trim(engfiles(i)), &
                   " is: ", major
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
    end do
  end subroutine check_version

  subroutine read_attr_integer(file_id, attr_name, buf)
    implicit none
    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: attr_name
    integer, intent(out) :: buf
    integer(hid_t) :: attr_id
    integer(hsize_t) :: dims(1)

    dims = -1  ! value of dims does not seem to have any effect
    call h5aopen_f(file_id, attr_name, attr_id, ierr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, ierr)
    call h5aclose_f(attr_id, ierr)
  end subroutine read_attr_integer

  subroutine sort_engfiles()
    use utility, only: swap
    implicit none
    integer :: i, j, buf
    integer(hid_t) :: file_id

    do i = 1, num_engfiles
      call h5fopen_f(engfiles(i), H5F_ACC_RDONLY_F, file_id, ierr)
      call read_attr_integer(file_id, 'sltspec', buf)
      sltspec_list(i) = buf
      call h5fclose_f(file_id, ierr)
    end do

    ! sort engfiles according to sltspec
    do i = 1, num_engfiles
      do j = i, num_engfiles
        if (sltspec_list(j) == i) then
          if (j /= i) then
            call swap(sltspec_list, i, j)
            call swap(engfiles, i, j)
          end if
          exit
        end if
      end do
    end do
  end subroutine sort_engfiles

  subroutine read_attrs()
    implicit none
    integer :: i
    integer(hid_t) :: file_id, attr_id
    integer(hsize_t) :: dims(1)
    integer :: buf

    do i = 1, num_engfiles
      call h5fopen_f(engfiles(i), H5F_ACC_RDONLY_F, file_id, ierr)

      call read_attr_integer(file_id, 'sltfirsttag', buf)
      sltfirsttag_list(i) = buf

      call read_attr_integer(file_id, 'numslt', buf)
      numslt_list(i) = buf

      call read_attr_integer(file_id, 'nummol', buf)
      nummol_list(i) = buf

      call read_attr_integer(file_id, 'numpair', buf)
      numpair_list(i) = buf

      call read_attr_integer(file_id, 'numframe', buf)
      numframe_list(i) = buf

      call h5fclose_f(file_id, ierr)
    end do
  end subroutine read_attrs

  subroutine check_sanity()
    implicit none
    integer :: i

    do i = 1, num_engfiles
      if (sltspec_list(i) /= i) then
        write(*,*) "sltspec_list is not sorted correctly. it should consist of &
                   &continuous integers starting from 1"
        write(*,*) "sltspec_list =", sltspec_list
        call exit(1)
      end if
    end do

    do i = 1, num_engfiles - 1
      if (sltfirsttag_list(i) >= sltfirsttag_list(i+1)) then
        write(*,*) "sltfirsttag_list should be an ascending list"
        write(*,*) "sltfirsttag_list =", sltfirsttag_list
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
    end do

    if (.not. all(nummol_list == nummol_list(1))) then
      write(*,*) "nummol's should be all the same"
      write(*,*) "nummol_list =", nummol_list
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    if (totnummol /= nummol_list(1)) then
      write(*,*) "nummol inconsistent:"
      write(*,*) "nummol from topology:", totnummol
      write(*,*) "nummol from energy trajectories =", nummol_list
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    if (.not. all(numframe_list == numframe_list(1))) then
      write(*,*) "numframe's should be all the same"
      write(*,*) "numframe_list =", numframe_list
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    if (sum(numslt_list) /= nummol_list(1)) then
      write(*,*) "sum of numslt's should equal to nummol"
      write(*,*) "numslt_list =", numslt_list
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    if (sum(numpair_list) /= nummol_list(1) * (nummol_list(1) - 1) / 2) then
      write(*,*) "numpair is not consistent with nummol"
      write(*,*) "numpair_list =", numpair_list
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end subroutine check_sanity

  subroutine bcast_attrs()
    implicit none
    call mpi_bcast(engfiles, num_engfiles * line_len, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
    call mpi_bcast(sltfirsttag_list, num_engfiles, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
  end subroutine bcast_attrs

  subroutine readEng()
    implicit none
    integer :: i, r, c, loc, lastLoc, locMax, stat, begin_time
    integer, allocatable :: eBinIndex_single(:)
    real(8), allocatable :: eng(:)
    real(8) :: engMax, engMin, engMax_node, engMin_node
    logical :: isFirstRun

    if (myrank == root) then
      call check_version()
      call sort_engfiles()
      call read_attrs()
      call check_sanity()
    end if

    call bcast_attrs()
    call openEngtraj()
    call makeEngPairLookupTable(r_start, r_end, c_start, c_end, locMax)
    allocate(eBinIndexAll(numframe, locMax), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: eBinIndexAll"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    allocate(eng(numframe), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: eng"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    allocate(eBinIndex_single(numframe), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: eBinIndex_single"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    ! determine max and min of engtraj records
    if (myrank == root) then
      write(*,*) "determining min and max of engtraj"
      begin_time = MPI_Wtime()
    end if
    engMin_node = 1.0
    engMax_node = 1.0
    lastLoc = 0
    isFirstRun = .true.
    do c = c_start, c_end
      do r = r_start, r_end
        if (r /= c) then
          loc = engLocLookupTable(r, c)
          if (loc > lastLoc) then
            ! new location, read new data
            call readPairEng(eng, r, c, numframe, totnummol)
            engMin = minval(eng)
            engMax = maxval(eng)
            if (isFirstRun) then
              engMin_node = engMin
              engMax_node = engMax
              isFirstRun = .false.
            else
              if (engMin < engMin_node) engMin_node = engMin
              if (engMax > engMax_node) engMax_node = engMax
            end if
            lastLoc = loc
          end if
        end if
      end do
    end do
    if (myrank == root) write(*,*) "time for determining min and max (sec):", &
                                   &MPI_Wtime() - begin_time
    call mpi_allreduce(engMin_node, engMin_global, 1, MPI_DOUBLE_PRECISION, &
                      &MPI_MIN, MPI_COMM_WORLD, ierr)
    call mpi_allreduce(engMax_node, engMax_global, 1, MPI_DOUBLE_PRECISION, &
                      &MPI_MAX, MPI_COMM_WORLD, ierr)

    if (myrank == root) then
      write(*, *) "energy min:", engMin_global, " max:", engMax_global
    end if

    eBinIndex_absolute_max = getAbsoluteIndex(engMax_global)
    eBinIndex_absolute_min = getAbsoluteIndex(engMin_global)
    num_ebin = eBinIndex_absolute_max - eBinIndex_absolute_min + 1

    if (myrank == root) write(*,*) "reading engtraj data"
    lastLoc = 0
    do c = c_start, c_end
      do r = r_start, r_end
        if (r /= c) then
          loc = engLocLookupTable(r, c)
          if (loc > lastLoc) then
            ! new loc, read new data
            call readPairEng(eng, r, c, numframe, totnummol)
            call eng2BinIndex(eng, eBinIndex_single)
            eBinIndexAll(:, loc) = eBinIndex_single
            lastLoc = loc
          end if
        end if
      end do
    end do

    deallocate(eng, eBinIndex_single)

    do i = 1, num_engfiles
      call H5Fclose_f(engfileids(i), ierr)
    end do
  end subroutine readEng

  subroutine ed_prep_corrmemory(maxlag, nummoltype, numframe)
    implicit none
    integer, intent(in) :: maxlag, nummoltype, numframe
    integer :: stat
    integer :: num_moltypepair

    num_moltypepair = nummoltype * (nummoltype + 1) / 2
    allocate(edcorr(maxlag+1, num_ebin, num_moltypepair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: edcorr"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    edcorr = 0d0

    allocate(edpaircount(num_ebin, num_moltypepair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: edpaircount"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    edpaircount = 0d0

    allocate(ed_binIndex(numframe), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: ed_binIndex"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end subroutine ed_prep_corrmemory

  subroutine openEngtraj()
    use utility, only: parse_version
    use H5LT
    implicit none
    character(len=11) :: engtrajVer
    integer :: engtrajVer_major, numpair
    integer(hid_t) :: plist_id  ! property list identifier for HDF5
    integer :: buf(1), i

    ! setup file access property list with parallel I/O access.
    call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)

    do i = 1, num_engfiles
      ! open the existing engtraj file.
      call H5Fopen_f(engfiles(i), H5F_ACC_RDONLY_F, engfileids(i), ierr, &
                    &access_prp = plist_id)
      if (ierr /= 0) then
        write(*,*) "Failed to open HDF5 file: ", trim(adjustl(engfiles(i)))
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
    end do

    ! close property list
    call H5Pclose_f(plist_id, ierr)
  end subroutine openEngtraj

  subroutine readPairEng(eng, r, c, numframe, nummol)
    use hdf5
    use utility, only: get_pairindex_upper_nodiag
    implicit none
    real(8), intent(out) :: eng(numframe)
    integer, intent(in) :: r, c, numframe, nummol
    integer(hid_t) :: engfile_id
    integer(hid_t) :: dset_id
    integer(hid_t) :: filespace     ! Dataspace identifier in file
    integer(hid_t) :: memspace      ! Dataspace identifier in memory
    integer(hid_t) :: plist_id      ! Property list identifier
    integer(hsize_t) :: offset(2), count(2), stride(2)
    integer :: pairindex, fileidx, first_pairindex_in_file, rel_pairindex

    fileidx = getmoltype(min(r, c))
    engfile_id = engfileids(fileidx)
    first_pairindex_in_file = get_pairindex_upper_nodiag(&
        &sltfirsttag_list(fileidx), sltfirsttag_list(fileidx)+1, nummol)

    ! open dataset in file and get the filespace (dataspace)
    call H5Dopen_f(engfile_id, ENGDSET_NAME, dset_id, ierr)
    call H5Dget_space_f(dset_id, filespace, ierr)

    ! select data slab (column) for pair (r, c) in filespace
    !          c
    !    | 1  2  3  4
    !  --+------------
    !  1 |    1  2  3
    !    |
    !  2 |       4  5
    !r   |
    !  3 |          6
    !    |
    !  4 |
    !
    !  pairindex(r, c) = (r - 1) * n + c - (r + 1) * r / 2
    !  n = 4 in this example
    pairindex = get_pairindex_upper_nodiag(r, c, nummol)
    rel_pairindex = pairindex - first_pairindex_in_file + 1
    offset = [skipEng - 1, rel_pairindex - 1]
    count = [numframe, 1]
    stride = [skipEng, 1]
    call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
                              &ierr, stride)

    ! create memory space
    call H5Screate_simple_f(1, [int(numframe, kind=hsize_t)], memspace, ierr)

    ! Create property list for reading dataset
    call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, eng, &
                  &[int(numframe, kind=hsize_t)], ierr, &
                  &file_space_id = filespace, mem_space_id = memspace, &
                  &xfer_prp = plist_id)

    call H5Pclose_f(plist_id, ierr)

    call H5Sclose_f(filespace, ierr)
    call H5Sclose_f(memspace, ierr)
    call H5Dclose_f(dset_id, ierr)
  end subroutine readPairEng

  ! index = 0 : self-interaction energy (ignore this record)
  ! index = N : N is the location of the stored data (in the eBinIndexAll(:, N))
  subroutine makeEngPairLookupTable(r_start, r_end, c_start, c_end, locMax)
    implicit none
    integer, intent(in) :: r_start, r_end, c_start, c_end
    integer, intent(out) :: locMax
    integer :: loc, pairHash, cacheIndex, stat, r, c
    integer :: pairHashCache((r_end - r_start + 1) * (c_end - c_start + 1))

    allocate(engLocLookupTable(r_start:r_end, c_start:c_end), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: engLocLookupTable"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    pairHashCache = -1
    engLocLookupTable = 0
    loc = 0
    cacheIndex = 0

    do c = c_start, c_end
      do r = r_start, r_end
        if (r /= c) then
          pairHash = getUnorderedCantorPairHash(r, c)
          ! check if the new pairHash exists in the cache
          if (.not. any(pairHashCache == pairHash)) then
            ! this is a new pair
            ! increase the counters
            cacheIndex = cacheIndex + 1
            loc = loc + 1

            ! make a engLocLookupTable so that we can know
            ! where the data for a certain pair are stored later
            pairHashCache(cacheIndex) = pairHash
            engLocLookupTable(r, c) = loc
          else
            ! this is an old pair, set the index to the previous symmetric one
            engLocLookupTable(r, c) = engLocLookupTable(c, r)
          end if
        end if
      end do
    end do

    locMax = loc
  end subroutine makeEngPairLookupTable

  integer function getCantorPairHash(k1, k2)
    ! ref: https://en.wikipedia.org/wiki/Pairing_function#/Cantor_pairing_function
    implicit none
    integer, intent(in) :: k1, k2
    integer :: k

    k = k1 + k2
    getCantorPairHash = k * (k + 1) / 2 + k2
  end function getCantorPairHash

  integer function getUnorderedCantorPairHash(k1, k2)
    implicit none
    integer, intent(in) :: k1, k2

    if (k1 <= k2) then
      getUnorderedCantorPairHash = getCantorPairHash(k1, k2)
    else
      getUnorderedCantorPairHash = getCantorPairHash(k2, k1)
    end if
  end function getUnorderedCantorPairHash

  !Absolute energy indexes are centered at 0
  !The bin grid is totally determined by ebinwidth
  integer elemental function getAbsoluteIndex(eng)
    implicit none
    real(8), intent(in) :: eng
      getAbsoluteIndex = floor(eng / ebinwidth + 0.5)
  end function getAbsoluteIndex

  subroutine eng2BinIndex(eng, eBinIndex_single)
    implicit none
    real(8), intent(in) :: eng(:)
    integer, intent(out) :: eBinIndex_single(size(eng))

    eBinIndex_single = getAbsoluteIndex(eng) - eBinIndex_absolute_min + 1
  end subroutine eng2BinIndex

  subroutine ed_getbinindex(r, c, eBin)
    implicit none
    integer, intent(in) :: r, c
    real(8), intent(out) :: eBin(:)

    eBin = eBinIndexAll(:, engLocLookupTable(r, c))
  end subroutine ed_getbinindex

  subroutine ed_collectcorr()
    implicit none
    if (myrank == root) then
      call mpi_reduce(MPI_IN_PLACE, edcorr, size(edcorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    else
      call mpi_reduce(edcorr, dummy_null, size(edcorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    end if
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (myrank == root) then
      call mpi_reduce(MPI_IN_PLACE, edpaircount, size(edpaircount), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    else
      call mpi_reduce(edpaircount, dummy_null, size(edpaircount), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    end if
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end subroutine ed_collectcorr

  subroutine ed_average(numframe, nummoltype, framecount)
    use utility, only: get_pairindex_upper_diag
    implicit none
    integer, intent(in) :: numframe, nummoltype, framecount(:)
    integer :: i, t1, t2, n, moltypepair_idx

    edpaircount = edpaircount / numframe
    do n = 1, nummoltype * (nummoltype + 1) / 2
      do i = 1, num_ebin
        edcorr(:,i,n) = edcorr(:,i,n) / framecount / edpaircount(i, n)
      end do
    end do

    do t2 = 1, nummoltype
      do t1 = t2, nummoltype
        if (t1 /= t2) then
          moltypepair_idx = get_pairindex_upper_diag(t1, t2, nummoltype)
          edpaircount(:, moltypepair_idx) = edpaircount(:, moltypepair_idx) / 2d0
        end if
      end do
    end do
  end subroutine ed_average

  subroutine ed_make_ebins()
    implicit none
    integer :: i, stat
    allocate(eBins(num_ebin), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: eBins"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    eBins = [(i * ebinwidth, i = eBinIndex_absolute_min, eBinIndex_absolute_max)]
  end subroutine ed_make_ebins

  subroutine ed_finish()
    implicit none
    deallocate(engfiles, engfileids, sltfirsttag_list)
    deallocate(eBinIndexAll, engLocLookupTable, edcorr)
    if (myrank == root) then
      deallocate(sltspec_list, numslt_list, nummol_list)
      deallocate(numpair_list, numframe_list)
    end if
    ! ed_binIndex has been deallocated in the main program
  end subroutine ed_finish

  integer function getmoltype(molidx)
    implicit none
    integer, intent(in) :: molidx
    integer :: i

    do i = size(sltfirsttag_list), 1, -1
      if (molidx >= sltfirsttag_list(i)) then
        getmoltype = i
        return
      end if
    end do
    write(*,*) "Unable to determine moltype for index:", molidx
    call mpi_abort(MPI_COMM_WORLD, 1, ierr);
    call exit(1)
  end function getmoltype
end module energy_dec
