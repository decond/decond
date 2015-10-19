module energy_dec
  use mpiproc
  implicit none
  integer :: num_eBin, skipEng
  integer, parameter :: MIN_ENGTRJ_VER_MAJOR = 0
  character(len=*), parameter :: ENGDSET_NAME = "energy"
  character(len=128) :: engtrajFilename
  real(8), allocatable :: eBinIndexAll(:, :)  !eBinIndexAll(numFrame, uniqueNumMolPair)
  real(8), allocatable :: ed_binIndex(:)  !ed_binIndex(numFrame)
  real(8), allocatable :: edCorr(:, :, :)  !edCorr(maxLag+1, num_eBin, numMolType*numMolType)
  real(8), allocatable :: edPairCount(:, :)  !edPairCount(num_eBin, numMolType*numMolType)
  real(8), allocatable :: eBins(:)  !eBins(num_eBin)
  real(8) :: engMin_global, engMax_global
  integer, allocatable :: engLocLookupTable(:, :)
  real(8) :: eBinWidth
  integer :: nummol, eBinIndex_absolute_max, eBinIndex_absolute_min

contains
  subroutine ed_init()
    implicit none
    eBinWidth = 0.5
    skipEng = 1
  end subroutine ed_init

  subroutine ed_readEng(engtrajFilename, numFrame)
    use HDF5
    implicit none
    character(len=*), intent(in) :: engtrajFilename
    integer, intent(in) :: numFrame
    integer(hid_t) :: engtrajFileid
    integer :: r, c, loc, lastLoc, locMax, stat, begin_time
    integer, allocatable :: eBinIndex_single(:)
    real(8), allocatable :: eng(:)
    real(8) :: engMax, engMin, engMax_node, engMin_node
    logical :: isFirstRun

    call openEngtraj(engtrajFileid)
    call makeEngPairLookupTable(r_start, r_end, c_start, c_end, locMax)
    allocate(eBinIndexAll(numFrame, locMax), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: eBinIndexAll"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    allocate(eng(numFrame), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: eng"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    allocate(eBinIndex_single(numFrame), stat=stat)
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
            call readPairEng(eng, r, c, engtrajFileid, numFrame)
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
    if (myrank == root) write(*,*) "time for determining min and max (sec):", MPI_Wtime() - begin_time
    call mpi_allreduce(engMin_node, engMin_global, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call mpi_allreduce(engMax_node, engMax_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    if (myrank == root) write(*, *) "energy min:", engMin_global, " max:", engMax_global

    eBinIndex_absolute_max = getAbsoluteIndex(engMax_global)
    eBinIndex_absolute_min = getAbsoluteIndex(engMin_global)
    num_eBin = eBinIndex_absolute_max - eBinIndex_absolute_min + 1

    if (myrank == root) write(*,*) "reading engtraj data"
    lastLoc = 0
    do c = c_start, c_end
      do r = r_start, r_end
        if (r /= c) then
          loc = engLocLookupTable(r, c)
          if (loc > lastLoc) then
            ! new loc, read new data
            call readPairEng(eng, r, c, engtrajFileid, numFrame)
            call eng2BinIndex(eng, eBinIndex_single)
            eBinIndexAll(:, loc) = eBinIndex_single
            lastLoc = loc
          end if
        end if
      end do
    end do

    deallocate(eng, eBinIndex_single)
    call H5Fclose_f(engtrajFileid, ierr)
  end subroutine ed_readEng

  subroutine ed_prepCorrMemory(maxLag, numMolType, numFrame)
    implicit none
    integer, intent(in) :: maxLag, numMolType, numFrame
    integer :: stat
    integer :: numMolTypePair

    numMolTypePair = numMolType * (numMolType + 1) / 2
    allocate(edCorr(maxLag+1, num_eBin, numMolTypePair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: edCorr"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    edCorr = 0d0

    allocate(edPairCount(num_eBin, numMolTypePair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: edPairCount"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    edPairCount = 0d0

    allocate(ed_binIndex(numFrame), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: ed_binIndex"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
  end subroutine ed_prepCorrMemory

  subroutine openEngtraj(engtrajFileid)
    use utility, only: parse_version
    use HDF5
    use H5LT
    implicit none
    integer(hid_t), intent(out) :: engtrajFileid
    character(len=11) :: engtrajVer
    integer :: engtrajVer_major, numpair
    integer(hid_t) :: plist_id  ! property list identifier for HDF5
    integer :: buf(1)

    ! initialize HDF5 Fortran predefined datatypes
    call H5open_f(ierr)

    ! setup file access property list with parallel I/O access.
    call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)

    ! open the existing engtraj file.
    call H5Fopen_f(engtrajFilename, H5F_ACC_RDONLY_F, engtrajFileid, ierr, access_prp = plist_id)
    if (ierr /= 0) then
      write(*,*) "Failed to open HDF5 file: ", trim(adjustl(engtrajFilename))
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    ! close property list
    call H5Pclose_f(plist_id, ierr)

    ! read version
    call H5LTget_attribute_string_f(engtrajFileid, "/", "version", engtrajVer, ierr)
    ! check compatibility at root
    if (myrank == root) then
      call parse_version(engtrajVer, engtrajVer_major)
      if (engtrajVer_major < MIN_ENGTRJ_VER_MAJOR) then
        write(*,*) "Only engtraj with a major version number greater than ", MIN_ENGTRJ_VER_MAJOR, " is supported."
        write(*,*) "The major version of the file '", trim(engtrajFilename), " is: ", engtrajVer_major
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
    end if

    ! check data consistency
    call H5LTget_attribute_int_f(engtrajFileid, "/", "nummol", buf, ierr)
    nummol = buf(1)
    call H5LTget_attribute_int_f(engtrajFileid, "/", "numpair", buf, ierr)
    numpair = buf(1)
    if (myrank == root) then
      if ((nummol * (nummol - 1) / 2) /= numpair) then
        write(*,*) "Error: engtraj data inconsistent. numpair should be equal to nummol*(nummol-1)/2"
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
    end if
  end subroutine openEngtraj

  subroutine readPairEng(eng, r, c, engtrajFileid, numFrame)
    use HDF5
    use H5LT
    implicit none
    real(8), intent(out) :: eng(numFrame)
    integer, intent(in) :: r, c, numFrame
    integer(hid_t), intent(in) :: engtrajFileid
    integer(hid_t) :: dset_id
    integer(hid_t) :: filespace     ! Dataspace identifier in file
    integer(hid_t) :: memspace      ! Dataspace identifier in memory
    integer(hid_t) :: plist_id      ! Property list identifier
    integer(hsize_t) :: offset(2), count(2), stride(2)
    integer :: pairIndex

    ! open dataset in file and get the filespace (dataspace)
    call H5Dopen_f(engtrajFileid, ENGDSET_NAME, dset_id, ierr)
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
    !  pairIndex(r, c) = (r - 1) * n + c - (r + 1) * r / 2
    !  n = 4 in this example
    pairIndex = getPairIndex(min(r, c), max(r, c), nummol)
    offset = [skipEng - 1, pairIndex - 1]
    count = [numFrame, 1]
    stride = [skipEng, 1]
    call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, ierr, stride)

    ! create memory space
    call H5Screate_simple_f(1, [int(numFrame, kind=hsize_t)], memspace, ierr)

    ! Create property list for reading dataset
    call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, ierr)
    call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, ierr)

    call H5Dread_f(dset_id, H5T_NATIVE_DOUBLE, eng, [int(numFrame, kind=hsize_t)], ierr, &
                    file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    call H5Pclose_f(plist_id, ierr)

    call H5Sclose_f(filespace, ierr)
    call H5Sclose_f(memspace, ierr)
    call H5Dclose_f(dset_id, ierr)
  end subroutine readPairEng

  integer function getPairIndex(r, c, n)
    implicit none
    integer, intent(in) :: r, c, n

    if (r < 0 .or. c < 0 .or. n < 0 .or. r >= c .or. r > n .or. c > n) then
      write(*,*) "Error: unreasonable parameters in getPairIndex(r, c, n): ", r, c, n
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    getPairIndex = (r - 1) * n + c - ((r + 1) * r) / 2
  end function getPairIndex

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
  !The bin grid is totally determined by eBinWidth
  integer elemental function getAbsoluteIndex(eng)
    implicit none
    real(8), intent(in) :: eng
      getAbsoluteIndex = floor(eng / eBinWidth + 0.5)
  end function getAbsoluteIndex

  subroutine eng2BinIndex(eng, eBinIndex_single)
    implicit none
    real(8), intent(in) :: eng(:)
    integer, intent(out) :: eBinIndex_single(size(eng))

    eBinIndex_single = getAbsoluteIndex(eng) - eBinIndex_absolute_min + 1
  end subroutine eng2BinIndex

  subroutine ed_getBinIndex(r, c, eBin)
    implicit none
    integer, intent(in) :: r, c
    real(8), intent(out) :: eBin(:)

    eBin = eBinIndexAll(:, engLocLookupTable(r, c))
  end subroutine ed_getBinIndex

  subroutine ed_collectCorr()
    implicit none
    if (myrank == root) then
      call mpi_reduce(MPI_IN_PLACE, edCorr, size(edCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    else
      call mpi_reduce(edCorr, dummy_null, size(edCorr), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    end if
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (myrank == root) then
      call mpi_reduce(MPI_IN_PLACE, edPairCount, size(edPairCount), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    else
      call mpi_reduce(edPairCount, dummy_null, size(edPairCount), mpi_double_precision, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    end if
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end subroutine ed_collectCorr

  subroutine ed_average(numFrame, numMolType, frameCount)
    use utility, only: getMolTypePairIndexFromTypes
    implicit none
    integer, intent(in) :: numFrame, numMolType, frameCount(:)
    integer :: i, t1, t2, n, molTypePairIndex

    edPairCount = edPairCount / numFrame
    do n = 1, numMolType * (numMolType + 1) / 2
      do i = 1, num_eBin
        edCorr(:,i,n) = edCorr(:,i,n) / frameCount / edPairCount(i, n)
      end do
    end do

    do t2 = 1, numMolType
      do t1 = t2, numMolType
        if (t1 /= t2) then
          molTypePairIndex = getMolTypePairIndexFromTypes(t1, t2, numMolType)
          edPairCount(:, molTypePairIndex) = edPairCount(:, molTypePairIndex) / 2d0
        end if
      end do
    end do
  end subroutine ed_average

  subroutine ed_make_eBins()
    implicit none
    integer :: i, stat
    allocate(eBins(num_eBin), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: eBins"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if 
    eBins = [(i * eBinWidth, i = eBinIndex_absolute_min, eBinIndex_absolute_max)]
  end subroutine ed_make_eBins

  subroutine ed_finish()
    implicit none
    deallocate(eBinIndexAll, engLocLookupTable, edCorr)
    ! ed_binIndex has been deallocated in the main program
  end subroutine ed_finish
end module energy_dec
