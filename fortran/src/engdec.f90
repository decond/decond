module engdec
  use mpiproc
  implicit none
  integer :: num_eBin
  character(len=128) :: engtrjFilename
  real(8), allocatable :: eBinIndex(:, :)  !eBinIndex(numFrame, ?)
  integer, allocatable :: storageIndexTable(:, :)
  integer, parameter :: MIN_ENGTRJ_VER_MAJOR = 1

contains
  subroutine ed_init()
    implicit none
    num_eBin = 500  ! default value
  end subroutine sd_init

  subroutine ed_readEng(engtrjFilename, numFrame, skip)
    use HDF5
    implicit none
    character(len=*), intent(in) :: engtrjFilename
    integer, intent(in) :: numFrame, skip
    integer(hid_t) :: engtrjFileid
    integer :: r, c, loc, lastLoc, locMax
    real(8), allocatable :: eng(:)

    call openEngtrj(engtrjFileid)
    call makeStorageIndexTable(r_start, r_end, c_start, c_end, locMax)
    allocate(eBinIndex(numFrame, locMax), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: storageIndexTable"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    allocate(eng(numFrame), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: eng"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    lastLoc = 0
    do c = c_start, c_end
      do r = r_start, r_end
        if (r /= c) then
          loc = storageIndexTable(r, c)
          if (loc > lastLoc) then
            ! new loc, read new data
            call readPairEng(eng, r, c, engtrjFileid, numFrame, skip)
            call eng2BinIndex(eng, eBinIndex_single)
            eBinIndex(:, loc) = eBinIndex_single
            lastLoc = loc
          end if
        end if
      end do
    end do

    call H5Fclose_f(engtrjFileid, ierr)
  end subroutine ed_readEng

  subroutine openEngtrj(engtrjFileid)
    use HDF5
    implicit none
    integer(hid_t), intent(out) :: engtrjFileid
    character(len=11) :: engtrjVer
    integer :: engtrjVer_major
    integer(hid_t) :: plist_id  ! property list identifier for HDF5

    ! initialize HDF5 Fortran predefined datatypes
    call H5open_f(ierr)

    ! setup file access property list with parallel I/O access.
    call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierr)
    call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, ierr)

    ! open the existing engtrj file.
    call H5Fopen_f(engtrjFilename, H5F_ACC_RDONLY, engtrjFileid, ierr, access_prp = plist_id)
    if (ierr /= 0) then
      write(*,*) "Failed to open HDF5 file: ", engtrjFilename
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    ! close property list
    call H5Pclose_f(plist_id, ierr)

    ! check compatibility at root
    if (myrank == root) then
      ! read version
      H5LTget_attribute_string_f(engtrjFileid, "/", "version", engtrjVer)
      call parseVersion(engtrjVer, engtrjVer_major)
      if (engtrjVer_major < MIN_ENGTRJ_VER_MAJOR) then
        write(*,*) "Only engtrj with a major version number greater than ", MIN_ENGTRJ_VER_MAJOR, " is supported."
        write(*,*) "The major version of the file '", trim(engtrjFilename), " is: ", engtrjVer_major
        call mpi_abort(MPI_COMM_WORLD, 1, ierr);
        call exit(1)
      end if
    end if
  end subroutine openEngtrj

  subroutine parseVersion(ver, major, minor, patch)
    implicit none
    character(len=11) :: ver
    integer, intent(out) :: major
    integer, optional, intent(out) :: minor, patch
    p1 = scan(ver, '.')
    p2 = scan(ver, '.', .true.)
    read(ver(1:p1-1), *) major
    if (present(minor)) read(ver(p1+1:p2-1), *) minor
    if (present(patch)) read(ver(p2+1:), *) patch
  end subroutine parseVersion

  subroutine readPairEng(eng, r, c, engtrjFileid, numFrame, skip)
    use HDF5
    use H5LT
    implicit none
    real(8), intent(out) :: eng(numFrame)
    integer, intent(in) :: r, c, numFrame, skip
    integer(hid_t), intent(in) :: engtrjFileid
    integer :: numTotMol, pairIndex

    H5LTget_attribute_int_f(engtrjFileid, "/", "numTotalMolecule", numTotMol)
    pairIndex = getPairIndex(r, c, numTotMol)

    ! select the desired data slab
  end subroutine readPairEng

  integer function getPairIndex(r, c, n)
    implicit none
    integer, intent(in) :: r, c, n

    getPairIndex = (r - 1) * n + c - ((r + 1) * r) / 2
  end function getPairIndex

  ! index = 0 : self-interaction energy (ignore this record)
  ! index = N : N is the location of the stored data (in the eBinIndex(:, N))
  subroutine makeStorageIndexTable(r_start, r_end, c_start, c_end, locMax)
    implicit none
    integer, intent(in) :: r_start, r_end, c_start, c_end
    integer, intent(out) :: locMax
    integer :: loc, pairHash, cacheIndex
    integer :: pairHashCache((r_end - r_start + 1) * (c_end - c_start + 1))

    allocate(storageIndexTable(r_start:r_end, c_start:c_end), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: storageIndexTable"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if

    pairHashCache = -1
    storageIndexTable = 0
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

            ! make a storageIndexTable so that we can know
            ! where the data for a certain pair are stored later
            pairHashCache(cacheIndex) = pairHash
            storageIndexTable(r, c) = loc
          else
            ! this is an old pair, set the index to the previous symmetric one
            storageIndexTable(r, c) = storageIndexTable(c, r)
          end if
        end if
      end do
    end do

    locMax = loc
  end subroutine makeStorageIndexTable

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
end module engdec
