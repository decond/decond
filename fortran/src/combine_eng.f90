program combine_eng
  use hdf5
  implicit none
  integer, parameter :: FILENAME_LEN = 128
  character(len=*), parameter :: VERSION = "0.1.0"
  character(len=*), parameter :: OUTFILENAME = "engtraj-all.h5"
  integer :: num_files, hdferr
  character(len=FILENAME_LEN), allocatable :: engfiles(:)

  call init()
  call read_attr()
  call sort_attr()
  call check_sanity()
  call combine_engtraj()
  call finish()

contains
  subroutine print_usage()
    implicit none
    write(*, *) "usage: $ combine_eng <engfile1> <engfile2> ..."
  end subroutine print_usage

  subroutine init()
    implicit none
    integer :: i
    num_files = command_argument_count()
    if (num_files < 2) then
      write(*, *) "at least two engfiles must be given"
      call print_usage()
      call exit(1)
    end if
    
    allocate(engfiles(num_files))
    do i = 1, num_files
      call get_command_argument(number=i, value=engfile(i))
    end do

    allocate(sltspecList(num_files))
    allocate(sltfirsttagList(num_files))
    allocate(numsltList(num_files))
    allocate(nummolList(num_files))
    allocate(numpairList(num_files))
    allocate(numframeList(num_files))
    call h5open_f(hdferr)
  end subroutine init

  subroutine read_attr()
    implicit none
    integer :: i
    integer(hid_t) :: file_id, dset_id, attr_id
    integer(hsize_t) :: dims
    integer :: buf
    character(len=*), parameter :: GROUP_ROOT = '/'

    dims = 0  ! dims is ignored since buf is a scalar
    do i = 1, num_files
      call h5fopen_f(engfiles[i], H5F_ACC_RDONLY_F, file_id, hdferr)
      call h5dopen_f(file_id, GROUP_ROOT, dset_id, hdferr)

      call h5aopen_f(dset_id, 'sltspec', attr_id, hdferr)
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, hdferr)
      sltspecList[i] = buf
      call h5aclose_f(attr_id, hdferr)

      call h5aopen_f(dset_id, 'sltfirsttag', attr_id, hdferr)
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, hdferr)
      sltfirsttagList[i] = buf
      call h5aclose_f(attr_id, hdferr)

      call h5aopen_f(dset_id, 'numslt', attr_id, hdferr)
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, hdferr)
      numsltList[i] = buf
      call h5aclose_f(attr_id, hdferr)

      call h5aopen_f(dset_id, 'nummol', attr_id, hdferr)
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, hdferr)
      nummolList[i] = buf
      call h5aclose_f(attr_id, hdferr)

      call h5aopen_f(dset_id, 'numpair', attr_id, hdferr)
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, hdferr)
      numpairList[i] = buf
      call h5aclose_f(attr_id, hdferr)

      call h5aopen_f(dset_id, 'numframe', attr_id, hdferr)
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, hdferr)
      numframeList[i] = buf
      call h5aclose_f(attr_id, hdferr)

      call h5dclose_f(dset_id, hdferr)
      call h5fclose_f(file_id, hdferr)
    end do
  end subroutine read_attr

  subroutine sort_attr()
    ! sort attributes according to sltspec
    implicit none
    integer :: i, j

    do i = 1, num_files
      do j = i, num_files
        if (sltspecList[j] == i) then
          if (j /= i) then
            call swap_int(sltspecList, i, j)
            call swap_int(sltfirsttagList, i, j)
            call swap_int(numsltList, i, j)
            call swap_int(nummolList, i, j)
            call swap_int(numpairList, i, j)
            call swap_int(numframeList, i, j)
            call swap_char(engfiles, i, j)
          end if
          exit
        end if
      end do
    end do
  end subroutine sort_attr

  subroutine swap_int(list, i, j)
    ! TODO
    implicit none
    integer, intent(inout) :: list(:)
    integer, intent(in) :: i, j
  end subroutine swap_int

  subroutine swap_char(list, i, j)
    ! TODO
    implicit none
    character(len=FILENAME_LEN), intent(inout) :: list(:)
    integer, intent(in) :: i, j
  end subroutine swap_char

  subroutine check_version()
    ! TODO
    use utility, only: parse_version
    implicit none
  end subroutine check_version

  subroutine check_sanity()
    ! TODO
    implicit none
    call check_version()
  end subroutine check_sanity

  subroutine combine_engtraj()
    ! TODO
    implicit none
    integer :: outfileid
    ! create output HDF5 file
    call h5fcreate_f(OUTFILENAME, H5F_ACC_EXCL_F, outfileid, hdferr)
    if (ierr /= 0) then
      write(*,*) "Failed to create HDF5 file: ", trim(adjustl(outCorrFilename))
      write(*,*) "Probably the file already exists?"
      call mpi_abort(MPI_COMM_WORLD, 1, ierr);
      call exit(1)
    end if
    call h5fclose_f(outfileid, hdferr)
  end subroutine combine_engtraj

  subroutine finish()
    implicit none
    call h5close_f(hdferr)
  end subroutine finish
end program combine_eng
