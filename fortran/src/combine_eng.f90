program combine_eng
  use hdf5
  implicit none
  integer, parameter :: FILENAME_LEN = 128
  character(len=*), parameter :: VERSION = "0.4.0"
  character(len=*), parameter :: OUTFILENAME = "engtraj-all.h5"
  character(len=*), parameter :: GROUP_ROOT = '/'
  integer(hsize_t), parameter :: VERSION_LEN = 11
  integer :: num_files, hdferr
  character(len=FILENAME_LEN), allocatable :: engfiles(:)
  integer, allocatable :: sltspec_list(:), sltfirsttag_list(:), numslt_list(:)
  integer, allocatable :: nummol_list(:), numpair_list(:), numframe_list(:)
  integer(hid_t) :: outfile_id, ver_type_id
  integer :: numpair_total

  call init()
  call check_version()
  call read_attrs()
  call sort_attrs()
  call check_sanity()
  call write_attrs()
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
      call get_command_argument(number=i, value=engfiles(i))
    end do

    allocate(sltspec_list(num_files))
    allocate(sltfirsttag_list(num_files))
    allocate(numslt_list(num_files))
    allocate(nummol_list(num_files))
    allocate(numpair_list(num_files))
    allocate(numframe_list(num_files))

    call h5open_f(hdferr)

    ! create output HDF5 file
    call h5fcreate_f(OUTFILENAME, H5F_ACC_EXCL_F, outfile_id, hdferr)
    if (hdferr /= 0) then
      write(*,*) "Failed to create HDF5 file: ", trim(adjustl(OUTFILENAME))
      write(*,*) "Probably the file already exists?"
      call exit(1)
    end if

    call h5tcreate_f(H5T_STRING_F, VERSION_LEN, ver_type_id, hdferr)
    call h5tset_strpad_f(ver_type_id, H5T_STR_NULLTERM_F, hdferr)
  end subroutine init

  subroutine check_version()
    use utility, only: parse_version
    implicit none
    integer, parameter :: MIN_ENGTRJ_VER_MAJOR = 0
    integer(hid_t) :: file_id, attr_id, type_id
    integer(hsize_t) :: dims(1)
    character(len=VERSION_LEN) :: buf
    integer :: i, major, minor, patch

    dims = -1  ! value of dims does not seem to have any effect
    do i = 1, num_files
      call h5fopen_f(engfiles(i), H5F_ACC_RDONLY_F, file_id, hdferr)

      call h5aopen_f(file_id, 'version', attr_id, hdferr)
      call h5aget_type_f(attr_id, type_id, hdferr)
      call h5aread_f(attr_id, type_id, buf, dims, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5fclose_f(file_id, hdferr)
      call parse_version(buf, major, minor, patch)
      if (major < MIN_ENGTRJ_VER_MAJOR) then
        write(*,*) "Only engtraj with a major version number greater than &
                   &or equal to ", MIN_ENGTRJ_VER_MAJOR, " is supported."
        write(*,*) "The major version of the file '", trim(engfiles(i)), &
                   " is: ", major
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
    call h5aopen_f(file_id, attr_name, attr_id, hdferr)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, hdferr)
    call h5aclose_f(attr_id, hdferr)
  end subroutine read_attr_integer

  subroutine read_attrs()
    implicit none
    integer :: i
    integer(hid_t) :: file_id, attr_id
    integer(hsize_t) :: dims(1)
    integer :: buf

    do i = 1, num_files
      call h5fopen_f(engfiles(i), H5F_ACC_RDONLY_F, file_id, hdferr)

      call read_attr_integer(file_id, 'sltspec', buf)
      sltspec_list(i) = buf

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

      call h5fclose_f(file_id, hdferr)
    end do
  end subroutine read_attrs

  subroutine sort_attrs()
    ! sort attributes according to sltspec
    implicit none
    integer :: i, j

    do i = 1, num_files
      do j = i, num_files
        if (sltspec_list(j) == i) then
          if (j /= i) then
            call swap_int(sltspec_list, i, j)
            call swap_int(sltfirsttag_list, i, j)
            call swap_int(numslt_list, i, j)
            call swap_int(nummol_list, i, j)
            call swap_int(numpair_list, i, j)
            call swap_int(numframe_list, i, j)
            call swap_char(engfiles, i, j)
          end if
          exit
        end if
      end do
    end do
  end subroutine sort_attrs

  subroutine swap_int(list, i, j)
    implicit none
    integer, intent(inout) :: list(:)
    integer, intent(in) :: i, j
    integer :: tmp
    tmp = list(i)
    list(i) = list(j)
    list(j) = tmp
  end subroutine swap_int

  subroutine swap_char(list, i, j)
    implicit none
    character(len=FILENAME_LEN), intent(inout) :: list(:)
    integer, intent(in) :: i, j
    character(len=FILENAME_LEN) :: tmp
    tmp = list(i)
    list(i) = list(j)
    list(j) = tmp
  end subroutine swap_char

  subroutine check_sanity()
    implicit none
    integer :: i

    do i = 1, num_files
      if (sltspec_list(i) /= i) then
        write(*,*) "sltspec_list is not sorted correctly. it should consist of &
                   &continuous integers starting from 1"
        write(*,*) "sltspec_list =", sltspec_list
        call exit(1)
      end if
    end do

    do i = 1, num_files - 1
      if (sltfirsttag_list(i) >= sltfirsttag_list(i+1)) then
        write(*,*) "sltfirsttag_list should be an ascending list"
        write(*,*) "sltfirsttag_list =", sltfirsttag_list
        call exit(1)
      end if
    end do

    if (.not. all(nummol_list == nummol_list(1))) then
      write(*,*) "nummol's should be all the same"
      write(*,*) "nummol_list =", nummol_list
      call exit(1)
    end if

    if (.not. all(numframe_list == numframe_list(1))) then
      write(*,*) "numframe's should be all the same"
      write(*,*) "numframe_list =", numframe_list
      call exit(1)
    end if

    if (sum(numslt_list) /= nummol_list(1)) then
      write(*,*) "sum of numslt's should equal to nummol"
      write(*,*) "numslt_list =", numslt_list
      call exit(1)
    end if

    if (sum(numpair_list) /= nummol_list(1) * (nummol_list(1) - 1) / 2) then
      write(*,*) "numpair is not consistent with nummol"
      write(*,*) "numpair_list =", numpair_list
      call exit(1)
    end if
  end subroutine check_sanity

  subroutine write_attr_version()
    implicit none
    integer(hid_t) :: space_id, attr_id
    integer(hsize_t) :: dims(1)

    dims = -1  ! value of dims does not affect scalar dataspace
    call h5screate_f(H5S_SCALAR_F, space_id, hdferr)
    call h5acreate_f(outfile_id, 'version', ver_type_id, space_id, attr_id, hdferr)
    call h5awrite_f(attr_id, ver_type_id, VERSION, dims, hdferr)
    call h5sclose_f(space_id, hdferr)
    call h5aclose_f(attr_id, hdferr)
  end subroutine write_attr_version

  subroutine write_attr_integers(attr_name, buf)
    implicit none
    character(len=*), intent(in) :: attr_name
    integer, intent(in) :: buf(:)
    integer(hid_t) :: file_id, space_id, attr_id
    integer(hsize_t) :: buf_size
    integer(hsize_t) :: dims(1)

    dims = -1  ! value of dims does not seem to have any effect
    buf_size = size(buf)
    call h5screate_simple_f(1, [buf_size], space_id, hdferr)
    call h5acreate_f(outfile_id, attr_name, H5T_NATIVE_INTEGER, space_id, attr_id, hdferr)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, buf, dims, hdferr)
    call h5sclose_f(space_id, hdferr)
    call h5aclose_f(attr_id, hdferr)
  end subroutine write_attr_integers

  subroutine write_attrs()
    implicit none

    numpair_total = sum(numpair_list)
    call write_attr_version()
    call write_attr_integers('sltspec', sltspec_list)
    call write_attr_integers('sltfirsttag', sltfirsttag_list)
    call write_attr_integers('numslt', numslt_list)
    call write_attr_integers('nummol', [nummol_list(1)])
    call write_attr_integers('numpair', [numpair_total])
    call write_attr_integers('numframe', [numframe_list(1)])
  end subroutine write_attrs

  subroutine combine_engtraj()
    implicit none
    character(len=*), parameter :: ENGDSET_NAME = "energy"
    integer(hid_t) :: infile_id, indset_id, outdset_id, memspace_id, outspace_id
    real(8), allocatable :: energy(:, :)
    integer(hsize_t) :: idx_offset, offset(2), block(2), stride(2), count(2)
    integer(hsize_t) :: out_energy_dims(2)
    integer :: i

    out_energy_dims = [numframe_list(1), numpair_total]
    call h5screate_simple_f(2, out_energy_dims, outspace_id, hdferr)
    call h5dcreate_f(outfile_id, ENGDSET_NAME, H5T_NATIVE_DOUBLE, outspace_id, outdset_id, hdferr)

    idx_offset = 0
    do i = 1, num_files
      allocate(energy(numframe_list(1), numpair_list(i)))
      call h5screate_simple_f(2, shape(energy, kind=hsize_t), memspace_id, hdferr)

      call h5fopen_f(engfiles(i), H5F_ACC_RDONLY_F, infile_id, hdferr)
      call h5dopen_f(infile_id, ENGDSET_NAME, indset_id, hdferr)
      call h5dread_f(indset_id, H5T_NATIVE_DOUBLE, energy, &
                     shape(energy, kind=hsize_t), hdferr)

      offset = [0_hsize_t, idx_offset]
      block = [numframe_list(1), numpair_list(i)]
      stride = [1, 1]
      count = [1, 1]
      ! select data slab in filespace
      call h5sselect_hyperslab_f(outspace_id, H5S_SELECT_SET_F, offset, count, hdferr, stride, block)
      call h5dwrite_f(outdset_id, H5T_NATIVE_DOUBLE, energy, shape(energy, kind=hsize_t), hdferr, &
                      file_space_id=outspace_id, mem_space_id=memspace_id)

      call h5dclose_f(indset_id, hdferr)
      call h5fclose_f(infile_id, hdferr)
      call h5sclose_f(memspace_id, hdferr)
      deallocate(energy)

      idx_offset = idx_offset + numpair_list(i)
    end do

    call h5dclose_f(outdset_id, hdferr)
    call h5sclose_f(outspace_id, hdferr)
  end subroutine combine_engtraj

  subroutine finish()
    implicit none
    call h5fclose_f(outfile_id, hdferr)
    call h5close_f(hdferr)
    write(*,*) "Output ", OUTFILENAME
  end subroutine finish
end program combine_eng
