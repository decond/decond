! -*- F90 -*-

! Writing octave compatible save format 
module octave_save
  character(len=*), parameter :: f_fmt = "EN30.16"

  type handle
     integer :: iohandle
  end type handle
  
contains
  ! Open file and returns handle as htraj. 
  ! Should open fails, the program abends.
  subroutine create_octave(htraj, fname)
    use utility
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname

    ! use form="FORMATTED" for human-redable file format
    open(unit=newunit(htraj%iohandle), file=fname, action="WRITE", form="FORMATTED")
  end subroutine create_octave

  ! Close trajectory specified by handle
  subroutine close_octave(htraj)
    implicit none
    type(handle), intent(inout) :: htraj

    close(htraj%iohandle)

    ! release extra resources if necessary
  end subroutine close_octave

  subroutine write_octave_scalar(htraj, name, data_s)
    implicit none
    type(handle), intent(in) :: htraj
    character(len=*) :: name
    real(8) :: data_s
    
    write(htraj%iohandle, *) "# name: " // trim(adjustl(name))
    write(htraj%iohandle, *) "# type: scalar" 
    write(htraj%iohandle, "("//f_fmt//")") data_s
    write(htraj%iohandle, *) 
    write(htraj%iohandle, *) 
  end subroutine write_octave_scalar

  subroutine write_octave_mat2(htraj, name, data_m)
    implicit none
    type(handle), intent(in) :: htraj
    character(len=*), intent(in) :: name
    real(8), intent(in) :: data_m(:,:)
    integer :: i, j
    
    write(htraj%iohandle, *) "# name: " // trim(adjustl(name))
    write(htraj%iohandle, *) "# type: matrix" 
    write(htraj%iohandle, *) "# rows: ", size(data_m, 1)
    write(htraj%iohandle, *) "# columns: ", size(data_m, 2)
    do i = 1, size(data_m, 1)
      do j = 1, size(data_m, 2)
        write(htraj%iohandle, "("//f_fmt//",1X)", advance="no") data_m(i,j)
      end do
      write(htraj%iohandle, *) 
    end do
    write(htraj%iohandle, *) 
    write(htraj%iohandle, *) 
  end subroutine write_octave_mat2

  subroutine write_octave_vec(htraj, name, data_v)
    implicit none
    type(handle), intent(in) :: htraj
    character(len=*), intent(in) :: name
    real(8), intent(in) :: data_v(:)
    real(8)  :: data_m(size(data_v), 1)
    
    data_m(:, 1) = data_v
    call write_octave_mat2(htraj, name, data_m)
  end subroutine write_octave_vec

  subroutine write_octave_mat3(htraj, name, data_m)
    implicit none
    type(handle), intent(in) :: htraj
    character(len=*), intent(in) :: name
    real(8) :: data_m(:,:,:)
    integer :: i, j, k

    write(htraj%iohandle, *) "# name: " // trim(adjustl(name))
    write(htraj%iohandle, *) "# type: matrix" 
    write(htraj%iohandle, *) "# ndims: 3"
    write(htraj%iohandle, *) size(data_m, 1), size(data_m, 2), size(data_m, 3)
    do k = 1, size(data_m, 3)
      do j = 1, size(data_m, 2)
        do i = 1, size(data_m, 1)
          write(htraj%iohandle, "("//f_fmt//")") data_m(i,j,k)
        end do
      end do
    end do
    write(htraj%iohandle, *) 
    write(htraj%iohandle, *) 
  end subroutine write_octave_mat3

end module octave_save

