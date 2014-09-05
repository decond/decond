module xdr
  use, intrinsic::iso_c_binding
  use utility, only : handle, f2c_string

  private 
  public :: open_trajectory, close_trajectory, read_trajectory, write_trajectory, get_natom

  integer, parameter :: line_len = 128, XDRFILE_MAX = 100
  type(C_PTR), dimension(XDRFILE_MAX) :: xdr_iohandle = c_null_ptr
  
  interface
    subroutine read_trr_natoms(filename, natoms) bind(C)
      use, intrinsic::iso_c_binding
      implicit none
      character(len=1, kind=C_CHAR), dimension(*), intent(in) :: filename
      integer(C_INT), intent(out) :: natoms
    end subroutine

    integer(C_INT) function read_trr(xd, natoms, step, t, lambda, box, x, v, f) bind(C)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: xd
      integer(C_INT), value :: natoms
      integer(C_INT) :: step
      real(C_DOUBLE) :: t, lambda
      real(C_DOUBLE) :: box(*), x(*), v(*)
      type(C_PTR), value :: f
    end function

    integer(C_INT) function write_trr(xd, natoms, step, t, lambda, box, x, v, f) bind(C)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), value :: xd
      integer(C_INT), value :: natoms
      integer(C_INT), value :: step
      real(C_DOUBLE), value :: t, lambda
      real(C_DOUBLE) :: box(*), x(*), v(*)
      type(C_PTR), value :: f
    end function

    type(C_PTR) function xdrfile_open(path, mode) bind(C)
      use, intrinsic :: iso_c_binding
      implicit none
      character(len=1, kind=C_CHAR), dimension(*) :: path, mode
    end function

    integer function xdrfile_close(xfp) bind(C)
      use, intrinsic :: iso_c_binding
      implicit none
      type(C_PTR), intent(in), value :: xfp
    end function
  end interface

contains
  integer function newunit(unit)
    use, intrinsic :: iso_c_binding
    implicit none
    integer, intent(out), optional :: unit
    integer, parameter :: LUN_MIN=1, LUN_MAX=XDRFILE_MAX
    integer :: lun
    ! begin
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
       if (.not. c_associated(xdr_iohandle(lun))) then
          newunit=lun
          exit
       end if
    end do
    if (present(unit)) unit=newunit
  end function
  
  ! Open trajectory and returns handle as htraj. 
  ! Should open fails, the program abends.
  subroutine open_trajectory(htraj, fname, mode)
    implicit none
    type(handle), intent(inout) :: htraj
    character(len=*), intent(in) :: fname
    character(len=1), optional, intent(in) :: mode
    character(len=1) :: fmode

    if (present(mode)) then
      fmode = mode
    else
      fmode = 'r'
    end if

    xdr_iohandle(newunit(htraj%iohandle)) = xdrfile_open(f2c_string(fname), fmode)
  end subroutine open_trajectory

  ! Close trajectory specified by handle
  subroutine close_trajectory(htraj)
    implicit none
    type(handle), intent(inout) :: htraj
    integer(C_INT) :: ret

    ret = xdrfile_close(xdr_iohandle(htraj%iohandle)) 
    xdr_iohandle(htraj%iohandle) = c_null_ptr
  end subroutine close_trajectory

  ! Read trajectory and returns [crd] as a coordinates, and [cell] as a
  ! periodic cell, represented in nanometer.
  ! [status] is non-zero if any error occurs. In such a case, [crd] and
  ! [cell] may be an arbitrary value if the trajectory does not contain
  ! cell information.
  ! The coordinate is not guaranteed to be within a unit cell.
  subroutine read_trajectory(htraj, natoms, is_periodic, crd, vel, cell, time, status)
    use, intrinsic :: iso_c_binding
    implicit none
    type(handle), intent(in) :: htraj
    integer(C_INT), intent(in) :: natoms
    logical, intent(in) :: is_periodic
    real(C_DOUBLE), intent(out) :: crd(3,natoms), vel(3,natoms)
    real(8), intent(out) :: cell(3)
    real(C_DOUBLE), intent(out) :: time
    integer, intent(out) :: status

    integer(C_INT) :: ret, step
    real(C_DOUBLE) :: lambda, box(3, 3)
    type(C_PTR), parameter :: f = C_NULL_PTR

    integer :: i

    ret = read_trr(xdr_iohandle(htraj%iohandle), natoms, step, time, lambda, box, crd, vel, f)

    do i = 1, 3
      cell(i) = box(i, i)
    end do

    status = ret
  end subroutine read_trajectory
  
  function get_natom(fname)
    implicit none
    integer :: get_natom
    character(len=*), intent(in):: fname

    get_natom = 0
    call read_trr_natoms(f2c_string(fname), get_natom)
  end function get_natom

  subroutine write_trajectory(htraj, natoms, is_periodic, crd, vel, cell, step, time, status)
    use, intrinsic :: iso_c_binding
    implicit none
    type(handle), intent(in) :: htraj
    integer(C_INT), intent(in) :: natoms
    logical, intent(in) :: is_periodic
    real(C_DOUBLE), intent(in) :: crd(3,natoms), vel(3,natoms)
    real(8), intent(in) :: cell(3)
    integer(C_INT) :: step
    real(C_DOUBLE), intent(in) :: time
    integer, intent(out) :: status

    integer(C_INT) :: ret
    real(C_DOUBLE) :: lambda, box(3, 3)
    type(C_PTR), parameter :: f = C_NULL_PTR

    integer :: i

    lambda = 0d0
    box = 0d0
    do i = 1, 3
      box(i, i) = cell(i)
    end do

    ret = write_trr(xdr_iohandle(htraj%iohandle), natoms, step, time, lambda, box, crd, vel, f)

    status = ret
  end subroutine write_trajectory

end module xdr

