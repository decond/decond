module xyz
  use varpars, only: dp
  implicit none
  private
  integer, parameter :: line_len = 1024
  character(len=*), parameter :: format_begin = "(A, 3F", &
                                 format_end = ")", format_r = "F"
  public :: open_xyz, write_xyz, read_xyz, read_natom_xyz, close_xyz

contains
  subroutine open_xyz(unit, file, status)
    implicit none
    integer, intent(out) :: unit
    character(len=*), intent(in) :: file, status

    open(newunit=unit, file=file, status=status)
  end subroutine

  subroutine write_xyz(unit, coord, num_atom, info, atom_name, opt_data)
    implicit none
    integer, intent(in) :: unit
    real(dp), intent(in) :: coord(:, :)  !(3, num_atom)
    integer, intent(in), optional :: num_atom
    character(len=*), intent(in), optional :: info, atom_name(:)
    real(dp), intent(in), optional :: opt_data(:, :)  ! (:, num_atom)
    integer :: num_atom_
    character(len=line_len) :: info_, atom_name_(size(coord, 2))
    character(len=line_len) :: format_opt, num_opt_c, format
    integer :: i

    num_atom_ = size(coord, 2)
    if (present(num_atom)) then
      if (num_atom /= num_atom_) then
        stop "num_atom is not consistent with coord"
      end if
    end if

    if (present(info)) then
      info_ = info
    else
      info_ = "XYZ"
    end if

    if (present(atom_name)) then
      atom_name_ = atom_name
    else
      atom_name_ = "X"
    end if

    write(unit, *) num_atom_
    write(unit, "(A)") trim(info_)

    if (present(opt_data)) then
      write(num_opt_c, *) size(opt_data, 1)
      format_opt = ", " // trim(adjustl(num_opt_c)) // format_r
      format = format_begin // trim(format_opt) // format_end
      do i = 1, num_atom_
        write(unit, format) trim(atom_name_(i)), coord(:, i), opt_data(:, i)
      end do
    else
      format = format_begin // format_end
      do i = 1, num_atom_
        write(unit, format) trim(atom_name_(i)), coord(:, i)
      end do
    end if
  end subroutine write_xyz

  subroutine read_xyz(unit, coord, num_atom, info, atom_name, opt_data)
    implicit none
    integer, intent(in) :: unit
    real(dp), intent(out) :: coord(:, :)  !(3, num_atom)
    integer, intent(out), optional :: num_atom
    character(len=*), intent(out), optional :: info, atom_name(:)
    real(dp), intent(out), optional :: opt_data(:, :)  ! (:, num_atom)
    integer :: num_atom_
    character(len=line_len) :: info_, atom_name_(size(coord, 2))
    integer :: i

    read(unit, *) num_atom_
    if (present(num_atom)) num_atom = num_atom_

    read(unit, "(A)") info_
    if (present(info)) info = info_

    if (present(opt_data)) then
      do i = 1, num_atom_
        read(unit, *) atom_name_(i), coord(:, i), opt_data(:, i)
      end do
    else
      do i = 1, num_atom_
        read(unit, *) atom_name_(i), coord(:, i)
      end do
    end if

    if (present(atom_name)) atom_name = atom_name_
  end subroutine read_xyz

  integer function read_natom_xyz(file) result(natom)
    implicit none
    character(len=*), intent(in) :: file
    integer :: unit

    call open_xyz(unit, file, 'old')
    read(unit, *) natom
    call close_xyz(unit)
  end function read_natom_xyz

  subroutine close_xyz(unit)
    implicit none
    integer, intent(in) :: unit

    close(unit)
  end subroutine close_xyz
end module xyz
