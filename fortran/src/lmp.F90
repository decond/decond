module lmp
  use varpars, only: rk
#ifndef FORTRAN2008
  use utility, only: newunit
#endif

  implicit none
  private
  public :: open_lmp, read_lmp, read_natom_lmp, close_lmp

contains
  subroutine open_lmp(unit, file, mode)
    implicit none
    integer, intent(out) :: unit
    character(len=*), intent(in) :: file
    character(len=1), optional, intent(in) :: mode
    character(len=1) :: fmode
    character(len=7) :: fstatus

    if (present(mode)) then
      fmode = mode
    else
      fmode = 'r'
    end if

    select case (fmode)
    case ('r')
      fstatus = 'old'
    case ('w')
      fstatus = 'new'
    end select

#ifdef FORTRAN2008
    open(newunit=unit, file=file, status=fstatus)
#else
    open(newunit(unit), file=file, status=fstatus)
#endif
  end subroutine

  subroutine read_lmp(unit, coord, num_atom, step, cell, opt_data)
    use varpars, only: world_dim
    implicit none
    integer, intent(in) :: unit
    real(rk), intent(out) :: coord(:, :)  !(3, num_atom)
    integer, intent(out), optional :: num_atom
    integer, intent(out), optional :: step
    real(rk), intent(out), optional :: cell(world_dim)
    real(rk), intent(out), optional :: opt_data(:, :)  ! (:, num_atom)
    integer :: num_atom_, step_
    real(rk) :: cell_(world_dim)
    integer :: i, dum_i

    read(unit, *)
    read(unit, *) step_
    read(unit, *)
    read(unit, *) num_atom_
    read(unit, *)
    do i = 1, world_dim
      read(unit, *) dum_i, cell_(i)
    end do
    if (present(num_atom)) num_atom = num_atom_
    if (present(step)) step = step_
    if (present(cell)) cell = cell_

    read(unit, *)
    if (present(opt_data)) then
      do i = 1, num_atom_
        read(unit, *) dum_i, coord(:, i), opt_data(:, i)
      end do
    else
      do i = 1, num_atom_
        read(unit, *) dum_i, coord(:, i)
      end do
    end if
  end subroutine read_lmp

  integer function read_natom_lmp(file) result(natom)
    implicit none
    character(len=*), intent(in) :: file
    integer :: unit

    call open_lmp(unit, file, 'r')
    read(unit, *)
    read(unit, *)
    read(unit, *)
    read(unit, *) natom
    call close_lmp(unit)
  end function read_natom_lmp

  subroutine close_lmp(unit)
    implicit none
    integer, intent(in) :: unit

    close(unit)
  end subroutine close_lmp
end module lmp
