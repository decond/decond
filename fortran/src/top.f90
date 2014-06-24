module top
  use utility, only : handle, newunit
  implicit none
  
  integer, parameter :: LINE_LEN = 128
  character(len=*), parameter :: LINE_LEN_STR = "128"

  type atom
    real(8) :: mass
    real(8) :: charge
  end type atom

  type molecule
    character(len=LINE_LEN) :: name
    integer :: num  ! number of molecules of the same type
    integer :: num_atom  ! number of atoms per molecule
    type(atom), dimension(:), allocatable :: atom
  end type molecule

  type topology
    integer :: num_moltype
    type(molecule), dimension(:), allocatable :: mol
  end type topology
  

contains
  subroutine open_top(htop, filename)
    implicit none
    type(handle), intent(inout) :: htop
    character(len=*), intent(in) :: filename
    open(unit=newunit(htop%iohandle), file=filename, status='old', action="READ", form="FORMATTED")
  end subroutine open_top

  subroutine close_top(htop)
    implicit none
    type(handle), intent(in) :: htop
    close(htop%iohandle)
  end subroutine close_top

  subroutine read_top(htop, tp)
    implicit none
    type(handle), intent(in) :: htop
    type(topology), intent(out) :: tp
    integer :: stat
    character(len=LINE_LEN) :: line

    stat = 0
    do while(stat == 0)
      read(htop%iohandle, "(A"//LINE_LEN_STR//")", iostat=stat) line
      line = adjustl(trim(line))
      if (line(1:1) == ';') then
        cycle
      end if
      write(*,*) adjustl(trim(line))
    end do
  end subroutine read_top

!  integer function get_top_num_moltype(htop)
!    implicit none
!    type(handle), intent(in) :: htop
!    integer :: stat
!    character(len=LINE_LEN) :: line
!
!    stat = 0
!    do while(stat >= 0)
!      read(htop%iohandle, "()", iostat=stat)
!    end do
!    get_top_num_moltype
!  end subroutine read_top_molecules
!
!  subroutine read_top_molecules(htop, cmp_name, num)
!    implicit none
!    type(handle), intent(in) :: htop
!    character(len=LINE_LEN), dimension(:), allocatable :: 
!  end subroutine read_top_molecules

end module top

program test
  use top
  implicit none
  character(len=LINE_LEN) :: top_filename
  type(handle) :: htop
  type(topology) :: tp
  call get_command_argument(1, top_filename)
  call open_top(htop, top_filename)
  call read_top(htop, tp)
  call close_top(htop)
end program test

