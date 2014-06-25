module top
  use utility, only : handle, newunit
  implicit none

  private
  public :: open_top, close_top, read_top, topology
  
  integer, parameter :: LINE_LEN = 128
  character(len=*), parameter :: LINE_LEN_STR = "128"
  character(len=*), parameter :: COMMENT_CHAR= ";"
  
  character(len=LINE_LEN) :: current_directive = ''

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
    integer :: num_moltype

    num_moltype = count_moltype(htop)
    write(*,*) "num_moltype=", num_moltype

    allocate(tp%mol(num_moltype))

    stat = 0
    do while(.true.)
      call read_line(htop, line, status=stat)
      if (stat < 0) then
        exit  ! EOF
      end if
    end do
  end subroutine read_top

  integer function count_moltype(htop)
    implicit none
    type(handle), intent(in) :: htop
    integer :: stat
    character(len=LINE_LEN) :: line

    count_moltype = 0
    stat = 0

    do while(.true.)
      call read_line(htop, line, status=stat)
      if (stat < 0) then
        exit  ! EOF
      end if

      if (current_directive == 'molecules' .and. stat == 0) then
        count_moltype = count_moltype + 1
      end if
    end do
    
    rewind(htop%iohandle)
  end function count_moltype

  subroutine read_line(htop, line, status)
    implicit none
    type(handle), intent(in) :: htop
    character(len=*), intent(out) :: line
    integer, intent(out) :: status
    integer :: stat
    character(len=LINE_LEN) :: directive
    character(len=1) :: dum_c
    
    status = 0
    do while(.true.)
      read(htop%iohandle, "(A"//LINE_LEN_STR//")", iostat=stat) line
      if (stat > 0) then
        write(*,*) "Error reading line"
        call exit(1)
      else if (stat < 0) then
        status = -1  ! EOF
        return
      end if

      line = adjustl(trim(line))
      call remove_comment(line)

      if (line /= '') then
        if (line(1:1) == '[') then
          ! set current directive
          read(line, *, iostat=stat) dum_c, directive, dum_c
          if (stat /= 0) then
            write(*,*) "Error reading [ directive ]"
            call exit(1)
          end if
          current_directive = directive
          status = 1  ! tell caller current_directive is changed at this line
        end if
        return
      end if  
    end do
  end subroutine read_line

  subroutine remove_comment(line)
    implicit none
    character(len=*), intent(inout) :: line
    integer :: idx

    idx = index(line, COMMENT_CHAR)
    if (idx > 0) then
      line = line(:idx - 1)
    end if
  end subroutine remove_comment
end module top

program test
  use utility, only : handle
  use top, only: open_top, close_top, read_top, topology
  implicit none
  integer, parameter :: LINE_LEN = 128
  character(len=LINE_LEN) :: top_filename
  type(handle) :: htop
  type(topology) :: tp
  call get_command_argument(1, top_filename)
  call open_top(htop, top_filename)
  call read_top(htop, tp)
  call close_top(htop)
end program test

