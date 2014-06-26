module top
  use utility, only : handle, newunit, LINE_LEN, count_record_in_string
  implicit none

  private
  public :: open_top, close_top, read_top, system
  
  character(len=*), parameter :: LINE_LEN_STR = "128"
  character(len=*), parameter :: COMMENT_CHAR= ";"
  
  character(len=LINE_LEN) :: current_directive = ''

  type atom
    character(len=LINE_LEN) :: type
    real(8) :: mass
    real(8) :: charge
  end type atom

  type molecule
    character(len=LINE_LEN) :: name
    integer :: num  ! number of molecules of the same type
    type(atom), dimension(:), allocatable :: atom
  end type molecule

  type system
    type(molecule), dimension(:), allocatable :: mol
  end type system
  
contains
  type(handle) function open_top(filename)
    implicit none
    character(len=*), intent(in) :: filename
  
    open_top%filename = filename
    open(unit=newunit(open_top%iohandle), file=filename, status='old', action="READ", form="FORMATTED")
  end function open_top

  subroutine close_top(htop)
    implicit none
    type(handle), intent(in) :: htop
    close(htop%iohandle)
  end subroutine close_top

  subroutine read_top(htop, tp)
    implicit none
    type(handle), intent(in) :: htop
    type(system), intent(out) :: tp
    integer :: i

    call read_top_system(htop, tp)
    call read_top_molecule(htop, tp%mol)

    write(*,*) "[ molecules ]"
    do i = 1, size(tp%mol)
      write(*, *) adjustl(trim(tp%mol(i)%name)), tp%mol(i)%num
    end do
  end subroutine read_top

  subroutine read_top_system(htop, tp)
    implicit none
    type(handle), intent(in) :: htop
    type(system), intent(out) :: tp
    integer :: stat
    character(len=LINE_LEN) :: line
    integer :: num_moltype, i

    num_moltype = count_moltype(htop)
    allocate(tp%mol(num_moltype))

    ! read name and number of molecules
    call read_to_next_directive(htop, 'molecules', status=stat)
    if (stat < 0) then
      write(*,*) "Error: no [ molecules ] in the top file"
      call exit(1)
    end if

    do i = 1, num_moltype
      call read_line(htop, line, status=stat)
      if (stat /= 0) then
        write(*,*) "Error reading name and number of molecules"
        call exit(1)
      end if
      read(line, *) tp%mol(i)%name, tp%mol(i)%num
    end do
    rewind(htop%iohandle)
  end subroutine read_top_system

  subroutine read_top_molecule(htop, mol)
    implicit none
    type(handle), intent(in) :: htop
    type(molecule), dimension(:), intent(inout) :: mol
    integer :: stat
    character(len=LINE_LEN) :: line, dum_s
    integer :: i, j, num_atom, num_rec, dum_i

    do i = 1, size(mol)
      num_atom = count_molatom(htop, mol(i))
      if (num_atom == 0) then
        write(*,*) "Error: no atom records for ", adjustl(trim(mol(i)%name))
        call exit(1)
      end if
      allocate(mol(i)%atom(num_atom))
      call read_to_moltype(htop, mol(i), stat)
      call read_line(htop, line, status=stat) ! skip [ atom ] line
write(*,*) "[ moleculetype ]"
write(*,*) trim(mol(i)%name)
write(*,*) "[ atom ]"
write(*,*) "type  mass  charge"
      do j = 1, num_atom
        call read_line(htop, line, status=stat)
        num_rec = count_record_in_string(line)
        select case (num_rec)
!          case (7)
!            read(line, *) dum_i, mol(i)%atom(j)%type, dum_i, dum_s, dum_s, dum_i, &
!                          mol(i)%atom(j)%charge
          case (8)
            read(line, *) dum_i, mol(i)%atom(j)%type, dum_i, dum_s, dum_s, dum_i, &
                          mol(i)%atom(j)%charge, mol(i)%atom(j)%mass
          case default
            write(*,*) "Error: I cannot deal with num_rec=", num_rec, " for molecule ", mol(i)%name
            call exit(1)
        end select
write(*,*) trim(mol(i)%atom(j)%type), " ", mol(i)%atom(j)%mass, mol(i)%atom(j)%charge
      end do
write(*,*)
    end do
    rewind(htop%iohandle)
  end subroutine read_top_molecule

  integer function count_molatom(htop, mol)
    implicit none
    type(handle), intent(in) :: htop
    type(molecule), intent(inout) :: mol
    integer :: stat
    character(len=LINE_LEN) :: line

    count_molatom = 0
    stat = 0

    call read_to_moltype(htop, mol, stat)
    if (stat < 0) then
      return  ! EOF
    end if
    
    call read_next_directive(htop, stat)
    if (stat < 0) then
      return  ! EOF
    else if (current_directive /= 'atoms') then
      write(*,*) "Error: no [ atoms ] found after [ moleculetype ] ", trim(mol%name)
      write(*,*) "Note: [ atoms ] must be the next directive of [ moleculetype ]"
      call exit(1)
    end if
    do while(.true.)
      call read_line(htop, line, status=stat)
      if (current_directive == 'atoms' .and. stat == 0) then
        count_molatom = count_molatom + 1
      else
        rewind(htop%iohandle)
        return
      end if
    end do
  end function count_molatom

  subroutine read_to_moltype(htop, mol, status)
    implicit none
    type(handle), intent(in) :: htop
    type(molecule), intent(in) :: mol
    integer, intent(out) :: status
    character(len=LINE_LEN) :: line, name
    
    do while(.true.)
      call read_to_next_directive(htop, 'moleculetype', status)
      if (status < 0) then
        return  ! EOF
      else if (status == 1) then
        call read_line(htop, line, status=status)
        if (status < 0) then
          return ! EOF
        else if (status == 0) then
          read(line, *) name
          if (name == mol%name) then
            return 
          end if
        end if
      end if
    end do
  end subroutine read_to_moltype
  
  subroutine read_to_next_directive(htop, directive, status)
    implicit none
    type(handle), intent(in) :: htop
    character(len=*), intent(in) :: directive
    integer, intent(out) :: status
    character(len=LINE_LEN) :: line
    do while(.true.)
      call read_line(htop, line, status=status)
      if (status < 0) then
        return ! EOF
      end if
      if (current_directive == directive .and. status == 1) then
        ! status == 1 means successfully read to a new directive
        return
      end if
    end do
  end subroutine read_to_next_directive

  subroutine read_next_directive(htop, status)
    implicit none
    type(handle), intent(in) :: htop
    integer, intent(out) :: status
    character(len=LINE_LEN) :: line
    do while(.true.)
      call read_line(htop, line, status=status)
      if (status < 0) then
        return ! EOF
      end if
      if (status == 1) then
        return  ! next directive returned as current_directive
      end if
    end do
  end subroutine read_next_directive

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

!program test
!  use utility, only : handle
!  use top, only: open_top, close_top, read_top, system
!  implicit none
!  integer, parameter :: LINE_LEN = 128
!  character(len=LINE_LEN) :: top_filename
!  type(handle) :: htop
!  type(system) :: tp
!
!  call get_command_argument(1, top_filename)
!  htop = open_top(top_filename)
!  call read_top(htop, tp)
!  call close_top(htop)
!end program test

