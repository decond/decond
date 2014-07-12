module top
  use utility, only : handle, newunit, LINE_LEN, count_record_in_string
  implicit none

  private
  public :: open_top, close_top, read_top, system, print_sys
  
  character(len=*), parameter :: LINE_LEN_STR = "128"
  character(len=1), dimension(2), parameter :: COMMENT_CHAR = [";", "#"]
  
  character(len=LINE_LEN) :: current_directive = ''

  type atom
    character(len=LINE_LEN) :: type
    real(8) :: mass
    real(8) :: charge
  end type atom

  type molecule
    character(len=LINE_LEN) :: type
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

  subroutine read_top(htop, sys)
    implicit none
    type(handle), intent(in) :: htop
    type(system), intent(inout) :: sys
    type(system) :: sys_top
    integer :: i, j

    call read_top_system(htop, sys_top)
    call read_top_molecule(htop, sys_top%mol)

    do i = 1, size(sys%mol)
      do j = 1, size(sys_top%mol)
        if (sys%mol(i)%type == sys_top%mol(j)%type) then
          sys%mol(i) = sys_top%mol(j)
        end if
      end do
    end do
  end subroutine read_top

  subroutine read_top_system(htop, sys)
    implicit none
    type(handle), intent(in) :: htop
    type(system), intent(out) :: sys
    integer :: stat
    character(len=LINE_LEN) :: line
    integer :: num_moltype, i

    num_moltype = count_moltype(htop)
    allocate(sys%mol(num_moltype))

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
      read(line, *) sys%mol(i)%type, sys%mol(i)%num
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
        write(*,*) "Error: no atom records for ", adjustl(trim(mol(i)%type))
        call exit(1)
      end if
      allocate(mol(i)%atom(num_atom))

      ! read atomtypes from moleculetype
      call read_to_moltype(htop, mol(i), stat)
      call read_line(htop, line, status=stat) ! skip [ atom ] line
      do j = 1, num_atom
        call read_line(htop, line, status=stat)
        read(line, *) dum_i, mol(i)%atom(j)%type
      end do
    end do

    ! read data from [ atomtypes ] first
    do i = 1, size(mol)
      do j = 1, size(mol(i)%atom)
        call read_top_atype(htop, mol(i)%atom(j))
      end do
    end do

    ! read data from [ moleculetype ] -> [ atoms ]
    do i = 1, size(mol)
      call read_to_moltype(htop, mol(i), stat)
      call read_line(htop, line, status=stat) ! skip [ atoms ] line
      do j = 1, size(mol(i)%atom)
        call read_line(htop, line, status=stat)
        num_rec = count_record_in_string(line)
        select case (num_rec)
          case (7)
            read(line, *) dum_i, mol(i)%atom(j)%type, dum_i, dum_s, dum_s, dum_i, &
                          mol(i)%atom(j)%charge
          case (8)
            read(line, *) dum_i, mol(i)%atom(j)%type, dum_i, dum_s, dum_s, dum_i, &
                          mol(i)%atom(j)%charge, mol(i)%atom(j)%mass
          case default
            write(*,*) "Error: I cannot deal with num_rec=", num_rec, " for molecule ", mol(i)%type
            call exit(1)
        end select
      end do
    end do
    rewind(htop%iohandle)
  end subroutine read_top_molecule

  subroutine read_top_atype(htop, at)
    implicit none
    type(handle), intent(in) :: htop
    type(atom), intent(inout) :: at
    integer :: stat
    character(len=LINE_LEN) :: line, atype, dum_s
    integer :: i, j, num_atom, num_rec, dum_i
    real(8) :: dum_r, mass, charge
    logical :: is_read

    is_read = .false.
    rewind(htop%iohandle)
    do while(.true.)
      call read_to_next_directive(htop, "atomtypes", stat)
      if (stat < 0) exit  ! EOF
      do while(.true.)
        call read_line(htop, line, status=stat)
        if (stat < 0) then
          exit  ! EOF
        else if (stat == 1) then  ! new directive is met
          if (current_directive == "atomtypes") then
            cycle
          else
            exit
          end if
        end if
        num_rec = count_record_in_string(line)
        select case (num_rec)
          case (6)
            read(line, *, iostat=stat) atype, mass, charge, dum_s, dum_r, dum_r
          case (7)
            read(line, *, iostat=stat) atype, dum_i, mass, charge, dum_s, dum_r, dum_r
          case default
            write(*,*) "Error: I cannot deal with num_rec=", num_rec, " for atomtype ", trim(at%type)
            write(*,*) trim(line)
            call exit(1)
        end select
        if (stat /= 0) then
          write(*,*) "Error: something went wrong when reading atomtype ", trim(at%type)
          write(*,*) trim(line)
          call exit(1)
        end if 
        if (atype == at%type) then
          is_read = .true.
          at%mass = mass
          at%charge = charge
        end if
      end do
    end do

    if (.not. is_read) then
      write(*,*) "Error: no [ atomtypes ] for ", trim(at%type), " is read"
      call exit(1)
    end if
    rewind(htop%iohandle)
  end subroutine read_top_atype

  integer function count_molatom(htop, mol)
    implicit none
    type(handle), intent(in) :: htop
    type(molecule), intent(inout) :: mol
    integer :: stat
    character(len=LINE_LEN) :: line

    rewind(htop%iohandle)
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
      write(*,*) "Error: no [ atoms ] found after [ moleculetype ] ", trim(mol%type)
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
    
    rewind(htop%iohandle)
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
          if (name == mol%type) then
            return 
          end if
        end if
      end if
    end do
    rewind(htop%iohandle)
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
    integer :: idx, i

    do i = 1, size(COMMENT_CHAR)
      idx = index(line, COMMENT_CHAR(i))
      if (idx > 0) then
        line = line(:idx - 1)
      end if
    end do
  end subroutine remove_comment

  subroutine print_sys(sys)
    implicit none
    type(system) :: sys
    integer :: i, j
    write(*,*) "[ molecules ]"
    do i = 1, size(sys%mol)
      write(*, *) adjustl(trim(sys%mol(i)%type)), sys%mol(i)%num
    end do
    do i = 1, size(sys%mol)
      write(*,*) "[ moleculetype ]"
      write(*, *) adjustl(trim(sys%mol(i)%type))
      write(*,*) "[ atoms ]"
      write(*,*) "type  mass  charge"
      do j = 1, size(sys%mol(i)%atom)
        write(*,*) trim(sys%mol(i)%atom(j)%type), " ", sys%mol(i)%atom(j)%mass, sys%mol(i)%atom(j)%charge
      end do
    end do
  end subroutine print_sys
end module top

!program test
!  use utility, only : handle
!  use top, only: open_top, close_top, read_top, system
!  implicit none
!  integer, parameter :: LINE_LEN = 128
!  character(len=LINE_LEN) :: top_filename
!  type(handle) :: htop
!  type(system) :: sys
!
!  call get_command_argument(1, top_filename)
!  htop = open_top(top_filename)
!  call read_top(htop, sys)
!  call close_top(htop)
!end program test

