module spatial_dec
  use mpiproc
  use varpars, only: rk, numframe, totnummol, world_dim, &
                     do_diagonal, do_orthogonal
  implicit none
  private
  public sd_init, com_pos, sd_prep_corrmemory, sd_getbinindex, &
         sd_cal_num_rbin, sd_broadcastpos, sd_prep, &
         sd_collectcorr, sd_average, sd_make_rbins, sd_finish
  real(rk), public :: rbinwidth
  !sdcorr: spatially decomposed correlation (lag, rBin, moltypepair_idx)
  !sdpaircount: (num_rbin, moltypepair_idx)
  real(rk), public, allocatable :: sdpaircount(:, :), sdcorr(:, :, :), pos(:, :, :)
  integer, public :: num_rbin
  integer, public, allocatable :: sd_binIndex(:)
  real(rk), public, allocatable :: rbins(:)
  !MPI variables
  real(rk), public, allocatable :: pos_r(:, :, :), pos_c(:, :, :)
  real(rk), public :: od_tol

  integer :: stat

contains
  subroutine sd_init()
    rbinwidth = 0.01
  end subroutine sd_init

  subroutine sd_prep()
    allocate(pos_r(world_dim, numframe, num_r), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_r"
      call mpi_abend()
    end if

    allocate(pos_c(world_dim, numframe, num_c), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: pos_c"
      call mpi_abend()
    end if

    if (myrank == root) then
      allocate(pos(world_dim, numframe, totnummol), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: pos"
        call mpi_abend()
      end if
    else
      !not root, allocate dummy pos to inhibit error messages
      allocate(pos(1, 1, 1), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: dummy pos on rank", myrank
        call mpi_abend()
      end if
    end if
  end subroutine sd_prep

  subroutine com_pos(com_p, pos, start_index, sys, cell)
    use top, only : system
    implicit none
    real(rk), dimension(:, :), intent(out) :: com_p
    real(rk), dimension(:, :), intent(in) :: pos
    real(rk), dimension(:, :), allocatable :: pos_gathered
    integer, dimension(:), intent(in) :: start_index
    type(system), intent(in) :: sys
    real(rk), dimension(world_dim), intent(in) :: cell
    integer :: d, i, j, idx_begin, idx_end, idx_com, num_atom

    idx_com = 0
    do i = 1, size(sys%mol)
      num_atom = size(sys%mol(i)%atom)
      allocate(pos_gathered(world_dim, num_atom), stat=stat)
      if (stat /=0) then
        write(*,*) "Allocation failed: pos_gathered"
        call mpi_abend()
      end if
      do j = 1, sys%mol(i)%num
        idx_begin = start_index(i) + (j-1) * num_atom
        idx_end = idx_begin + num_atom - 1
        idx_com = idx_com + 1
        call gatherMolPos(pos_gathered, pos(:, idx_begin:idx_end), cell)
        do d = 1, world_dim
          com_p(d, idx_com) = sum(pos_gathered(d, :) * sys%mol(i)%atom(:)%mass) / sum(sys%mol(i)%atom(:)%mass)
        end do
      end do
      deallocate(pos_gathered)
    end do
  end subroutine com_pos

  ! make sure the gathered positions are from the same moleucle, instead of different images
  subroutine gatherMolPos(pos_gathered, pos, cell)
    implicit none
    real(rk), dimension(:, :), intent(out) :: pos_gathered
    real(rk), dimension(:, :), intent(in) :: pos
    real(rk), dimension(world_dim), intent(in) :: cell
    real(rk), dimension(world_dim) :: ref_pos
    integer :: d

    ref_pos = pos(:, 1)
    do d = 1, world_dim
      pos_gathered(d, :) = pos(d, :) - ref_pos(d)
      pos_gathered(d, :) = ref_pos(d) + pos_gathered(d, :) - &
                           nint(pos_gathered(d, :) / cell(d)) * cell(d)
    end do
  end subroutine gatherMolPos

  subroutine sd_broadcastpos()
    if (r_group_idx == 0) then
      call mpi_scatterv(pos, scounts_c_world_dim, displs_c_world_dim, &
                        mpi_double_precision, pos_c, &
                        scounts_c_world_dim(c_group_idx + 1), &
                        mpi_double_precision, root, row_comm, ierr)
    end if
    call mpi_bcast(pos_c, scounts_c_world_dim(c_group_idx + 1), &
                   mpi_double_precision, root, col_comm, ierr)

    if (c_group_idx == 0) then
      call mpi_scatterv(pos, scounts_r_world_dim, displs_r_world_dim, &
                        mpi_double_precision, pos_r,&
                        scounts_r_world_dim(r_group_idx + 1), &
                        mpi_double_precision, root, col_comm, ierr)
    end if
    call mpi_bcast(pos_r, scounts_r_world_dim(r_group_idx + 1), &
                   mpi_double_precision, root, row_comm, ierr)
    deallocate(pos)
  end subroutine sd_broadcastpos

  subroutine sd_cal_num_rbin(cell)
    implicit none
    real(rk), intent(in) :: cell(world_dim)

    num_rbin = ceiling(cell(1) / 2d0 * sqrt(3d0) / rbinwidth)
    ! *sqrt(3) to accommodate the longest distance inside a cubic (diagonal)
    if (myrank == root) write(*,*) "num_rbin = ", num_rbin
  end subroutine sd_cal_num_rbin

  subroutine sd_prep_corrmemory(maxlag, nummoltype, numframe)
    implicit none
    integer, intent(in) :: maxlag, nummoltype, numframe
    integer :: num_moltypepair

    num_moltypepair = nummoltype * (nummoltype + 1) / 2
    allocate(sdcorr(maxlag+1, num_rbin, num_moltypepair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: sdcorr"
      call mpi_abend()
    end if
    sdcorr = 0d0

    allocate(sdpaircount(num_rbin, num_moltypepair), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: sdpaircount"
      call mpi_abend()
    end if
    sdpaircount = 0d0

    allocate(sd_binIndex(numframe), stat=stat)
    if (stat /= 0) then
      write(*,*) "Allocation failed: sd_binIndex"
      call mpi_abend()
    end if
  end subroutine sd_prep_corrmemory

  subroutine sd_getbinindex(r, c, cell, sd_binIndex)
    implicit none
    integer, intent(in) :: r, c
    real(rk), intent(in) :: cell(world_dim)
    integer, intent(out) :: sd_binIndex(:)
    real(rk) :: pp(world_dim, size(sd_binIndex))
    real(rk) :: ppd(size(sd_binIndex))
    integer :: d

    pp = pos_r(:,:,r) - pos_c(:,:,c)
    do d = 1, world_dim
      pp(d, :) = pp(d, :) - nint(pp(d, :) / cell(d)) * cell(d)
    end do
    ppd = sqrt(sum(pp*pp, 1))
    sd_binIndex = ceiling(ppd / rbinwidth)
    where (sd_binIndex == 0)
      sd_binIndex = 1
    end where
    if (do_diagonal) then
        where(not_diag(pp, ppd, od_tol)) sd_binIndex = -1
    else if (do_orthogonal) then
        where(not_orth(pp, ppd, od_tol)) sd_binIndex = -1
    end if
  end subroutine sd_getbinindex

  subroutine sd_collectcorr()
    implicit none
    if (myrank == root) then
      write(*,*) "collecting sdcorr"
      call mpi_reduce(MPI_IN_PLACE, sdcorr, size(sdcorr), mpi_double_precision, MPI_SUM, root, mpi_comm_world, ierr)
    else
      call mpi_reduce(sdcorr, dummy_null, size(sdcorr), mpi_double_precision, MPI_SUM, root, mpi_comm_world, ierr)
    end if
    call mpi_barrier(mpi_comm_world, ierr)

    if (myrank == root) then
      write(*,*) "collecting sdpaircount"
      call mpi_reduce(MPI_IN_PLACE, sdpaircount, size(sdpaircount), mpi_double_precision, MPI_SUM, root, mpi_comm_world, ierr)
    else
      call mpi_reduce(sdpaircount, dummy_null, size(sdpaircount), mpi_double_precision, MPI_SUM, root, mpi_comm_world, ierr)
    end if
    call mpi_barrier(mpi_comm_world, ierr)
  end subroutine sd_collectcorr

  subroutine sd_average(numframe, nummoltype, framecount)
    use utility, only: get_pairindex_upper_diag
    implicit none
    integer, intent(in) :: numframe, nummoltype, framecount(:)
    integer :: i, t1, t2, n, moltypepair_idx

    sdpaircount = sdpaircount / numframe
    do n = 1, nummoltype * (nummoltype + 1) / 2
      do i = 1, num_rbin
        sdcorr(:,i,n) = sdcorr(:,i,n) / framecount / sdpaircount(i, n)
      end do
    end do

    do t2 = 1, nummoltype
      do t1 = t2, nummoltype
        if (t1 /= t2) then
          moltypepair_idx = get_pairindex_upper_diag(t1, t2, nummoltype)
          sdpaircount(:, moltypepair_idx) = sdpaircount(:, moltypepair_idx) / 2d0
        end if
      end do
    end do
  end subroutine sd_average

  subroutine sd_make_rbins()
    implicit none
    integer :: i, stat
    allocate(rbins(num_rbin), stat=stat)
    if (stat /=0) then
      write(*,*) "Allocation failed: rbins"
      call mpi_abend()
    end if
    rbins = [ (i - 0.5d0, i = 1, num_rbin) ] * rbinwidth
  end subroutine sd_make_rbins

  subroutine sd_finish()
    implicit none
    !deallocate(pos_r, pos_c, sdcorr, sdpaircount)
    ! sd_binIndex has been deallocated in the main program
  end subroutine sd_finish

  function not_diag(r, rd, tol)
      implicit none
      real(rk), intent(in) :: r(:, :)    !(world_dim, num_frame)
      real(rk), intent(in) :: rd(:) !(num_frame)
      real(rk), intent(in) :: tol        !tolerance = sin^2(tol_angle)
      logical :: not_diag(size(r, 2))
      real(rk) :: c(world_dim, size(r, 2))
      real(rk) :: nr(world_dim, size(r, 2))

      nr = abs(r) / spread(rd, 1, world_dim)
      c(1, :) = nr(2, :) - nr(3, :)
      c(2, :) = nr(3, :) - nr(1, :)
      c(3, :) = nr(1, :) - nr(2, :)
      c = c / sqrt(3.0_rk)  !normalize to sin(theta)
      not_diag = sum(c * c, 1) > tol
  end function

  function not_orth(r, rd, tol)
      implicit none
      real(rk), intent(in) :: r(:, :)    !(world_dim, num_frame)
      real(rk), intent(in) :: rd(:) !(num_frame)
      real(rk), intent(in) :: tol        !tolerance = sin^2(tol_angle)
      logical :: not_orth(size(r, 2))
      real(rk) :: c(world_dim, size(r, 2))
      real(rk) :: nr(world_dim, size(r, 2))
      real(rk) :: orth(world_dim, world_dim)
      logical :: is_orth(size(r, 2))
      integer :: i

      nr = r / spread(rd, 1, world_dim)
      orth = reshape([1_rk, 0_rk, 0_rk, 0_rk, 1_rk, 0_rk, 0_rk, 0_rk, 1_rk], [3, 3])
      is_orth = .false.
      do i = 1, 3
          c = cross_prod(nr, orth(i, :))
          is_orth = is_orth .or. sum(c * c, 1) < tol
      end do
      not_orth = .not. is_orth
  end function

  pure function cross_prod(a, v)
    implicit none
    real(rk), intent(in) :: a(:, :), v(world_dim)
    real(rk) :: cross_prod(world_dim, size(a, 2))

    cross_prod(1, :) = a(2, :) * v(3) - a(3, :) * v(2)
    cross_prod(2, :) = a(3, :) * v(1) - a(1, :) * v(3)
    cross_prod(3, :) = a(1, :) * v(2) - a(2, :) * v(1)
  end function
end module spatial_dec
