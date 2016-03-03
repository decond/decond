! variables and parameters
module varpars
  use, intrinsic :: iso_fortran_env, only: real64
  use top, only: system
  implicit none
  public
  integer, parameter:: rk = real64
  integer, parameter :: line_len = 1024
  character(len=*), parameter :: line_len_str = "1024"
  character(len=*), parameter :: decond_version = "0.4.6"
  character(len=*), parameter :: dec_mode_ec0 = "ec0"
  character(len=*), parameter :: dec_mode_ec1 = "ec1"
  character(len=*), parameter :: dec_mode_vsc = "vsc"
  integer, parameter :: world_dim = 3
  integer :: totnummol, sysnumatom
  character(len=line_len) :: corrfile, trjfile, dec_mode
  integer :: numframe, maxlag, nummoltype, skiptrj, num_moltypepair_all
  integer, allocatable :: charge(:), framecount(:)
  real(rk) :: cell(world_dim), timestep, temperature
  logical :: do_sd, do_ed
  type(system) :: sys
  integer :: qnt_dim
end module varpars
