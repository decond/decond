! variables and parameters
module varpars
  use, intrinsic :: iso_fortran_env, only: real64
  use top, only: system
  implicit none
  public
  integer, parameter:: rk = real64
  integer, parameter :: line_len = 1024
  character(len=*), parameter :: line_len_str = "1024"
  character(len=*), parameter :: decond_version = "0.5.2"
  character(len=*), parameter :: dec_mode_ec0 = "ec0"
  character(len=*), parameter :: dec_mode_ec1 = "ec1"
  character(len=*), parameter :: dec_mode_vsc = "vsc"
  character(len=*), parameter :: dec_mode_vel = "vel"
  integer, parameter :: world_dim = 3
  integer :: totnummol, sysnumatom
  character(len=line_len) :: corrfile, trjfile, dec_mode
  integer :: numframe, maxlag, nummoltype, skiptrj, num_moltypepair_all
  integer, allocatable :: charge(:), framecount(:)
  real(rk) :: cell(world_dim), timestep, temperature
  logical :: do_sd, do_ed, do_minus_avg
  type(system) :: sys
  integer :: qnt_dim
end module varpars
