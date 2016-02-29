! variables and parameters
module varpars
  use top, only: system
  implicit none
  public
  integer, parameter:: dp = kind(0.d0)  
  integer, parameter :: line_len = 1024
  character(len=*), parameter :: line_len_str = "1024"
  character(len=*), parameter :: decond_version = "0.4.6"
  character(len=*), parameter :: dec_mode_ec0 = "ec0"
  character(len=*), parameter :: dec_mode_ec1 = "ec1"
  character(len=*), parameter :: dec_mode_vsc = "vsc"
  integer :: totnummol, sysnumatom
  character(len=line_len) :: corrfile, trjfile, dec_mode
  integer :: numframe, maxlag, nummoltype, skiptrj, num_moltypepair_all
  integer, allocatable :: charge(:), framecount(:)
  real(dp) :: cell(3), timestep, temperature
  logical :: do_sd, do_ed
  type(system) :: sys
end module varpars
