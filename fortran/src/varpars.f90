! variables and parameters
module varpars
  use top, only: system
  implicit none
  public
  integer, parameter:: dp = kind(0.d0)  
  integer, parameter :: line_len = 1024
  character(len=*), parameter :: line_len_str = "1024"
  character(len=*), parameter :: decond_version = "0.4.6"
  integer :: totnummol, sysnumatom
  character(len=line_len) :: corrfile, datafile
  integer :: numframe, maxlag, nummoltype, skiptrr, num_moltypepair_all 
  integer, allocatable :: charge(:), framecount(:)
  real(dp) :: cell(3), timestep, temperature
  logical :: is_sd, is_ed
  type(system) :: sys
end module varpars
