! variables and parameters
module varpars
  implicit none
  public
  integer, parameter:: dp = kind(0.d0)  
  integer, parameter :: line_len = 128
  character(len=*), parameter :: line_len_str = "128"
  character(len=*), parameter :: decond_version = "0.4.6"
  integer, parameter :: num_positional_arg = 3, least_required_num_arg = 7

  integer :: num_arg, num_subarg, num_arg_per_moltype
  integer :: totnummol, sysnumatom
  character(len=line_len) :: corrfile, datafile, logfile, topfile, arg
  integer :: numframe, maxlag, nummoltype, skiptrr, num_moltypepair_all 
  integer, allocatable :: charge(:), framecount(:), start_index(:)
  real(dp) :: cell(3), timestep, temperature
  logical :: is_periodic, is_pa_mode, is_pm_mode, is_sd, is_ed

end module varpars
