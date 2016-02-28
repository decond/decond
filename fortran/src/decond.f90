program decond
  use mpiproc
  !use top, only: system
  !use correlation, only: corr
  !use spatial_dec, only: rbinwidth, pos_r, pos_c, sdpaircount, &
                         !sdcorr, pos, num_rbin, com_pos, sd_binIndex, &
                         !sd_prep_corrmemory, sd_getbinindex, &
                         !sd_cal_num_rbin, sd_broadcastpos, sd_prepPosMemory, &
                         !sd_collectcorr, sd_average, sd_make_rbins, sd_finish
  !use energy_dec, only: engfiles, ed_readEng, ed_getbinindex, &
                        !num_ebin, ed_binIndex, ed_prep_corrmemory, &
                        !skipEng, ed_collectcorr, ed_average, &
                        !edpaircount, edcorr, ebinwidth, &
                        !ed_make_ebins, ed_finish, num_engfiles, &
                        !ed_prep_engfiles
  !use varpars, only: dp, decond_version, &
                     !totnummol, sysnumatom, &
                     !corrfile, datafile, &
                     !line_len_str, numframe, maxlag, nummoltype, &
                     !skiptrr, num_moltypepair_all, charge, framecount, &
                     !cell, timestep, temperature, &
                     !is_sd, is_ed, sys
  use manager, only: is_pa_mode, is_pm_mode, init_config, read_config, &
                     prepare, decompose, output, finish
  implicit none

  call mpi_setup('init')

  prog_starttime = MPI_Wtime()

  call init_config()

  call read_config()

  call prepare()

  call decompose()

  call output()

  if (myrank == root) write(*,*)
  if (myrank == root) write(*,*) "time for the whole program (sec):", MPI_Wtime() - prog_starttime

  call finish()
  call mpi_setup('stop')
  stop
end program decond
