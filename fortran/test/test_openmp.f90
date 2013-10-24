program test_openmp
  use omp_lib
  implicit none
  integer :: id, nprocs, nthrds

!  call OMP_set_num_threads(4)
  !$OMP PARALLEL PRIVATE(id)
    nprocs = OMP_get_num_procs()
    nthrds = OMP_get_num_threads()
    id = OMP_get_thread_num()
    if (id == 0) then
      write(*,*) "nprocs=", nprocs
      write(*,*) "nthrds=", nthrds
    end if
    write(*,*) "Hello from thread id ", id
  !$OMP END PARALLEL
end program test_openmp
