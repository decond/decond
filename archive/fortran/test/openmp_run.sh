#PBS -l nodes=1:ppn=16
#PBS -N test_openmp

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8
./a.out
