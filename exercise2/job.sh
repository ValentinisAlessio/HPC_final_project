#!/bin/bash
#SBATCH --job-name=job
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --get-user-env
#SBATCH -p EPYC
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ALESSIO.VALENTINIS@studenti.units.it

date
pwd
hostname

module load openMPI/4.1.5/gnu/12.2.1
module load architecture/AMD


make

OMP_THREADS=8

export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "Serial run"
export OMP_NUM_THREADS=1
./main 10000000

echo "OMP run"
export OMP_NUM_THREADS=$(OMP_THREADS)
mpirun -np 1 ./main 10000000

echo "MPI run"
mpirun -np 16 ./main 10000000


make clean
