#!/bin/bash
#SBATCH --job-name=job
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --get-user-env
#SBATCH -p EPYC
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=00:30:00
#SBATCH --mem=10000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ALESSIO.VALENTINIS@studenti.units.it

date
pwd
hostname

module load openMPI/4.1.5/gnu/12.2.1
module load architecture/AMD


make

export OMP_NUM_THREADS=16
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mpirun -np 8 ./main 100000

make clean
