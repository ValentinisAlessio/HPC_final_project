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

##exec="./main"
N=320000000
MPI_procs=32
OMP_threads=8

module load openMPI/4.1.5/gnu/12.2.1
module load architecture/AMD


make

export OMP_NUM_THREADS=1
echo "Serial run"
./main $N

export OMP_NUM_THREADS=$OMP_threads
echo "OMP run with $OMP_threads threads"
mpirun -np 1 --map-by node --bind-to socket ./main $N

echo "MPI run"
mpirun -np $MPI_procs --map-by socket ./main $N

echo "Processes,Threads,Time" > mpi_scaling.csv

export OMP_NUM_THREADS=1
for iter in {1..10..1}
do
	for p in 1 2 4 8 16 32 64
	do
		mpirun -np $p --map-by socket ./main $N >> mpi_scaling.csv
	done
done

make clean
