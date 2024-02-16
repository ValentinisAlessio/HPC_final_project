#!/bin/bash
#SBATCH --job-name=omp_timings
#SBATCH --output=omp_timings.out
#SBATCH --error=omp_timings.err
#SBATCH --get-user-env
#SBATCH -p EPYC
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ALESSIO.VALENTINIS@studenti.units.it

date
pwd
hostname

module purge
module load architecture/AMD
module load openMPI/4.1.5/gnu/12.2.1

csv_file="data/timings$N.csv"


make


N = 64000000

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=4

./main $N

echo "Size,Threads,Time" > $csv_file

for i in {1..64}
do
    export OMP_NUM_THREADS=$i
    for j in {1..5}
    do 
        echo -n "$N,$i," >> $csv_file
        mpirun -np 1 --map-by node --bind-to socket ./main $N >> $csv_file
    done
done

make clean