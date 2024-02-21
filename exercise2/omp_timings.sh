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

N=240000000

csv_file="data/timings$N.csv"


make



export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=4


echo "Size,Threads,Time" > $csv_file

for i in {1..64}
do
    export OMP_NUM_THREADS=$i
    for j in {1..5}
    do 
	mpirun -np 1 --map-by node --bind-to socket ./main $N | tail -n 1 | awk -v N="$N" -v nthr="$i" '{printf "%s,%s,%s\n",N,nthr,$1}' >> $csv_file
        ##echo -n "$N,$i,$out" >> $csv_file
    done
done

make clean
