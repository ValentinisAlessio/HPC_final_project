#!/bin/bash
#SBATCH --job-name=wmpi_timings
#SBATCH --output=wmpi_timings.out
#SBATCH --error=wmpi_timings.err
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

csv_file="data/wmpi_timings$N.csv"

make


N = 10000000

export OMP_PLACES=threads
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=2

./main $N

echo "Size,Cores,Time" > $csv_file

for i in {1..64}
do
    for j in {1..5}
    do 
        let $size = $(($N * $i))
        echo -n "$size,$i," >> $csv_file
        mpirun -np $i --map-by socket ./main $size >> $csv_file
    done
done

make clean