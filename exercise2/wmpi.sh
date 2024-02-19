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

N=2400000

csv_file="data/wmpi2_timings$N.csv"

make

echo "Size,Cores,Time" > $csv_file


export OMP_NUM_THREADS=1

for i in {1..5}
do
	./main $N | tail -n 1 | awk -v N="$N" -v nproc="1" '{printf "%s,%s,%s\n",N,nproc,$1}' >> $csv_file
done

export OMP_NUM_THREADS=4

#for i in {1..64}
for i in {65..128}
do
    for j in {1..5}
    do 
	let size=$((N*i))
        mpirun -np $i --map-by socket ./main $size | tail -n 1 | awk -v N="$size" -v nproc="$i" '{printf "%s,%s,%s\n",N,nproc,$1}' >> $csv_file
    done
done

make clean
