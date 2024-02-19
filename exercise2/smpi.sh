#!/bin/bash
#SBATCH --job-name=smpi2_timings
#SBATCH --output=smpi2_timings.out
#SBATCH --error=smpi2_timings.err
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
csv_file="data/smpi2_timings$N.csv"

make

echo "Size,Cores,Time" > $csv_file

export OMP_NUM_THREADS=1
for i in {1..5}
do
	./main $N | tail -n 1 | awk -v N="$N" -v nproc="1" '{printf "%s,%s,%s\n",N,nproc,$1}' >> $csv_file
done

export OMP_NUM_THREADS=4


#for i in {2..64}
for i in {65..128}
do
    for j in {1..5}
    do 
        mpirun -np $i --map-by socket ./main $N | tail -n 1 | awk -v N="$N" -v nproc="$i" '{printf "%s,%s,%s\n",N,nproc,$1}' >> $csv_file
    done
done

make clean
