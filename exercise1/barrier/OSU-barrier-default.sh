#!/bin/bash
# Script to automatize OSU-benchmark routines for boocking broadcast operation, varying the number of processes and the message size
#SBATCH --job-name=hpcex1
#SBATCH --output=barrier-default.out
#SBATCH --error=barrier-default.err
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH -p EPYC
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ALESSIO.VALENTINIS@studenti.units.it

# Load the MPI module
module load openMPI/4.1.5/gnu/12.2.1

# Specify the path to the result file
csv_file="../results_def/barrier_default_results_wu.csv"

# Go to the directory where the benchmark is located
src_path="../../../osu-micro-benchmarks-7.3/c/mpi/collective/blocking/"


# Define variables
#MESSAGE_SIZES=(1024 2048 4096 8192)  # Example message sizes
# np_values=(2 $(seq 8 8 256))  # Example number of processes

# Define different process map_values to evaluate
map_values="core socket node"

# Create CSV file and add headers
echo "Algorithm,Allocation,Processes,MessageSize,Avg Latency(us)" > $csv_file

# Loop through number of processes
for mapping in $map_values; do
    for np in {1..256}; do
        echo "Running MPI Bcast benchmark: map=$mapping, np=$np, broadcast_algo=$broadcast_algo ..."
        # Run MPI Bcast benchmark and capture output
        mpirun -np $np --map-by $mapping --mca coll_tuned_use_dynamic_rules true --mca coll_tuned_barrier_algorithm 0 $src_path/osu_barrier -x 100 -i 10000 |\
        # Append results to CSV file
        tail -n 1 | awk -v mapping="$mapping" -v np="$np" '{printf "default,%s,%s,%s,%s\n",mapping,np,$1,$2}' | sed 's/,$//' >> $csv_file
    done
done
