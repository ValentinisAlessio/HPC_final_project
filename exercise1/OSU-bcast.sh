#!/bin/bash
# Script to automatize OSU-benchmark routines for boocking broadcast operation, varying the number of processes and the message size
#SBATCH --job-name=OSU-bcast
#SBATCH --output=OSU-bcast.out
#SBATCH --error=OSU-bcast.err
#SBATCH --nodes=2
#SBATCH --ntasks=256
#SBATCH -p EPYC
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL

# Load the MPI module
module load openMPI/4.1.5/gnu/12.2.1

# Specify the path to the result file
csv_file="results/bcast_results.csv"

# Go to the directory where the benchmark is located
cd ../../osu-micro-benchmarks-7.3/c/mpi/collective/blocking/


# Define variables
#MESSAGE_SIZES=(1024 2048 4096 8192)  # Example message sizes
np_values="2 4 8 16 32 64 128 256"  # Example number of processes

# Define the MPI broadcast algorithms to test
bcast_algorithm="0 1 2 3"  # Example broadcast algorithms

# Define different process map_values to evaluate
map_values="core socket node"

# Create CSV file and add headers
echo "Algorithm,Allocation,Processes,MessageSize,Avg Latency(us)" > csv_file

# Loop through process map_values
for mapping in $map_values; do
    # Loop through number of processes
    for np in $np_values; do
        # Loop through broadcast algorithms
        for broadcast_algo in $bcast_algorithm; do
            echo "Running MPI Bcast benchmark: map=$mapping, np=$np, broadcast_algo=$broadcast_algo ..."
            # Run MPI Bcast benchmark and capture output
            output=$(mpirun -np $np --map-by $mapping --mca coll_tuned_use_dynamic_rules true --mca coll_tuned_bcast_algorithm $broadcast_algo osu_bcast)
            # Append results to CSV file
            echo "$output" | awk -v mapping="$mapping" -v np="$np" -v algo="$broadcast_algo" 'NR>2 {print algo "," mapping "," np "," $1 "," $2}' >> $csv_file
        done
    done
done
