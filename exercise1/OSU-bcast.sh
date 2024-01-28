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

# Define variables
MESSAGE_SIZES=(1024 2048 4096 8192)  # Example message sizes
NUM_PROCESSES=(4 8)  # Example number of processes

# Define the MPI broadcast algorithms to test
BROADCAST_ALGORITHMS=(0 1 2)  # Example broadcast algorithms

# Define different process mappings to evaluate
MAPPINGS=("core" "socket" "numa" "board")

# Create CSV file and add headers
echo "Message Size,Processes,Mapping,Broadcast Algorithm,Avg Latency(us)" > bcast_results.csv

# Loop through message sizes
for size in "${MESSAGE_SIZES[@]}"; do
    # Loop through number of processes
    for np in "${NUM_PROCESSES[@]}"; do
        # Loop through process mappings
        for mapping in "${MAPPINGS[@]}"; do
            # Loop through broadcast algorithms
            for broadcast_algo in "${BROADCAST_ALGORITHMS[@]}"; do
                echo "Running MPI Bcast benchmark: Message size $size, Processes $np, Mapping $mapping, Broadcast Algorithm $broadcast_algo"
                # Run MPI Bcast benchmark and capture output
                output=$(mpirun -np $np --map-by $mapping --mca coll_tuned_use_dynamic_rules true --mca coll_tuned_bcast_algorithm $broadcast_algo osu_bcast)
                # Extract average latency from output
                avg_latency=$(echo "$output" | awk '{print $3}')
                # Append results to CSV file
                echo "$size,$np,$mapping,$broadcast_algo,$avg_latency" >> bcast_results.csv
            done
        done
    done
done