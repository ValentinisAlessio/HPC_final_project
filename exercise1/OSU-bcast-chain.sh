#!/bin/bash
# Script to automatize OSU-benchmark routines for boocking broadcast operation, varying the number of processes and the message size
#SBATCH --job-name=hpcex1
#SBATCH --output=OSU-chain-bcast.out
#SBATCH --error=OSU-chin-bcast.err
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
csv_file="results/bcast_chain_results_wu_1.csv"

# Go to the directory where the benchmark is located
src_path="../../osu-micro-benchmarks-7.3/c/mpi/collective/blocking/"


# Define variables
#MESSAGE_SIZES=(1024 2048 4096 8192)  # Example message sizes
np_values=(2) 
np_values+=$(seq 8 8 256)  # Example number of processes

# Define different process map_values to evaluate
map_values="core socket node"

# Create CSV file and add headers
echo "Algorithm,Allocation,Processes,MessageSize,Avg Latency(us)" > $csv_file

# Loop through process map_values
for mapping in $map_values; do
    # Loop through number of processes
    for np in $np_values; do
        echo "Running MPI Bcast benchmark: map=$mapping, np=$np, broadcast_algo=$broadcast_algo ..."
        # Run MPI Bcast benchmark and capture output
        mpirun -np $np --map-by $mapping --mca coll_tuned_use_dynamic_rules true --mca coll_tuned_bcast_algorithm 2 $src_path/osu_bcast -x 100 -i 10000 |\
        # Append results to CSV file
        tail -n 21 | awk -v mapping="$mapping" -v np="$np" '{printf "chain,%s,%s,%s,%s\n",mapping,np,$1,$2}' | sed 's/,$//' >> $csv_file
    done
done
