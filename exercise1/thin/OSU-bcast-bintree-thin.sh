#!/bin/bash
# Script to automatize OSU-benchmark routines for boocking broadcast operation, varying the number of processes and the message size
#SBATCH --job-name=hpcex1
#SBATCH --output=OSU-bintree-bcast-thin.out
#SBATCH --error=OSU-bintree-bcast-thin.err
#SBATCH --nodes=2
#SBATCH --ntasks=48
#SBATCH -p THIN
#SBATCH --exclusive
#SBATCH --exclude=fat[001-002]
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ALESSIO.VALENTINIS@studenti.units.it

date
pwd
hostname

# Load the MPI module
module load openMPI/4.1.5/gnu/12.2.1

# Specify the path to the result file
csv_file="results/bcast_bintree_thin_results.csv"

# Go to the directory where the benchmark is located
src_path="../../../osu-micro-benchmarks-thin/c/mpi/collective/blocking/"


# Define variables
#MESSAGE_SIZES=(1024 2048 4096 8192)  # Example message sizes

# Define different process map_values to evaluate
map_values="core socket node"

# Create CSV file and add headers
echo "Algorithm,Allocation,Processes,MessageSize,Avg Latency(us)" > $csv_file

# Loop through process map_values
for mapping in $map_values; do
    # Loop through number of processes
    for np in {1..48}; do
        echo "Running MPI Bcast benchmark: map=$mapping, np=$np, broadcast_algo=$broadcast_algo ..."
        # Run MPI Bcast benchmark and capture output
        mpirun -np $np --map-by $mapping --mca coll_tuned_use_dynamic_rules true --mca coll_tuned_bcast_algorithm 5 $src_path/osu_bcast -x 100 -i 10000 |\
        # Append results to CSV file
        tail -n 21 | awk -v mapping="$mapping" -v np="$np" '{printf "bin_tree,%s,%s,%s,%s\n",mapping,np,$1,$2}'\
        | sed 's/,$//' >> $csv_file
    done
done
