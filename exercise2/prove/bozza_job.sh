#!/bin/bash
#SBATCH --job-name=parallel_sorting
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00
#SBATCH --output=parallel_sorting_output.%J.out
#SBATCH --error=parallel_sorting_output.%J.err

# Load necessary modules
module load mpi

# Set the maximum number of processes
MAX_PROCESSES=8

# Define the range of problem sizes for weak scalability
PROBLEM_SIZES=(100 200 400 800 1600)

# Define the problem size for strong scalability
PROBLEM_SIZE=1000

# Define the output file
OUTPUT_FILE="sorting_results.csv"

# Write the header to the CSV file
echo "Size,Processes,Threads,Time" > $OUTPUT_FILE

# Compile your program
mpicc -o parallel_sorting your_program.c

# Run the weak scalability test
for SIZE in "${PROBLEM_SIZES[@]}"
do
    echo "Weak Scalability Test - Problem Size: $SIZE"
    mpirun -np $MAX_PROCESSES ./parallel_sorting $SIZE >> $OUTPUT_FILE
done

# Run the strong scalability test
echo "Strong Scalability Test - Problem Size: $PROBLEM_SIZE"
echo "1,$PROBLEM_SIZE,1" >> $OUTPUT_FILE
mpirun -np 1 ./parallel_sorting $PROBLEM_SIZE >> $OUTPUT_FILE
for ((i=2; i<=$MAX_PROCESSES; i*=2))
do
    echo "$i,$PROBLEM_SIZE,$i" >> $OUTPUT_FILE
    mpirun -np $i ./parallel_sorting $PROBLEM_SIZE >> $OUTPUT_FILE
done
