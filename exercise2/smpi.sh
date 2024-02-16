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

csv_file="data/smpi_timings$N.csv"