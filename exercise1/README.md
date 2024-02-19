# ASSIGNMENT 1

This folder contains the scripts and data collected for the study of exercise 1 assignment.
The aim of the assignmento was to study the different types of implementations of the `MPI_Bcast` and one other collective algorithm of the `MPI` library. The report is available in its own folder.

## Folder organization
```
📂 exercise1/
│ 
├── 📂 barrier/
│   └── 📄 barrier_scripts.sh
│ 
├── 📂 results_def/
│   └── 📊 collected_data.csv
│
├── 📂 bcast/	
│   └── 📄 bcast_scripts.sh
│
├── 📂 source/
│   ├── 📝 plots_barrier.ipynb
│   └── 📄 plots_bcast.ipynb
│
├── 📂 models/
│   ├── 📄 barrier_model.R
│   └── 📄 bcast_model.R
│
└── 📰 README.md
```

## How to run the tests

First, download the OSU-microbenchmark library

```bash
wget https://mvapich.cse.ohio-state.edu/download/mvapich/osu-micro-benchmarks-7.3.tar.gz
tar -xzvf osu-micro-benchmarks-7.3.tar.gz
```

Then go into the directory and
```bash
./configure CC=/path/to/mpicc
make
make install
```

To make this process onto the target machine (mandatory in order to achieve targeted compilation), just put `srun` before each command after having allocated the appropriate resources.

The data collection process is fully authomatize via SBATCH scripts available into the `bcast` and `barrier` folders.