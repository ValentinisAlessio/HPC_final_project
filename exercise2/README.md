# Exercise 2b

This folder contains the necessary scripts and source files useful to run tests on my parallel implementation of the quicksort algorithm.

## Folder organization
```
📂 exercise2/
│ 
├── 📂 data/
│
├── 📂 objects/   
│ 
├── 📂 headers/
│   └── 📄 par_quicksort.h
│
├── 📂 source/
│   ├── 📄 main.c
│   └── 📄 par_quicksort.c
│
├── 📂 plots/
│   ├── 📄 omp.ipynb
│   ├── 📄 wmpi.ipynb
│   └── 📄 smpi.ipynb
│
├── 📄 Makefile
├── 📄 omp_timings.sh
├── 📄 wmpi.sh
├── 📄 smpi.sh
└── 📰 README.md
```

Each folder name has a significant meaning:
- `data`: contains `.csv` files obtained by testing
- `objects`: contains objects files generated during compilation
- `plots`: contains jupyter notebook used to obtain the plots
- `headers` and `source`: contain headers and source files used to develop my program.

## How to compile

Authomatic compilation is provided through the use of a `Makefile`. Tests and data collection are authomatized through the various `.sh` SBATCH scripts.

If you want to run a test on the ORFEO cluster, just perform `sbatch <job>.sh`, while if you want just to perform compilation of the script, just do `make`.
