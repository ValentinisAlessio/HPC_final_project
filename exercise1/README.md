# ASSIGNMENT 1

Download library

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