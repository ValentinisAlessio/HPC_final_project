/*------------------------------------------------------------------
    par_quicksort.c - Parallel quicksort using OpenMP and MPI
    Author: Valentinis Alessio
    Last modified: February 2024

    Compile: mpicc -fopenmp -o par_quicksort par_quicksort.c
    Run: mpirun -np 4 ./par_quicksort
    Run with 4 threads: mpirun -np 4 ./par_quicksort 4
    Run with 8 threads: mpirun -np 4 ./par_quicksort 8

    The program will generate a random array of size 1000000 and sort it
    using quicksort. The array will be divided in 4 parts and each part
    will be sorted by a different process. The sorted array will be
    printed to the screen.

    The program will also measure the time it takes to sort the array
    and print it to the screen.

------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>

// ================================================================
//  MACROS and DATATYPES
// ================================================================


#if defined(_OPENMP)

// measure the wall-clock time
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)

// measure the cpu thread time
#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec +     \
                     (double)myts.tv_nsec * 1e-9)

#else

// measure ther cpu process time
#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)
#endif


#if defined(DEBUG)
#define VERBOSE
#endif

#if defined(VERBOSE)
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...)
#endif

#if !defined(DATA_SIZE)
#define DATA_SIZE 8
#endif
#define HOT       0

// let's define the default amount of data
//
#if (!defined(DEBUG) || defined(_OPENMP))
#define N_dflt    100000
#else
#define N_dflt    10000
#endif

// Define the data_t struct
typedef struct {
    double data[DATA_SIZE];
} data_t;

// Define macros for min and max between data_t objects
#define MIN(a,b) ( (a)->data[HOT] < (b)->data[HOT]? (a) : (b));
#define MAX(a,b) ( (a)->data[HOT] > (b)->data[HOT]? (a) : (b));

// ================================================================
//  FUNCTION PROTOTYPES
// ================================================================

// Define compare function that will be used by qsort
typedef int compare_t(const void *, const void *);

// Define verifying function type, used to test results
typedef int verify_t(data_t *, int, int, int);

// Declare the functions
compare_t compare;      // compare function
compare_t compare_ge;   // compare function for "greater or equal"
verify_t verify_partitioning;
verify_t verify_sorting;
verify_t show_array;

// Declare partitioning and sorting functions
int partition(data_t *, int, int, compare_t);
void quicksort(data_t *, int, int, compare_t);

void par_quicksort(data_t *, int, int, compare_t);

data_t* merge(data_t *, int, data_t *, int, compare_t);