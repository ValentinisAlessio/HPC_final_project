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
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <mpi.h>

// ================================================================
//  MACROS and DATATYPES
// ================================================================

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

#define L1_CACHE  65500 // Number of array elements that fit in a L1 cache of an Epyc Core

// Define the data_t struct
typedef struct {
    double data[DATA_SIZE];
} data_t;

// Define macros for min and max between data_t objects
#define MIN(a,b) ( (a)->data[HOT] < (b)->data[HOT]? (a) : (b));
#define MAX(a,b) ( (a)->data[HOT] > (b)->data[HOT]? (a) : (b));
#define SWAP(A,B,SIZE) do {int sz = (SIZE); char *a = (A); char *b = (B); \
do { char _temp = *a;*a++ = *b;*b++ = _temp;} while (--sz);} while (0)

// ================================================================
//  FUNCTION PROTOTYPES
// ================================================================

// Define compare function that will be used by qsort
typedef int compare_t(const void *, const void *);

// Define verifying function type, used to test results
typedef int verify_t(data_t *, int, int, int);

// Declare the functions
compare_t compare_ge;   // compare function for "greater or equal"
verify_t verify_sorting;
verify_t show_array;


int partitioning(data_t *, int, int, compare_t);


int mpi_partitioning(data_t *, int, int, compare_t, void*);

// Quicksort in distributed memory
void mpi_quicksort(data_t**, int*, MPI_Datatype, MPI_Comm, compare_t);

// Serial quicksort
void quicksort(data_t *, int, int, compare_t);

// Quicksort in shared memory
void par_quicksort(data_t *, int, int, compare_t);

// Global verification function
int verify_global_sorting(data_t*, int, int, MPI_Datatype, int, int, int);