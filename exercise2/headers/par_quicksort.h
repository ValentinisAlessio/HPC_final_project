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


// #if defined(_OPENMP)

// // measure the wall-clock time
// #define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
//                   (double)ts.tv_nsec * 1e-9)

// // measure the cpu thread time
// #define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec +     \
//                      (double)myts.tv_nsec * 1e-9)

// #else

// // measure ther cpu process time
// #define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
//                   (double)ts.tv_nsec * 1e-9)
// #endif


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
compare_t compare;      // compare function
compare_t compare_ge;   // compare function for "greater or equal"
verify_t verify_partitioning;
verify_t verify_sorting;
verify_t show_array;

// Declare partitioning and sorting functions
// static inline int partitioning(data_t * data, int start, int end, compare_t cmp_ge){
//     // Pick the median of the [0], [mid] and [end] element as pivot
//     int mid = (start + end-1) / 2;
//     if (cmp_ge((void*)&data[start], (void*)&data[mid]))
//         SWAP((void*)&data[start], (void*)&data[mid], sizeof(data_t));
//     if (cmp_ge((void*)&data[mid], (void*)&data[end-1]))
//         SWAP((void*)&data[mid], (void*)&data[end-1], sizeof(data_t));
//     if (cmp_ge((void*)&data[mid], (void*)&data[start]))
//         SWAP((void*)&data[start], (void*)&data[mid], sizeof(data_t));

//     // Pick the first element as pivot
//     void* pivot = (void*)&data[start];

//     // Partition around the pivot
//     int pointbreak = start + 1;

//     for (int i = start + 1; i < end; ++i){
//         if (!cmp_ge((void*)&data[i], pivot)){
            
//             // Move elements less than pivot to the left side
//             SWAP((void*)&data[i], (void*)&data[pointbreak], sizeof(data_t));

//             ++ pointbreak;
            
//         }
//     }

//     // Put the pivot in the right place
//     SWAP((void*)&data[start], (void*)&data[pointbreak - 1], sizeof(data_t));

//     // Return the pivot position
//     return pointbreak - 1;
// }

int partitioning(data_t *, int, int, compare_t);

// static inline int mpi_partitioning(data_t * data, int start, int end, compare_t cmp_ge, void* pivot){
//     // Function that partitions the array into two parts given a pivot
//     // and returns the position of the last element of the first part

//     // Partition around the pivot
//     int pointbreak = start;

//     // This can't be done in parallel because of possible data races in the exchanges and pointbreak increment
//     // Could be done with a parallel for loop, synchronized by an atomic increment of pointbreak, but it would be slower
//     for (int i = start; i < end; ++i){
//         if (!cmp_ge((void*)&data[i], pivot)){
            
//             // Move elements less than pivot to the left side
//             SWAP((void*)&data[i], (void*)&data[pointbreak], sizeof(data_t));

//             ++ pointbreak;
            
//         }
//     }

//     // We don't need to Put the pivot in the right place since the mpi pivot might not contain it!
//     // SWAP((void*)&data[start], (void*)&data[pointbreak - 1], sizeof(data_t));

//     // Return the pivot position
//     return pointbreak - 1;
// }

int mpi_partitioning(data_t *, int, int, compare_t, void*);

// Quicksort in distributed memory
void mpi_quicksort(data_t**, int*, MPI_Datatype, MPI_Comm, compare_t);

// Serial quicksort
void quicksort(data_t *, int, int, compare_t);

// Quicksort in shared memory
void par_quicksort(data_t *, int, int, compare_t);

// Global verification function
int verify_global_sorting(data_t*, int, int, MPI_Datatype, int, int, int);