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
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>

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
#define N_dflt    10
#else
#define N_dflt    10
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

int main ( int argc, char **argv )
{

    
    // ---------------------------------------------
    //  get the arguments
    //


    int N          = N_dflt;
    
    /* check command-line arguments */
    {
        int a = 0;
        
        if ( argc > ++a ) N = atoi(*(argv+a));
    }
    
    // ---------------------------------------------
    //  generate the array
    //
    
    data_t *data = (data_t*)malloc(N*sizeof(data_t));
    long int seed;
    #if defined(_OPENMP)
    #pragma omp parallel
    {
        int me             = omp_get_thread_num();
        short int seed     = time(NULL) % ( (1 << sizeof(short int))-1 );
        short int seeds[3] = {seed-me, seed+me, seed+me*2};

    #pragma omp for
        for ( int i = 0; i < N; i++ )
        data[i].data[HOT] = erand48( seeds );
    }
    #else
    {
        seed = time(NULL);
        srand48(seed);
        
        PRINTF("ssed is % ld\n", seed);
        
        for ( int i = 0; i < N; i++ )
        data[i].data[HOT] = drand48();
    }    
    #endif

    #if defined(_OPENMP)
    // Try to sort the array
    // printf("Sorting array of size %d\n", N);
    // show_array(data, 0, N, 0);
    // printf("Size of data: %ld\n", sizeof(data));
    par_quicksort(data, 0, N, compare_ge);
    // printf("Sorted array:\n");
    // show_array(data, 0, N, 0); 
    #else
    // Try to sort the array
    quicksort(data, 0, N, compare_ge);
    #endif 
    
    int sorted = verify_sorting(data, 0, N, 0);
    //printf("Array is sorted: %s\n", sorted ? "true" : "false");


    printf("Array is sorted: %s\n", sorted ? "true" : "false");
    free( data );


    return 0;
}

#define SWAP(A,B,SIZE) do {int sz = (SIZE); char *a = (A); char *b = (B); \
do { char _temp = *a;*a++ = *b;*b++ = _temp;} while (--sz);} while (0)

int par_partition(data_t* data, int start, int end, compare_t cmp_ge){
    
    // // Pick the median of the [0], [mid] and [end] element as pivot
    // int mid = (start + end) / 2;
    // if(cmp_ge(data[start], data[mid]))
    //     SWAP(&data[start], &data[mid], sizeof(data_t));
    // if(cmp_ge(data[mid], data[end]))
    //     SWAP(&data[mid], &data[end], sizeof(data_t));
    // if(cmp_ge(data[start], data[end]))
    //     SWAP(&data[start], &data[end], sizeof(data_t));

    // Pick the first element as pivot
    void* pivot = (void*)&data[start];

    // Partition around the pivot
    int pointbreak = start + 1;

    // Parallelize this loop
    #pragma omp parallel for shared(data, pivot, pointbreak)
    for (int i = start + 1; i <= end; ++i){
        if (!cmp_ge((void*)&data[i], pivot)){
            
            // Move elements less than pivot to the left side
            SWAP((void*)&data[i], (void*)&data[pointbreak], sizeof(data_t));

            #pragma omp atomic
            ++ pointbreak;
            
        }
    }

    // Put the pivot in the right place
    SWAP((void*)&data[start], (void*)&data[pointbreak - 1], sizeof(data_t));

    // Return the pivot position
    return pointbreak - 1;
}

int partition(data_t* data, int start, int end, compare_t cmp_ge){
    
    // // Pick the median of the [0], [mid] and [end] element as pivot
    // int mid = (start + end) / 2;
    // if(cmp_ge(data[start], data[mid]))
    //     SWAP(&data[start], &data[mid], sizeof(data_t));
    // if(cmp_ge(data[mid], data[end]))
    //     SWAP(&data[mid], &data[end], sizeof(data_t));
    // if(cmp_ge(data[start], data[end]))
    //     SWAP(&data[start], &data[end], sizeof(data_t));

    // Pick the first element as pivot
    void* pivot = (void*)&data[start];

    // Partition around the pivot
    int pointbreak = start + 1;

    // Parallelize this loop
    for (int i = start + 1; i <= end; ++i){
        if (!cmp_ge((void*)&data[i], pivot)){
            
            // Move elements less than pivot to the left side
            SWAP((void*)&data[i], (void*)&data[pointbreak], sizeof(data_t));

            ++ pointbreak;
            
        }
    }

    // Put the pivot in the right place
    SWAP((void*)&data[start], (void*)&data[pointbreak - 1], sizeof(data_t));

    // Return the pivot position
    return pointbreak - 1;
}

// Quicksort algorithm
void quicksort(data_t* data, int start, int end, compare_t cmp_ge){

    #if defined(DEBUG)
    #define CHECK{\
        if (verify_partitioning(data, start, end, mid)){\
            printf("Partitioning error!\n");\
            printf("%4d, %4d (%4d, &g) -> %4d, %4d + %4d, %4d\n",\
            start, end, mid, data[mid].data[HOT], start, mid, mid+1, end);\
            show_array(data, start, end, 0); }}
    #else
    #define CHECK
    #endif


    int size = end - start;

    if(start < end){
        // Partition the array
        int pivot = partition(data, start, end, cmp_ge);

        CHECK;  // Verify partitioning

        // Sort the left and right side
        quicksort(data, start, pivot, cmp_ge);
        quicksort(data, pivot + 1, end, cmp_ge);
    }
    else{
        
        // Sort the array of size 2
        if ( (size == 2) && cmp_ge((void*)&data[start], (void*)&data[end-1]) )
            SWAP((void*)&data[start], (void*)&data[end-1], sizeof(data_t));
    }
}

// Parallel quicksort algorithm
void par_quicksort(data_t* data, int start, int end, compare_t cmp_ge){

    #if defined(DEBUG)
    #define CHECK{\
        if (verify_partitioning(data, start, end, mid)){\
            printf("Partitioning error!\n");\
            printf("%4d, %4d (%4d, &g) -> %4d, %4d + %4d, %4d\n",\
            start, end, mid, data[mid].data[HOT], start, mid, mid+1, end);\
            show_array(data, start, end, 0); }}
    #else
    #define CHECK
    #endif

    int size = end - start;
    if (size > 2){    
        // Partition the array
        int pivot = partition(data, start, end, cmp_ge);
        // int pivot = par_partition(data, start, end, cmp_ge);

        CHECK;  // Verify partitioning

        // Sort the left and right side
        #pragma omp task 
        {
        par_quicksort(data, start, pivot, cmp_ge);
        }
        #pragma omp task 
        {
        par_quicksort(data, pivot + 1, end, cmp_ge);
        }

    } else {
        // Sort the array of size 2
        if ( (size == 2) && cmp_ge((void*)&data[start], (void*)&data[end-1]) )
            SWAP((void*)&data[start], (void*)&data[end-1], sizeof(data_t));
    
    }
}

data_t* merge(data_t* data1, int n1, data_t* data2, int n2, compare_t cmp_ge){

    // Allocate memory for the merged array
    data_t* merged = (data_t*)malloc((n1 + n2)*sizeof(data_t));
    int i =0, j = 0, k;

    // Merge the two arrays
    for (k = 0; k < n1 + n2; k++){
        if (i>=n1){
            merged[k] = data2[j];
            j++;
        }
        else if (j>=n2){
            merged[k] = data1[i];
            i++;
        }

        // Indeces are in bounds
        // i < n1 && j < n2
        else if (cmp_ge((void*)&data1[i], (void*)&data2[j])){
            merged[k] = data1[i];
            i++;
        }
        else{
            merged[k] = data2[j];
            j++;
        }
    }

    free(data1);
    free(data2);

    return merged;
}

// =================================================================================================
// Verification functions (for debugging)
// =================================================================================================


int verify_sorting(data_t* data, int start, int end, int not_used){

    // Check if the array is sorted
    int i = start;
    while (++i < end && (data[i].data[HOT] >= data[i-1].data[HOT]));
    return (i == end);
}

int verify_partitioning(data_t* data, int start, int end, int mid){
    int failure = 0;
    int fail = 0;

    for (int i = start; i< mid; i++)
        if (compare((void*)&data[i], (void*)&data[mid]) >= 0 )
            fail ++;

    failure += fail;
    if (fail){
        printf("Left side error\n");
        fail = 0;
    }

    for (int i = mid + 1; i <= end; i++)
        if (compare((void*)&data[i], (void*)&data[mid]) < 0 )
            fail ++;
    
    failure += fail;
    if (fail){
        printf("Right side error\n");
    }

    return failure;
}

int show_array(data_t* data, int start, int end, int not_used){
    for (int i = start; i <= end; i++)
        printf("%f ", data[i].data[HOT]);
    printf("\n");
    return 0;
}

int compare(const void* a, const void* b){
    data_t* A = (data_t*)a;
    data_t* B = (data_t*)b;
    double diff = A->data[HOT] - B->data[HOT];

    // return 1 if A > B, 0 if A == B, -1 if A < B
    return ((diff > 0) - (diff < 0));
}

int compare_ge(const void* a, const void* b){
    data_t* A = (data_t*)a;
    data_t* B = (data_t*)b;

    // return 1 if A >= B, 0 if A < B
    return (A->data[HOT] >= B->data[HOT]);
}