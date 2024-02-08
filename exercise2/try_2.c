#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <mpi.h>

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

#if !defined(DATA_SIZE)
#define DATA_SIZE 8
#endif
#define HOT       0

#define SIZE 100

#if defined(DEBUG)
#define VERBOSE
#endif

#if defined(VERBOSE)
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...)
#endif

typedef struct {
    double data[SIZE];
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

int partition(data_t*, int, int, compare_t);
void par_quicksort(data_t*, int, int, compare_t);
void divide(data_t*, int, int, MPI_Datatype, data_t*, int);
void exchange();
data_t* merge(data_t*, int, data_t*, int, compare_t);
// data_t* global_merge(data_t*, int, int, MPI_Datatype, data_t*, int);

int main(int argc, char** argv){

    // Default values
    int N = SIZE;
    int nthreads=1;
    
    /* check command-line arguments */
    {
        int a = 0;
        
        if ( argc > ++a ) N = atoi(*(argv+a));
    }

    char* env_var = getenv("OMP_NUM_THREADS");
    if (env_var != NULL) {
        nthreads = atoi(env_var);
    } else {
        printf("OMP_NUM_THREADS environment variable not set.\n");
    }

    int num_processes, rank;
    int mpi_err = MPI_Init(&argc, &argv);

    if (mpi_err != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, mpi_err);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ---------------------------------------------
    //  generate the array
    //
    if (rank == 0){
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
    
        printf("Generating array of size %d\n", N);
        printf("Array before sorting:\n");
        show_array(data, 0, N, 0);

        int pivot = partition(data, 0, N, compare_ge);
        printf("Pivot is %d\n", pivot);
    }
    struct timespec ts;
    double tstart = CPU_TIME;
    double tend;

if (num_processes > 1){ // Run this section only if there are more than 1 processes
    printf("Running on %d processes\n", num_processes);
    // ---------------------------------------------
    // Create custom MPI data type for data_t
    MPI_Datatype MPI_DATA_T;
    MPI_Type_contiguous(sizeof(data_t), MPI_BYTE, &MPI_DATA_T);
    MPI_Type_commit(&MPI_DATA_T);


    // ---------------------------------------------
    // Broadcast the Size to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Compute chunk size
    int chunk_size = (N % num_processes == 0) ? (N / num_processes) : N / (num_processes - 1);

    // Calculate total size of chunk according to bits
    data_t* chunk = (data_t*)malloc(chunk_size*sizeof(data_t));

    // Scatter the chuck size data to all process
    MPI_Scatter(data, chunk_size, MPI_DATA_T, chunk, chunk_size, MPI_DATA_T, 0, MPI_COMM_WORLD);
    free(data);

    // Compute size of own chunk and then sort them using quick sort

    int own_chunk_size = (N >= chunk_size * (rank + 1)) ? chunk_size : (N - chunk_size * rank);

    // Sorting array with quick sort for every chunk as called by process
    #pragma omp parallel
    {
        #pragma omp master
            par_quicksort(chunk, 0, own_chunk_size, compare_ge);
    }

    // ---------------------------------------------
    // Return the sorted data to the master process by recursively merging two sorted arrays
    for (int step = 1; step < num_processes; step = 2 * step) {

        if (rank % (2 * step) != 0) {
            MPI_Send(chunk, own_chunk_size, MPI_DATA_T, rank - step, 0, MPI_COMM_WORLD);
            break;
        }

        if (rank + step < num_processes) {
            int received_chunk_size = (N >= chunk_size * (rank + 2 * step)) ? chunk_size * step : (N - chunk_size * (rank + step));
            data_t* received_chunk = (data_t*)malloc(received_chunk_size*sizeof(data_t));

            MPI_Recv(received_chunk, received_chunk_size, MPI_DATA_T, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Merge the received array with own array
            data = merge(chunk, own_chunk_size, received_chunk, received_chunk_size, compare_ge);

            // Free the memory
            free(chunk);
            free(received_chunk);

            // Update the chunk pointer
            chunk = data;
            own_chunk_size += received_chunk_size;
        }
    }
    // Return the sorted array to the data pointer
    data = chunk;

    // Release chunk
    free(chunk);

    MPI_Type_free(&MPI_DATA_T);
    tend = CPU_TIME;
    // Print the sorted array
    if (rank == 0){
        printf("Array after sorting:\n");
        show_array(data, 0, N, 0);
    }

}else{  // Run this section if there is only 1 process
    printf("Running on 1 process\n");
    #pragma omp parallel
    {
        #pragma omp master
        par_quicksort(data, 0, N, compare_ge);
    }
    tend = CPU_TIME;
    printf("Array after sorting:\n");
    show_array(data, 0, N, 0);
}

    // ---------------------------------------------

    // ---------------------------------------------
    // Finalize MPI
    MPI_Finalize();
    // free(chunk);
    // free(data);


    if (rank == 0)
        if ( verify_sorting( data, 0, N, 0) )
            printf("%d\t%g sec\n", nthreads, tend-tstart);
        else
            printf("the array is not sorted correctly\n");


    return 0;
}

#define SWAP(A,B,SIZE) do {int sz = (SIZE); char *a = (A); char *b = (B); \
do { char _temp = *a;*a++ = *b;*b++ = _temp;} while (--sz);} while (0)

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
    for (int i = start + 1; i < end; ++i){
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

    if (start < end){    
        // Partition the array
        int pivot = partition(data, start, end, cmp_ge);
        // int pivot = par_partition(data, start, end, cmp_ge);

        CHECK;  // Verify partitioning

        // Sort the left and right side
        #pragma omp task shared(data)
        par_quicksort(data, start, pivot, cmp_ge);
        #pragma omp task shared(data)
        par_quicksort(data, pivot + 1, end, cmp_ge);

    }
}

int show_array(data_t* data, int start, int end, int not_used){
    for (int i = start; i < end; i++)
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

void divide(data_t* data, int start, int end, MPI_Datatype MPI_DATA_T, data_t* loc_data, int num_procs){
    // Function that, given an array, divides it into num_procs parts
    // and scatters them to the processes (including the master process)
    // taking care of the case where the number of processes is not a
    // multiple of the array size

    int size = end - start;
    int chunk_size = size / num_procs;
    int remainder = size % num_procs;

    int* sendcounts = (int*)malloc(num_procs*sizeof(int));
    int* displs = (int*)malloc(num_procs*sizeof(int));

    for (int i = 0; i < num_procs; i++){
        sendcounts[i] = chunk_size;
        if (remainder > 0){
            sendcounts[i] += 1;
            remainder -= 1;
        }
        displs[i] = start;
        start += sendcounts[i];
    }

    MPI_Scatterv(data, sendcounts, displs, MPI_DATA_T, loc_data, sendcounts[0], MPI_DATA_T, 0, MPI_COMM_WORLD);

    free(sendcounts);
    free(displs);
}

data_t* merge(data_t* left, int left_size, data_t* right, int right_size, compare_t cmp_ge){
    // Function that merges two sorted arrays into a single sorted array
    // Returns a pointer to the merged array

    data_t* merged = (data_t*)malloc((left_size + right_size)*sizeof(data_t));
    int i = 0, j = 0, k = 0;

    while (i < left_size && j < right_size) {
        // Compare the elements from the left and right arrays
        if (cmp_ge((void*)&left[i], (void*)&right[j])) {
            merged[k] = right[j];
            j++;
        } else {
            merged[k] = left[i];
            i++;
        }
        k++;
    }

    // Copy the remaining elements from the left array
    while (i < left_size) {
        merged[k] = left[i];
        i++;
        k++;
    }

    // Copy the remaining elements from the right array
    while (j < right_size) {
        merged[k] = right[j];
        j++;
        k++;
    }

    return merged;
}

int verify_sorting( data_t *data, int start, int end, int not_used )
{
    int i = start;
    while( (i <= end) && (data[i].data[HOT] >= data[i-1].data[HOT]) )
        i++;
    return ( i == end );
}