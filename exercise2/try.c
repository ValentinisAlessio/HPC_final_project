#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <mpi.h>

#if !defined(DATA_SIZE)
#define DATA_SIZE 8
#endif
#define HOT       0

#define SIZE 100

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
data_t* merge(data_t*, int, data_t*, int, compare_t);
data_t* global_merge(data_t*, int, int, MPI_Datatype, data_t*, int);

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
        int nthreads = atoi(env_var);
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
    
    if (rank == 0){
        printf("Generating array of size %d\n", N);
        printf("Array before sorting:\n");
        show_array(data, 0, N, 0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // ---------------------------------------------
    // Create custom MPI data type for data_t
    MPI_Datatype MPI_DATA_T;
    MPI_Type_contiguous(sizeof(data_t), MPI_BYTE, &MPI_DATA_T);
    MPI_Type_commit(&MPI_DATA_T);

    // ---------------------------------------------
    // Divide the array into chunks and scatter them to the processes
    // Generate local array to store the scattered data of the right size
    int chunk_size = (rank < N % num_processes) ? N / num_processes + 1 : N / num_processes;
    printf("Process %d has chunk size %d\n", rank, chunk_size);
    data_t* loc_data = (data_t*)malloc(chunk_size*sizeof(data_t));
    printf("Process %d has allocated %d bytes\n", rank, chunk_size*sizeof(data_t));
    divide(data, 0, N, MPI_DATA_T, loc_data, num_processes);

    // Print scattered arrays
    for (int i = 0; i < num_processes; i++){
        if (rank == i){
            printf("Process %d received:\n", rank);
            show_array(loc_data, 0, chunk_size, 0);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // ---------------------------------------------
    // Sort the local array
    par_quicksort(loc_data, 0, chunk_size, compare_ge);

    // Print sorted arrays
    for (int i = 0; i < num_processes; i++){
        if (rank == i){
            printf("Process %d sorted:\n", rank);
            show_array(loc_data, 0, chunk_size, 0);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // ---------------------------------------------
    // Try merge on two arrays
    if (rank == 0){
        data_t* left = (data_t*)malloc(4*sizeof(data_t));
        data_t* right = (data_t*)malloc(4*sizeof(data_t));
        for (int i = 0; i < 4; i++){
            left[i].data[HOT] = 2*i +1;
            right[i].data[HOT] = 2*i +2;
        }
        printf("Left array:\n");
        show_array(left, 0, 4, 0);
        printf("Right array:\n");
        show_array(right, 0, 4, 0);
        data_t* merged = merge(left, 4, right, 4, compare_ge);
        printf("Merged array:\n");
        show_array(merged, 0, 8, 0);
    }

    // data_t* merged = (data_t*)malloc((N)*sizeof(data_t));
    // for (int i = 0; i < N; i++)
    //     merged[i].data[HOT] = 1.;
    // for (int i = 0; i < num_processes; i++){
    //     if (rank == i)
    //     merged = merge(merged, N, loc_data, chunk_size, compare_ge);
    // }

    // Merge all the local arrays into a single array
    data_t* merged = (data_t*)malloc((N)*sizeof(data_t));
    merged = global_merge(merged, 0, N, MPI_DATA_T, loc_data, num_processes);

    // ---------------------------------------------
    if (rank == 0){
        printf("Merged array:\n");
        show_array(merged, 0, N, 0);
    }

    MPI_Type_free(&MPI_DATA_T);
    MPI_Finalize();
    free(data);
    free(loc_data);
    free(merged);

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

data_t* global_merge(data_t* data, int start, int end, MPI_Datatype MPI_DATA_T, data_t* loc_data, int num_procs){
    // Function that merges the sorted arrays of all the processes
    // into a single sorted array
    // Returns a pointer to the merged array

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

    MPI_Gatherv(loc_data, sendcounts[0], MPI_DATA_T, data, sendcounts, displs, MPI_DATA_T, 0, MPI_COMM_WORLD);

    free(sendcounts);
    free(displs);
    
    return data;
}

int verify_sorting( data_t *data, int start, int end, int not_used )
{
    int i = start;
    while( (i <= end) && (data[i].data[HOT] >= data[i-1].data[HOT]) )
        i++;
    return ( i == end );
}