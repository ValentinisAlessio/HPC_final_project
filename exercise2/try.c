#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <mpi.h>
#include <string.h>

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
void gather(data_t*, int, int, MPI_Datatype, data_t*, int, int);
int mpi_partition(data_t*, int, int, compare_t, void*);
void mpi_quicksort(data_t**, int*, int, int, int, MPI_Datatype, compare_t);
void mpi_quicksort1(data_t**, int*, MPI_Datatype, MPI_Comm);
// void mpi_quicksort(data_t*, int*, int, int, int, MPI_Datatype, compare_t);
int verify_global_sorting(data_t*, int, int, MPI_Datatype, int, int, int);
int run(data_t*, int, int, MPI_Datatype, int, int);

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
    MPI_Datatype MPI_DATA_T;
    MPI_Type_contiguous(sizeof(data_t), MPI_BYTE, &MPI_DATA_T);
    MPI_Type_commit(&MPI_DATA_T);

    // ---------------------------------------------
    //  generate the array
    //

    // start the clock
    int chunk_size = (rank < N % num_processes) ? N / num_processes + 1 : N / num_processes;
    data_t *data = (data_t*)malloc(chunk_size*sizeof(data_t));
    long int seed = rank;
    #if defined(_OPENMP)
    #pragma omp parallel
    {
        int me             = omp_get_thread_num();
        //short int seed     = time(NULL) % ( (1 << sizeof(short int))-1 );
        short int seeds[3] = {seed-me, seed+me, seed+me*2};

    #pragma omp for
        for ( int i = 0; i < N; i++ )
        data[i].data[HOT] = erand48( seeds );
    }
    #else
    {
        //seed = time(NULL);
        srand48(seed);
        
        PRINTF("ssed is % ld\n", seed);
        
        for ( int i = 0; i < N; i++ )
        data[i].data[HOT] = drand48();
    }    
    #endif
    
    // if (rank == 0){
    //     // printf("Generating array of size %d\n", N);
    //     // printf("Array before sorting:\n");
    //     // show_array(data, 0, N, 0);

    //     // int pivot = partition(data, 0, N, compare_ge);
    //     // printf("Pivot is %d\n", pivot);

    //     // void* pivot_ptr = (void*)&data[2];
    //     // float pivot_val = data[2].data[HOT];
    //     // int pivot_pos = mpi_partition(data, 0, N, compare_ge, pivot_ptr);
    //     // printf("Pivot position wrt %f is %d\n", pivot_val, pivot_pos);
    // }


    // ---------------------------------------------
    // Create custom MPI data type for data_t

    // ---------------------------------------------
    // Divide the array into chunks and scatter them to the processes
    // Generate local array to store the scattered data of the right size
    // int chunk_size = (rank < N % num_processes) ? N / num_processes + 1 : N / num_processes;
    // int* chunk_sizes = (int*)malloc(num_processes*sizeof(int));
    // for (int i = 0; i < num_processes; i++){
    //     chunk_sizes[i] = (i < N % num_processes) ? N / num_processes + 1 : N / num_processes;
    // }
    printf("Process %d has chunk size %d\n", rank, chunk_size);
    // data_t* loc_data = (data_t*)malloc(chunk_size*sizeof(data_t));
    // data_t* loc_data = (data_t*)malloc(chunk_sizes[rank]*sizeof(data_t));
    // printf("Process %d has allocated %d bytes\n", rank, chunk_size*sizeof(data_t));
    // printf("Process %d has allocated %d bytes\n", rank, chunk_sizes[rank]*sizeof(data_t));
    //divide(data, 0, N, MPI_DATA_T, loc_data, num_processes);

    //free(data);

    // Print scattered arrays
    for (int i = 0; i < num_processes; i++){
        if (rank == i){
            printf("Process %d received:\n", rank);
            show_array(data, 0, chunk_size, 0);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // data_t* try = (data_t*)malloc(0*sizeof(data_t));


    // ---------------------------------------------
    // try on mpi_qicksort
    // mpi_quicksort(&loc_data, &chunk_size, 0, num_processes - 1, rank, MPI_DATA_T, compare_ge);
    // if (rank <= 1)
    // mpi_quicksort(&loc_data, &chunk_size, 0, 1, rank, MPI_DATA_T, compare_ge);
    // else
    // mpi_quicksort(&loc_data, &chunk_size, 2, 3, rank, MPI_DATA_T, compare_ge);
    double start_time = MPI_Wtime();

    mpi_quicksort1(&data, &chunk_size, MPI_DATA_T, MPI_COMM_WORLD);
    
    double end_time = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);

    // ---------------------------------------------
    // Show the sorted array
    for (int i = 0; i < num_processes; i++){
        if (rank == i){
            printf("Process %d has sorted:\n", rank);
            show_array(data, 0, chunk_size, 0);
        }
        //MPI_Barrier(MPI_COMM_WORLD);
    }

    // ---------------------------------------------
    // Run the parallel quicksort algorithm
    // data_t* sorted = run(data, 0, N, MPI_DATA_T, num_processes, rank);
    // if (rank == 0){
    //     printf("Array after sorting:\n");
    //     show_array(sorted, 0, N, 0);
    // }

    // for (int i = 0; i < num_processes; i++){
    //     if (rank == i){
    //         printf("Process %d has finished\n", rank);
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    // // ---------------------------------------------
    // // Gather the local arrays into a single array
    // printf("Process %d has arrived to gather\n", rank);
    // data_t* merged = (data_t*)malloc(N*sizeof(data_t));
    // gather(merged, 0, N, MPI_DATA_T, loc_data, chunk_size, num_processes);
    // printf("Process %d has gathered\n", rank);

    // if (rank == 0){
    //     printf("Array after sorting:\n");
    //     show_array(merged, 0, N, 0);
    // }
    double time = end_time - start_time;
    // ---------------------------------------------
    // Verify the results
    int test = verify_global_sorting(data, 0, chunk_size, MPI_DATA_T, rank, num_processes, 0);

    // int global_test = 0;
    // MPI_Reduce(&test, &global_test, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


    if (rank == 0){
        //if (global_test == num_processes)
        printf("Test result is %d\n", global_test);
        printf("Time elapsed is %f\n", time);
        // else
        // printf("Test failed\n");
    }

    // printf("Rank %d has loc_data of size %d\n", rank, chunk_size);
    //free(data);
    free(data);
    // free(sorted);
    // free(sorted);
    // free(merged);
    MPI_Type_free(&MPI_DATA_T);
    printf("Data freed\n");
    //MPI_Finalize();
    int finalize_retcode = MPI_Finalize();
    fprintf(stderr, "Process, return_code\n");
    fprintf(stderr, "%i, %i\n", rank, finalize_retcode);
    // free(try);
    // free(chunk_sizes);

    return 0;
}

#define SWAP(A,B,SIZE) do {int sz = (SIZE); char *a = (A); char *b = (B); \
do { char _temp = *a;*a++ = *b;*b++ = _temp;} while (--sz);} while (0)

int partition(data_t* data, int start, int end, compare_t cmp_ge){

    // Pick the median of the [0], [mid] and [end] element as pivot
    int mid = (start + end-1) / 2;
    if (cmp_ge((void*)&data[start], (void*)&data[mid]))
        SWAP((void*)&data[start], (void*)&data[mid], sizeof(data_t));
    if (cmp_ge((void*)&data[mid], (void*)&data[end-1]))
        SWAP((void*)&data[mid], (void*)&data[end-1], sizeof(data_t));
    if (cmp_ge((void*)&data[mid], (void*)&data[start]))
        SWAP((void*)&data[start], (void*)&data[mid], sizeof(data_t));

    // Pick the first element as pivot
    void* pivot = (void*)&data[start];

    // Partition around the pivot
    int pointbreak = start + 1;

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

        #pragma omp taskwait
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

void gather(data_t* data, int start, int end, MPI_Datatype MPI_DATA_T, data_t* loc_data, int chunk_size, int num_procs){
    // Function that gathers the local arrays of the processes
    // into a single array

    // int size = end - start;4
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int* sendcounts = (int*)malloc(num_procs*sizeof(int));
    int* displs = (int*)malloc(num_procs*sizeof(int));

    // Calculate sendcounts and displs
    int displacement = -1;
    for (int i = 0; i < num_procs; i++) {
        sendcounts[i] = chunk_size;
        displs[i] = displacement;
        displacement += chunk_size;
    }

    MPI_Gatherv(loc_data, sendcounts[0], MPI_DATA_T, data, sendcounts, displs, MPI_DATA_T, 0, MPI_COMM_WORLD);

    free(sendcounts);
    free(displs);
}


int mpi_partition(data_t* data, int start, int end, compare_t cmp_ge, void* pivot){
    // Function that partitions the array into two parts given a pivot
    // and returns the position of the last element of the first part

    // Partition around the pivot
    int pointbreak = start;

    // This can't be done in parallel because of possible data races in the exchanges and pointbreak increment
    // Could be done with a parallel for loop, synchronized by an atomic increment of pointbreak, but it would be slower
    for (int i = start; i < end; ++i){
        if (!cmp_ge((void*)&data[i], pivot)){
            
            // Move elements less than pivot to the left side
            SWAP((void*)&data[i], (void*)&data[pointbreak], sizeof(data_t));

            ++ pointbreak;
            
        }
    }

    // // Put the pivot in the right place
    // SWAP((void*)&data[start], (void*)&data[pointbreak - 1], sizeof(data_t));

    // Return the pivot position
    return pointbreak - 1;
}


void mpi_quicksort1 (data_t** loc_data, int* chunk_size, MPI_Datatype MPI_DATA_T, MPI_Comm comm){
    int rank, num_procs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_procs);
    
    if (num_procs > 1){
        int pivot_rank = (num_procs - 1) / 2;
        printf("Pivot rank is %d\n", pivot_rank);
        data_t* pivot = (data_t*)malloc(sizeof(data_t));
        data_t* pivots = (data_t*)malloc((num_procs+1)*sizeof(data_t));
        // MPI_Gather(&(*loc_data)[(*chunk_size -1)/2], 1, MPI_DATA_T, pivots, 1, MPI_DATA_T, 0, comm);
        // Generate a random index within each chunk
        srand(time(NULL));
        int random_index = rand() % *chunk_size;

        // Select the random element from the local data
        data_t local_pivot;
        memcpy(&local_pivot, &(*loc_data)[random_index], sizeof(data_t));

        // Gather the randomly selected elements from all processes
        MPI_Gather(&local_pivot, 1, MPI_DATA_T, pivots, 1, MPI_DATA_T, 0, comm);
        if (rank == 0){
            par_quicksort(pivots, 0, num_procs, compare_ge);
            memcpy(pivot, &pivots[(num_procs / 2)], sizeof(data_t));
        }
        MPI_Bcast(pivot, 1, MPI_DATA_T, 0, comm);
        (void*)pivot;
        int pivot_pos = mpi_partition(*loc_data, 0, *chunk_size, compare_ge, pivot);
        free(pivot);

        MPI_Comm left_comm, right_comm;
        MPI_Comm_split(comm, rank <= pivot_rank, rank, &left_comm);

        // To manage the case where the number of processes is odd
        (num_procs % 2 == 0) ? MPI_Comm_split(comm, rank > pivot_rank, rank, &right_comm) : MPI_Comm_split(comm, rank >= pivot_rank, rank, &right_comm);

        // To avoid deadlock, the pivot process will call the function recursively
        if (num_procs % 2 != 0 && rank == pivot_rank){
            mpi_quicksort1(loc_data, chunk_size, MPI_DATA_T, left_comm);
            mpi_quicksort1(loc_data, chunk_size, MPI_DATA_T, right_comm);
        }

        if (rank < pivot_rank || (num_procs % 2 == 0 && rank == pivot_rank)){
            int elements_to_send = *chunk_size - (pivot_pos + 1);
            MPI_Send(&elements_to_send, 1, MPI_INT, rank + pivot_rank + 1, 0, comm);
            int recv_elements;
            MPI_Recv(&recv_elements, 1, MPI_INT, rank + pivot_rank + 1, 0, comm, MPI_STATUS_IGNORE);
            recv_elements = (recv_elements > 0) ? recv_elements : 0;
            int new_chunk_size = *chunk_size - elements_to_send + recv_elements;
            data_t* merged = (data_t*)malloc((new_chunk_size)*sizeof(data_t));
            for (int i = 0; i <= pivot_pos; i++){
                merged[i] = (*loc_data)[i];
            }
            MPI_Send(&(*loc_data)[pivot_pos + 1], elements_to_send, MPI_DATA_T, rank + pivot_rank + 1, 0, comm);
            free(*loc_data);
            MPI_Recv(&merged[pivot_pos + 1], recv_elements, MPI_DATA_T, rank + pivot_rank + 1, 0, comm, MPI_STATUS_IGNORE);
            *chunk_size = new_chunk_size;
            mpi_quicksort1(&merged, chunk_size, MPI_DATA_T, left_comm);
            *loc_data = merged;
        }
        if (rank > pivot_rank){
            int recv_elements;
            MPI_Recv(&recv_elements, 1, MPI_INT, rank - (pivot_rank + 1), 0, comm, MPI_STATUS_IGNORE);
            int elements_to_send = pivot_pos +1;
            MPI_Send(&elements_to_send, 1, MPI_INT, rank - (pivot_rank +1), 0, comm);
            recv_elements = (recv_elements > 0) ? recv_elements : 0;
            int new_chunk_size = *chunk_size - elements_to_send + recv_elements;
            data_t* merged = (data_t*)malloc((new_chunk_size)*sizeof(data_t));
            for (int i = pivot_pos + 1; i < *chunk_size; i++){
                merged[i - (pivot_pos + 1) + recv_elements] = (*loc_data)[i];
            }
            MPI_Recv(&merged[0], recv_elements, MPI_DATA_T, rank - (pivot_rank +1), 0, comm, MPI_STATUS_IGNORE);
            MPI_Send(&(*loc_data)[0], pivot_pos +1, MPI_DATA_T, rank - (pivot_rank +1), 0, comm);
            free(*loc_data);
            *chunk_size = new_chunk_size;
            mpi_quicksort1(&merged, chunk_size, MPI_DATA_T, right_comm);
            *loc_data = merged;
        }
        MPI_Comm_free(&left_comm);
        MPI_Comm_free(&right_comm);
    } else {
        #pragma omp parallel
        {
            #pragma omp single
            par_quicksort(*loc_data, 0, *chunk_size, compare_ge);
        }
    }
}

// Single chunk version
int verify_sorting( data_t *data, int start, int end, int not_used )
{
    int i = start;
    while( (i <= end) && (data[i].data[HOT] >= data[i-1].data[HOT]) )
        i++;
    return ( i == end );
}

// Global version
int verify_global_sorting( data_t *loc_data, int start, int end, MPI_Datatype MPI_DATA_T, int rank, int num_procs, int not_used )
{
    // First I check that the local array is sorted
    verify_sorting( loc_data, start, end, not_used );

    // Then I check that the last element of the local array is less than or equal to the first element of the next process
    if (rank >= 0 && rank < num_procs - 1) {
        // Send the last element of loc_data to the next process (rank + 1)
        if (end - start > 0)
            MPI_Send(&loc_data[end - start - 1], 1, MPI_DATA_T, rank + 1, 0, MPI_COMM_WORLD);
        else{
            // If the array is empty, send a dummy element
            data_t dummy;
            // To access the index's element of the dummy array
            // dummy[index].data[HOT]
            dummy.data[HOT] = -1;
            MPI_Send(&dummy, 1, MPI_DATA_T, rank + 1, 0, MPI_COMM_WORLD);
        }
    }

    if (rank >0 && rank <= num_procs - 1) {
        data_t prev_last;
        // Receive the last element from the previous process (rank - 1)
        MPI_Recv(&prev_last, 1, MPI_DATA_T, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Check if the first element of the current process is greater than or equal to the last element of the previous process
        if (loc_data[0].data[HOT] > prev_last.data[HOT]) {
            // If not sorted, return 0
            return 0;
        }
    }

    // If everything is fine, return 1 to indicate global sorting
    return 1;

    // A LIVELLO DI SINGOLO PROCESSO FUNZIONA, MA BISOGNEREBBE TROVARE IL MODO DI ACCUMULARE I CODICI DI USCITA SU UN PROCESSO
    // (ad esempio il master) E POI FARE UNA VERIFICA GLOBALE (IDEA: USARE UNA MPI_REDUCE E VEDERE SE LA SOMMA DEI CODICI DI USCITA 
    // E' UGUALE AL NUMERO DI PROCESSI).

    // OGNUNO AVRA IL SUO CODICE DI USCITA, CHE SARA' 1 SE TUTTO E' ANDATO BENE, 0 ALTRIMENTI. POI RANK0 FARA' UNA REDUCE SUI CODICI
    // DI USCITA E SE LA SOMMA E' UGUALE AL NUMERO DI PROCESSI, ALLORA TUTTO E' ANDATO BENE.
}


    

int run(data_t* data, int start, int end, MPI_Datatype MPI_DATA_T, int num_procs, int rank){
    // Function that runs the parallel quicksort algorithm
    // on the given array

    // Create custom MPI data type for data_t
    // MPI_Type_contiguous(sizeof(data_t), MPI_BYTE, &MPI_DATA_T);
    // MPI_Type_commit(&MPI_DATA_T);
     int N = end - start;
     int num_processes = num_procs;

    // Divide the array into chunks and scatter them to the processes
    int chunk_size = (rank < N % num_processes) ? N / num_processes + 1 : N / num_processes;
    data_t* loc_data = (data_t*)malloc((chunk_size)*sizeof(data_t));
    divide(data, start, end, MPI_DATA_T, loc_data, num_procs);

    // Run the parallel quicksort algorithm
    mpi_quicksort(&loc_data, &chunk_size, 0, num_procs - 1, rank, MPI_DATA_T, compare_ge);

    // Verify sorting in a distributed manner
    int sorted = verify_global_sorting(loc_data, 0, chunk_size, MPI_DATA_T, rank, num_processes, 0);


    // Free the memory
    free(data);
    free(loc_data);

    return sorted;
}