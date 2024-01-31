#include "par_quicksort.h"

int main(int argc, char** argv){

    // Default values
    int N = N_dflt;
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


    // ---------------------------------------------
    // Initialize MPI
    int num_processes, rank;
    int mpi_err = MPI_Init(&argc, &argv);

    if (mpi_err != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, mpi_err);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ---------------------------------------------
    // Create custom MPI data type for data_t
    MPI_Datatype MPI_DATA_T;
    MPI_Type_contiguous(sizeof(data_t), MPI_BYTE, &MPI_DATA_T);
    MPI_Type_commit(&MPI_DATA_T);

    // ---------------------------------------------
    // Print the unsorted array
    if (rank == 0) {
        printf("Unsorted array:\n");
        show_array(data, 0, N, 0);
    }

    // ---------------------------------------------
    // Sinchronize all processes
    MPI_Barrier(MPI_COMM_WORLD);

    // ---------------------------------------------
    // Sort the array

    // Broadcast the dimension of the array to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Compute chunk size to be sent to each process
    // int chunk_size = (N % num_processes == 0) ? N / num_processes : N / (num_processes - 1);
    int chunk_size = (rank < N % num_processes) ? (N / num_processes) + 1 : N / num_processes;
    if (rank == 0) {
        printf("Chunk size: %d\n", chunk_size);
    }
    // Allocate memory for the chunk
    data_t* chunk = (data_t*)malloc(chunk_size*sizeof(data_t));

    // Compute sendcounts and displacements of each chunk
    int* sendcounts = (int*)malloc(num_processes*sizeof(int));
    int* displs = (int*)malloc(num_processes*sizeof(int));
    displs[0] = 0;
    for (int i=0; i<num_processes; i++) {
        sendcounts[i] = (i < N % num_processes) ? chunk_size+1 : chunk_size;
        if (i > 0) {
            displs[i] = displs[i-1] + sendcounts[i-1];
        }
    }

    MPI_Scatterv(data, sendcounts, displs, MPI_DATA_T, chunk, chunk_size, MPI_DATA_T, 0, MPI_COMM_WORLD);

    // MPI_Scatter(data, chunk_size, MPI_DATA_T, chunk, chunk_size, MPI_DATA_T, 0, MPI_COMM_WORLD);

    // Send the remaining elements to the first processes (if any)
    if (rank == 0) {
        int remaining = N % num_processes;
        if (remaining > 0) {
            #pragma omp parallel for
            for (int i=0; i<remaining; i++) {
                MPI_Send(&data[chunk_size*num_processes+i], 1, MPI_DATA_T, i, 0, MPI_COMM_WORLD);
            }
        }
    }
    if (rank < N % num_processes) {
        MPI_Recv(&chunk[chunk_size], 1, MPI_DATA_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    // Free the memory of the original array
    free(data);

    // Print the part of the array arrived at each process
    for (int i=0; i<num_processes; i++) {
        if (rank == i) {
            printf("Process %d received:\n", rank);
            show_array(chunk, 0, chunk_size, 0);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // ---------------------------------------------
    // Compute size of own chunk and
    // then sort them
    // using quick sort

    if (rank == 0)
        printf("Chunk size: %d\n", chunk_size);

    
    int own_chunk_size = (rank < N % num_processes) ? chunk_size+1 : chunk_size;

    for (int i=0; i<num_processes; i++) {
        if (rank == i) {
            printf("Process %d has %d elements\n", rank, own_chunk_size);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    // Sorting array with quick sort for every
    // chunk as called by process
    par_quicksort(chunk, 0, chunk_size, compare_ge);

    // ---------------------------------------------
    // Debug message
    for (int i=0; i<num_processes; i++) {
        if (rank == i) {
            printf("Process %d sorted:\n", rank);
            show_array(chunk, 0, chunk_size, 0);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    // ---------------------------------------------
    // Gather all the chunks to the root process
    // and merge them
    data_t* sorted_data = NULL;
    if (rank == 0) {
        sorted_data = (data_t*)malloc(N*sizeof(data_t));
    }

    // For gatherv we need to know the size of each chunk
    int* chunk_sizes = NULL;
    if (rank == 0) {
        chunk_sizes = (int*)malloc(num_processes*sizeof(int));
    }
    MPI_Gather(&own_chunk_size, 1, MPI_INT, chunk_sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // For gatherv we need to know the displacements of each chunk
    displs = NULL;
    if (rank == 0) {
        displs = (int*)malloc(num_processes*sizeof(int));
        displs[0] = 0;
        for (int i=1; i<num_processes; i++) {
            displs[i] = displs[i-1] + chunk_sizes[i-1];
        }
    }

    // Gather all the chunks to the root process
    MPI_Gatherv(chunk, own_chunk_size, MPI_DATA_T, sorted_data, chunk_sizes, displs, MPI_DATA_T, 0, MPI_COMM_WORLD);

    // Debug message
    if (rank == 0) {
        printf("Sorted array:\n");
        show_array(sorted_data, 0, N, 0);
    }

    MPI_Finalize();
    MPI_Type_free(&MPI_DATA_T);
    free(chunk);
    chunk = NULL;
    free(chunk_sizes);
    chunk_sizes = NULL;
    free(displs);
    displs = NULL;
    free(sorted_data);
    sorted_data = NULL;

    return 0;
}