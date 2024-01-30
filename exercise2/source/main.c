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


    char* env_var = getenv("OMP_NUM_THREADS");
    if (env_var != NULL) {
        int nthreads = atoi(env_var);
    } else {
        printf("OMP_NUM_THREADS environment variable not set.\n");
    }

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
    int chunk_size = (N % num_processes == 0) ? N / num_processes : N / (num_processes - 1);

    // Allocate memory for the chunk
    data_t* chunk = (data_t*)malloc(chunk_size*sizeof(data_t));

    // Define arrays to describe scatter
    int *sendcounts = NULL;
    int *displs = NULL;

    // Allocate memory for sendcounts and displs
    if (rank == 0) {
        sendcounts = (int*)malloc(num_processes * sizeof(int));
        displs = (int*)malloc(num_processes * sizeof(int));
        
        // Compute sendcounts and displs
        int elements_per_process = (N / num_processes)+1;
        printf("elements_per_process: %d\n", elements_per_process);
        int remaining_elements = N % num_processes;
        printf("remaining_elements: %d\n", remaining_elements);
        int current_displ = 0;
        for (int i = 0; i < num_processes; i++) {
            sendcounts[i] = elements_per_process + (i <= (remaining_elements-1) ? 0 : -1);
            displs[i] = current_displ;
            current_displ += sendcounts[i];
        }
    }

    // Scatter the array to all processes
    MPI_Scatterv(data, sendcounts, displs, MPI_DATA_T, chunk, chunk_size, MPI_DATA_T, 0, MPI_COMM_WORLD);

    // Free memory for sendcounts and displs
    if (rank == 0) {
        free(sendcounts);
        free(displs);
    }


    // Free the memory of the original array
    free(data);
    data = NULL;

    // Print the part of the array arrived at each process
    for (int i=0; i<num_processes; i++) {
        if (rank == i) {
            printf("Process %d received:\n", rank);
            show_array(chunk, 0, chunk_size, 0);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    // // Compute the start and end index of the chunk
    // int start = rank * chunk_size;
    // int end = start + chunk_size;

    // own_chunk_size = (number_of_elements
    //                   >= chunk_size * (rank_of_process + 1))
    //                      ? chunk_size
    //                      : (number_of_elements
    //                         - chunk_size * rank_of_process);

    // // Sort the chunk
    // par_quicksort(chunk, 0, own_chunk_size, compare_ge);

    MPI_Finalize();
    MPI_Type_free(&MPI_DATA_T);
    free(chunk);
    chunk = NULL;

    return 0;
}