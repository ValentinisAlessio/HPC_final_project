#include "par_quicksort.h"
#include <unistd.h>
#include <limits.h>

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
    //  generate the array
    //
    int chunk_size = N / num_processes;
    int remainder = N % num_processes;
    int my_chunk_size = (rank < remainder) ? chunk_size + 1 : chunk_size;
    
    data_t *data = (data_t*)malloc(my_chunk_size*sizeof(data_t));
    long int seed = time(NULL) + rank + getpid();
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);
    for (int i = 0; i < HOST_NAME_MAX && hostname[i] != '\0'; i++) {
        seed += (long int)hostname[i] << ((i % sizeof(long int)) * CHAR_BIT); // Include characters from hostname
    }
    #if defined(_OPENMP)
    #pragma omp parallel private(seed)
    {
        int me             = omp_get_thread_num();
        // short int seed     = time(NULL) % ( (1 << sizeof(short int))-1 );
        seed += me * rank;
        // short int seeds[3] = {seed-me, seed+me, seed+me*2};
        unsigned short int my_seed[3] = {
            (unsigned short int)(seed % ((1 << 16) - 1) + me),
            (unsigned short int)(seed % ((1 << 16) - 1) + me + 1),
            (unsigned short int)(seed % ((1 << 16) - 1) + me + 2)
        }; // Adjust the seed for each thread

    #pragma omp for
        for ( int i = 0; i < N; i++ )
        // data[i].data[HOT] = erand48( seeds );
        data[i].data[HOT] = erand48( my_seed );
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

    for (int i=0; i<num_processes; i++) {
        if (rank == i) {
            printf("Unordered generated array by rank %d:\n", rank);
            show_array(data, 0, N, 0);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // ---------------------------------------------
    free(data);
    MPI_Type_free(&MPI_DATA_T);
    MPI_Finalize();
    
    return 0;
}