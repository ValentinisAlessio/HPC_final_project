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

    // Print unordered array
    printf("Unordered array:\n");
    show_array(data, 0, N, 0);
    printf("\n");

    int num_processes, rank;
    char* env_var = getenv("OMP_NUM_THREADS");
    if (env_var != NULL) {
        int nthreads = atoi(env_var);
    } else {
        printf("OMP_NUM_THREADS environment variable not set.\n");
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (num_processes == 1){
        omp_set_num_threads(nthreads);
        par_quicksort(data, 0, N, compare_ge, nthreads);
    }else{
        int local_size = N/num_processes;
        int local_start = rank*local_size;
        int local_end = local_start + local_size;
        data_t *local_data = (data_t*)malloc(local_size*sizeof(data_t));
        MPI_Scatter(data, local_size, MPI_DOUBLE, local_data, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        par_quicksort(local_data, 0, local_size, compare_ge, nthreads);
        MPI_Gather(local_data, local_size, MPI_DOUBLE, data, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        free(local_data);
    }

    // Print ordered array
    printf("Ordered array:\n");
    show_array(data, 0, N, 0);
    free(data);
    MPI_Finalize();

    return 0;
}