#include "quick_mpi.h"
// ================================================================
//  MAIN 
// ================================================================
int main(int argc, char** argv){

    //.................................................................................................
    // (1) INITIALIZATION 
    //.................................................................................................

    int N = SIZE;
    //int nthreads;
    // {
    //     int a = 0;
    //     if ( argc > ++a ) N = atoi(*(argv+a));
    // }
    if ( argc > 1 ) N = atoi(*(argv+1));
    //char* env_var = getenv("OMP_NUM_THREADS");
    //if (env_var != NULL) {
    //    int nthreads = atoi(env_var);
    //} else {
    //    printf("OMP_NUM_THREADS environment variable not set :(\n");
    //}

    int n_processes, rank;
    int mpi_err = MPI_Init(&argc, &argv);

    if (mpi_err != MPI_SUCCESS) {
        printf("Error initializing MPI program :(\n");
        MPI_Abort(MPI_COMM_WORLD, mpi_err);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Datatype MPI_DATA_T;
    MPI_Type_contiguous(sizeof(data_t), MPI_BYTE, &MPI_DATA_T);
    MPI_Type_commit(&MPI_DATA_T);

    //.................................................................................................
    //  (2) GENERATE PROCESS CHUNK OF DATA
    //.................................................................................................

    int chunk_size = (rank < N % n_processes) ? N / n_processes + 1 : N / n_processes;
    data_t *data = (data_t*)malloc(chunk_size*sizeof(data_t));
    // set a seed depending on the process rank
    long int seed = rank;
    #if defined(_OPENMP)
    #pragma omp parallel
    {
        int me             = omp_get_thread_num();
        short int seeds[3] = {seed-me, seed+me, seed+me*2};

    #pragma omp for
        for ( int i = 0; i < chunk_size; i++ )
        data[i].data[HOT] = erand48( seeds );
    }
    #else
    {
        //seed = time(NULL);
        srand48(seed);
        
        printf("ssed is % ld\n", seed);
        
        for ( int i = 0; i < chunk_size; i++ )
        data[i].data[HOT] = drand48();
    }    
    #endif
    
    MPI_Barrier(MPI_COMM_WORLD);

    // Show the sorted array
    // for (int i = 0; i < n_processes; i++){
    //      if (rank == i){
    //          printf("---------------------------------------\n");
    //          printf("Process %d has generated:\n", rank);
    //          show_array(data, 0, chunk_size, 0);
    //          printf("---------------------------------------\n");
    //      }
    //      MPI_Barrier(MPI_COMM_WORLD);
    // }

    //.................................................................................................
    // (3) SORT THE DATA AND MEASURE TIME
    //.................................................................................................

    // Wait all processes to finish generating the data
    double t_start, t_end;
    MPI_Barrier(MPI_COMM_WORLD);
    t_start= MPI_Wtime();
    mpi_quicksort(&data, &chunk_size, MPI_DATA_T, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    t_end = MPI_Wtime();
    double time = t_end - t_start;
    
    //.................................................................................................
    // (4) VERIFY THE SORTING
    //.................................................................................................

    //MPI_Barrier(MPI_COMM_WORLD);
    // for (int i = 0; i < n_processes; i++){
    //      if (rank == i){
    //          printf("---------------------------------------\n");
    //          printf("Process %d has sorted:\n", rank);
    //          show_array(data, 0, chunk_size, 0);
    //          printf("---------------------------------------\n");
    //      }
    //      MPI_Barrier(MPI_COMM_WORLD);
    // }
    //MPI_Barrier(MPI_COMM_WORLD)
    // Verify the results
    int chunk_check = verify_global_sorting(data, 0, chunk_size, MPI_DATA_T, rank, n_processes, 0);
    MPI_Barrier(MPI_COMM_WORLD);
    int global_check = 0;
    MPI_Reduce(&chunk_check, &global_check, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0){
        printf("\n");
	if (global_check == n_processes){
	    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	    printf("Array sorted correctly :)\n");
	    printf("Execution time: %f\n", time);
	    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        } else {
	    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
            printf("The array has not been sorted correctly :(\n");
       	    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        }
    }

    //.................................................................................................
    // (4) FINALIZATION 
    //.................................................................................................
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
    MPI_Type_free(&MPI_DATA_T);
    int finalize_retcode = MPI_Finalize();
    //fprintf(stderr, "Process, return_code\n");
    //fprintf(stderr, "%i, %i\n", rank, finalize_retcode);
    return 0;
}