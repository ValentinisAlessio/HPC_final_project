#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

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

    printf("Size of data is %ld\n", sizeof(*data)/sizeof(data_t));

    show_array(data, 0, N, 0);

    par_quicksort(data, 0, N, compare_ge);

    show_array(data, 0, N, 0);

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