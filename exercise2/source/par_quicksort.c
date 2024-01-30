#include "par_quicksort.h"

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
    ++ start;
    void* pivot = (void*)&data[start];

    // Partition around the pivot
    int pointbreak = start + 1;
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

    if (size <= 0)
        return;
    
    if (size > 2){
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
void par_quicksort(data_t* data, int start, int end, compare_t cmp_ge, int n_threads){

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

    if (end - start <= 0)
        return;
    
    // Partition the array
    int pivot = partition(data, start, end, cmp_ge);

    // Sort the left and right side
    if(n_threads > 1) {
        #pragma omp task shared(data) if(end - start > (sizeof(data)/DATA_SIZE)/n_threads)
            par_quicksort(data, start, pivot-1, cmp_ge, n_threads);
        #pragma omp task shared(data) if(end - start > (sizeof(data)/DATA_SIZE)/n_threads)
            par_quicksort(data, pivot+1, end, cmp_ge, n_threads);
    } else {
        quicksort(data, start, pivot-1, cmp_ge);
        quicksort(data, pivot+1, end, cmp_ge);
    }
    
}


// Verification functions (for debugging)
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