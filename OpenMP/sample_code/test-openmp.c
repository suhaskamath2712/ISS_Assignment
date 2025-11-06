/******************************************************************************
 *  * * OpenMP Example - Loop Work-sharing - C/C++ Version
 *   * * FILE: test-openmp.c
 *    * * DESCRIPTION:
 *     * *   In this example, the iterations of a loop are scheduled dynamically
 *      * *   across the team of threads.  A thread will perform CHUNK iterations
 *       * *   at a time before being scheduled for the next CHUNK of work.
 *        * * SOURCE: Blaise Barney  5/99
 *         * * LAST REVISED: 01/09/04
 *          * ******************************************************************************/

#include <stdio.h>
#include <omp.h>
#define CHUNKSIZE   10
#define N       100

int main ()  {

    int nthreads, tid, i, chunk;
    float a[N], b[N], c[N];

    /* Some initializations */
    for (i=0; i < N; i++)
        a[i] = b[i] = i * 1.0;
    chunk = CHUNKSIZE;

    omp_set_num_threads(16);

    #pragma omp parallel shared(a,b,c,nthreads,chunk) private(i,tid) 
        {
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }
        printf("Thread %d starting...\n",tid);

        #pragma omp for schedule(static,chunk)
        for (i=0; i<N; i++)
        {
            c[i] = a[i] + b[i];
            printf("Thread %d: c[%d]= %f\n",tid,i,c[i]);
        }

        }  /* end of parallel section */
    return 0;
}
