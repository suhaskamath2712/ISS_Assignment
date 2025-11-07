/**
 * k_means_parallel.c
 * -----------------------------------------------------------------------------
 * Parallel (OpenMP) implementation of Lloyd's k-means for 2D integer points.
 * Provides both static- and dynamic-scheduling variants and records timing.
 *
 * INPUT FORMAT:
 *   CSV file (no header) with lines:  x,y\n  (both integers). Example:
 *       10,42
 *       -3,17
 *   File path configured via INPUT_FILE.
 *
 * ALGORITHM (LLOYD'S):
 *   1. Initialize K centroids.
 *   2. Repeat until convergence or MAX_ITERATIONS:
 *        a. ASSIGNMENT (parallel): assign each point to nearest centroid.
 *        b. ACCUMULATION (parallel): per-thread partial sums and counts;
 *           merge partials into global arrays once (no critical section).
 *        c. UPDATE: compute new centroid means.
 *        d. MOVEMENT CHECK (parallel reduction): compute movement across
 *           centroids; stop early when movement < EPSILON.
 *
 * PARALLELIZATION STRATEGY:
 *   - Assignment: parallel for over points. Static vs dynamic scheduling is
 *     selectable via argv[2].
 *   - Accumulation: each thread writes into its own partial arrays of size K
 *     and a single-thread serial merge combines them (O(threads * K)).
 *   - Movement: reduction over centroid indices. For efficiency, squared
 *     movement is summed and compared to K * EPSILON^2.
 *   - Thread count configurable via argv[1]; passed to omp_set_num_threads.
 *
 * EDGE CASES:
 *   - num_points < K -> abort with message.
 *   - Empty cluster -> centroid remains unchanged for that iteration.
 *
 * OUTPUT:
 *   Prints final centroids and a runtime summary, and appends a CSV line to
 *   OUTPUT_FILE: isDynamic,threads,elapsed_ms.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

// File to read data from 
#define INPUT_FILE "data_20k.csv"

//File to write results to
#define OUTPUT_FILE "output_parallel.csv"

// Number of clusters 
#define K 20

// Maximum number of iterations 
#define MAX_ITERATIONS 100

// Convergence threshold (in average centroid movement, Euclidean units).
// When avg movement < EPSILON, iterations stop early.
#define EPSILON 0.0001

#define NO_OF_POINTS 20000

/**
 * SQUARED_DISTANCE
 * ----------------------------------------
 * Squared Euclidean distance between two 2D points (avoids sqrt).
 * Usage: SQUARED_DISTANCE(a, b) where a and b are Point structs.
 */
#define SQUARED_DISTANCE(a, b) (((a).x - (b).x) * ((a).x - (b).x) + ((a).y - (b).y) * ((a).y - (b).y))

// 2D integer point
typedef struct {
    int x;
    int y;
} Point;

// Configured via argv
int number_of_threads = 1;
int use_dynamic_schedule = 0; // 0=static, 1=dynamic

/** read_csv
 *  Reads INPUT_FILE as CSV of x,y ints into a heap-allocated array.
 *  Returns number of points (>=0) or -1 on error. O(N).
 */
int read_csv(Point** points)
{
    FILE* file = fopen(INPUT_FILE, "r");

    if (!file)
    {
        perror("Could not open file");
        return -1;
    }

    int capacity = NO_OF_POINTS;       // Initial capacity for dynamic array growth.
    int size = 0;
    *points = (Point*)malloc(capacity * sizeof(Point));
    if (!*points)
    {
        fclose(file);
        fprintf(stderr, "Memory allocation failed for points array\n");
        return -1;
    }

    // Read until EOF
    while (fscanf(file, "%d,%d", &(*points)[size].x, &(*points)[size].y) == 2)
    {
        size++;
        if (size >= capacity) {
            capacity *= 2; // If more space needed, double capacity
            Point* tmp = (Point*)realloc(*points, capacity * sizeof(Point));
            if (!tmp) {
                fprintf(stderr, "Realloc failed when expanding to %d points\n", capacity);
                free(*points);
                fclose(file);
                return -1;
            }
            *points = tmp;
        }
    }

    fclose(file);
    return size;
}

// Append a results row to OUTPUT_FILE: isDynamic,threads,elapsed_ms
int write_csv(int isDynamic, int threadCount, double elapsed_ms)
{
    FILE* file = fopen(OUTPUT_FILE, "a");
    if (!file) {
        perror("Could not open output file");
        return -1;
    }

    fprintf(file, "%d,%d,%f\n", isDynamic, threadCount, elapsed_ms);
    fclose(file);
    return 0;
}

/**
 * k_means_static
 * -----------------------------------------------------------------------------
 * Static-scheduling variant of Lloyd's k-means with early stopping.
 * Phases: assignment (parallel-for static), accumulation (per-thread partials
 * merged once serially), update, and movement check (parallel reduction of
 * squared movement vs K*EPSILON^2).
 */
void k_means_static(Point* points, int num_points, Point* centroids, int* assignments)
{
    // Keep previous centroids for movement computation
    Point* old_centroids = (Point*)malloc(K * sizeof(Point));
    if (!old_centroids) {
        fprintf(stderr, "Failed to allocate old_centroids\n");
        return; // Abort clustering
    }

    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
    // Snapshot current centroids
        memcpy(old_centroids, centroids, K * sizeof(Point));

    // Assignment: nearest centroid per point (static scheduling)
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < num_points; i++) {
            int min_distance = SQUARED_DISTANCE(points[i], centroids[0]);
            int best = 0;
            for (int j = 1; j < K; j++) {
                int distance = SQUARED_DISTANCE(points[i], centroids[j]);
                if (distance < min_distance) {
                    min_distance = distance;
                    best = j;
                }
            }
            assignments[i] = best;
        }

    // Allocate global accumulators and per-thread partials
    Point* new_centroids = (Point*)calloc(K, sizeof(Point));
    int* counts = (int*)calloc(K, sizeof(int));

        if (!new_centroids || !counts)
        {
            fprintf(stderr, "Allocation failed during k-means update step\n");
            free(new_centroids);
            free(counts);
            // We must also free old_centroids before returning
            free(old_centroids); // Ensure all memory is freed on error
            return; 
        }

    // Accumulate into per-thread partials (no critical section)
        int P = omp_get_max_threads();
        Point* partial_sums = (Point*)calloc((size_t)P * K, sizeof(Point));
        int* partial_counts = (int*)calloc((size_t)P * K, sizeof(int));
        if (!partial_sums || !partial_counts) {
            fprintf(stderr, "Allocation failed for partial arrays\n");
            free(partial_sums); free(partial_counts);
            free(new_centroids); free(counts); free(old_centroids);
            return;
        }

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            Point* ls = &partial_sums[(size_t)tid * K];
            int* lc = &partial_counts[(size_t)tid * K];

            #pragma omp for schedule(static)
            for (int i = 0; i < num_points; i++) {
                int cid = assignments[i];
                ls[cid].x += points[i].x;
                ls[cid].y += points[i].y;
                lc[cid]++;
            }
        }

    // Merge per-thread partials serially
        for (int c = 0; c < K; ++c) {
            long sx = 0, sy = 0; long cnt = 0;
            for (int t = 0; t < P; ++t) {
                size_t idx = (size_t)t * K + c;
                sx  += partial_sums[idx].x;
                sy  += partial_sums[idx].y;
                cnt += partial_counts[idx];
            }
            new_centroids[c].x = (int)sx;
            new_centroids[c].y = (int)sy;
            counts[c] = (int)cnt;
        }
        free(partial_sums);
        free(partial_counts);

    // Update centroids (means); keep previous if cluster empty
        for (int j = 0; j < K; j++)
        {
            if (counts[j] > 0) {
                centroids[j].x = new_centroids[j].x / counts[j];
                centroids[j].y = new_centroids[j].y / counts[j];
            }
        }
        
    // Movement: parallel reduction on squared displacement, compare to K*EPS^2
        double total_sq = 0.0;
        #pragma omp parallel for reduction(+:total_sq) schedule(static)
        for (int j = 0; j < K; j++) {
            double dx = (double)centroids[j].x - (double)old_centroids[j].x;
            double dy = (double)centroids[j].y - (double)old_centroids[j].y;
            total_sq += dx*dx + dy*dy;
        }
        
    // Early stopping threshold reached
        if (total_sq < (double)K * EPSILON * EPSILON) {
            free(new_centroids);
            free(counts);
            break; // Exit the loop early
        }

        free(new_centroids);
        free(counts);
    }
    
    // Cleanup
    free(old_centroids);
}

/**
 * k_means_dynamic
 * -----------------------------------------------------------------------------
 * Dynamic-scheduling variant (uses schedule(dynamic)) for load balancing.
 * Same accumulation, update, and movement phases as the static variant.
 */
void k_means_dynamic(Point* points, int num_points, Point* centroids, int* assignments)
{
    // Keep previous centroids for movement computation
    Point* old_centroids = (Point*)malloc(K * sizeof(Point));
    if (!old_centroids) {
        fprintf(stderr, "Failed to allocate old_centroids\n");
        return; // Abort clustering
    }

    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
    // Snapshot current centroids
        memcpy(old_centroids, centroids, K * sizeof(Point));

    // Assignment: nearest centroid per point (dynamic scheduling)
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < num_points; i++) {
            int min_distance = SQUARED_DISTANCE(points[i], centroids[0]);
            int best = 0;
            for (int j = 1; j < K; j++) {
                int distance = SQUARED_DISTANCE(points[i], centroids[j]);
                if (distance < min_distance) {
                    min_distance = distance;
                    best = j;
                }
            }
            assignments[i] = best;
        }

    // Allocate global accumulators and per-thread partials
    Point* new_centroids = (Point*)calloc(K, sizeof(Point));
    int* counts = (int*)calloc(K, sizeof(int));

        if (!new_centroids || !counts)
        {
            fprintf(stderr, "Allocation failed during k-means update step\n");
            free(new_centroids);
            free(counts);
            // We must also free old_centroids before returning
            free(old_centroids); // Ensure all memory is freed on error
            return; 
        }

    // Accumulate into per-thread partials (no critical section)
        int P = omp_get_max_threads();
        Point* partial_sums = (Point*)calloc((size_t)P * K, sizeof(Point));
        int* partial_counts = (int*)calloc((size_t)P * K, sizeof(int));
        if (!partial_sums || !partial_counts) {
            fprintf(stderr, "Allocation failed for partial arrays\n");
            free(partial_sums); free(partial_counts);
            free(new_centroids); free(counts); free(old_centroids);
            return;
        }

        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            Point* ls = &partial_sums[(size_t)tid * K];
            int* lc = &partial_counts[(size_t)tid * K];

            #pragma omp for schedule(dynamic)
            for (int i = 0; i < num_points; i++) {
                int cid = assignments[i];
                ls[cid].x += points[i].x;
                ls[cid].y += points[i].y;
                lc[cid]++;
            }
        }

    // Merge per-thread partials serially
        for (int c = 0; c < K; ++c) {
            long sx = 0, sy = 0; long cnt = 0;
            for (int t = 0; t < P; ++t) {
                size_t idx = (size_t)t * K + c;
                sx  += partial_sums[idx].x;
                sy  += partial_sums[idx].y;
                cnt += partial_counts[idx];
            }
            new_centroids[c].x = (int)sx;
            new_centroids[c].y = (int)sy;
            counts[c] = (int)cnt;
        }
        free(partial_sums);
        free(partial_counts);

    // Update centroids (means); keep previous if cluster empty
        for (int j = 0; j < K; j++)
        {
            if (counts[j] > 0) {
                centroids[j].x = new_centroids[j].x / counts[j];
                centroids[j].y = new_centroids[j].y / counts[j];
            }
        }
        
    // Movement: parallel reduction on squared displacement, compare to K*EPS^2
        double total_sq = 0.0;
        #pragma omp parallel for reduction(+:total_sq) schedule(dynamic)
        for (int j = 0; j < K; j++) {
            double dx = (double)centroids[j].x - (double)old_centroids[j].x;
            double dy = (double)centroids[j].y - (double)old_centroids[j].y;
            total_sq += dx*dx + dy*dy;
        }
        
    // Early stopping threshold reached
        if (total_sq < (double)K * EPSILON * EPSILON) {
            free(new_centroids);
            free(counts);
            break; // Exit the loop early
        }

        free(new_centroids);
        free(counts);
    }
    
    // Cleanup
    free(old_centroids);
}

/**
 * main
 * -----------------------------------------------------------------------------
 * Usage: program <threads> <0|1>
 *   threads: number of OpenMP threads
 *   0|1   : 0 = static scheduling, 1 = dynamic scheduling
 * Steps:
 *   1) Parse args and configure OpenMP
 *   2) Read points from CSV
 *   3) Initialize centroids (first K points)
 *   4) Run selected k-means variant and time with omp_get_wtime
 *   5) Print centroids and timing; append CSV row to OUTPUT_FILE
 * Returns: elapsed ms on success (as float), non-zero error code on failure.
 */
int main(int argc, char** argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s <number_of_threads> <0|1>\n", argv[0]);
        return 1;
    }
    
    number_of_threads = atoi(argv[1]);

    if (number_of_threads < 1)
        number_of_threads = 1;

    omp_set_num_threads(number_of_threads);

    use_dynamic_schedule = atoi(argv[2]); // second arg: 0=static, 1=dynamic

    double start, end; // wall-clock timing using OpenMP

    Point* points = NULL;
    int num_points = read_csv(&points);

    if (num_points <= 0)
    {
        fprintf(stderr, "Error reading points from CSV\n");
        return 1;
    }

    // Check if K is too large (simple safeguard)
    if (num_points < K)
    {
        fprintf(stderr, "Error: K (%d) is larger than number of points (%d).\n", K, num_points);
        free(points);
        return 1;
    }

    // Initialize centroids (for simplicity, use first K points; better: random sample or k-means++).
    Point centroids[K];

    for (int i = 0; i < K; i++)
        centroids[i] = points[i];

    int* assignments = (int*) malloc(num_points * sizeof(int));
    if (!assignments)
    {
        fprintf(stderr, "Failed to allocate assignments array\n");
        free(points);
        return 1;
    }

    //Start clock
    start = omp_get_wtime();

    if (use_dynamic_schedule)
        k_means_dynamic(points, num_points, centroids, assignments);
    else
        k_means_static(points, num_points, centroids, assignments);
    
    //Stop clock
    end = omp_get_wtime();

    printf("Final centroids:\n");

    //Print each centroid with the number of points assigned to it
    for (int i = 0; i < K; i++)
    {
        int count = 0;
        for (int j = 0; j < num_points; j++)
            if (assignments[j] == i)
                count++;

        printf("Centroid %d: (%d, %d) with %d points in the cluster\n", i, centroids[i].x, centroids[i].y, count);
    }

    float elapsed_ms = (end - start) * 1000;

    if (use_dynamic_schedule)
        printf("K-means clustering completed in %f ms using %d thread(s) with DYNAMIC scheduling.\n", elapsed_ms, number_of_threads);
    else
        printf("K-means clustering completed in %f ms using %d thread(s) with STATIC scheduling.\n", elapsed_ms, number_of_threads);

    printf("----------------------------------------------\n");

    //Write results to CSV
    write_csv(use_dynamic_schedule, number_of_threads, elapsed_ms);

    // Cleanup
    free(assignments);
    free(points);

    return 0;
}