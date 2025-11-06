/**
 * k_means_clustering.c
 * -----------------------------------------------------------------------------
 * Parallel (OpenMP) implementation of Lloyd's k-means for 2D integer points.
 *
 * INPUT FORMAT:
 *   CSV file (no header) with lines:  x,y\n  (both integers). Example:
 *       10,42
 *       -3,17
 *   File path configured via FILE_NAME.
 *
 * ALGORITHM (LLOYD'S):
 *   1. Initialize K centroids (first K points; can substitute k-means++).
 *   2. Repeat until convergence or MAX_ITERATIONS:
 *        a. ASSIGNMENT (parallel): assign each point to nearest centroid.
 *        b. ACCUMULATION (parallel): per-thread local sums merged into global.
 *        c. UPDATE: compute new centroid means.
 *        d. MOVEMENT CHECK (parallel reduction): compute average Euclidean movement;
 *           stop early if movement < EPSILON.
 *
 * PARALLELIZATION STRATEGY:
 *   - Accumulation: each thread uses stack-local arrays (local_sums/local_counts) to avoid
 *     contention; merged once inside a critical region (O(K) merge per thread).
 *   - Movement: reduction over centroid indices using OpenMP reduction clause.
 *   - Thread count configurable via argv[1]; validated >=1; passed to omp_set_num_threads.
 *
 * COMPLEXITY:
 *   Sequential theoretical: O(T * N * K) distance ops, T <= MAX_ITERATIONS.
 *   Parallel ideal runtime: O(T * (N/threads) * K) + merge overhead O(threads * K).
 *   Memory: O(N) points + O(N) assignments + O(K) centroids + O(K) temporary arrays per thread.
 *
 * CONVERGENCE:
 *   Uses average centroid movement threshold (EPSILON). Movement is averaged Euclidean
 *   distance across centroids after each iteration.
 *
 * EDGE CASES:
 *   - num_points < K -> program aborts early with message.
 *   - Empty cluster: centroid remains at previous position (no reinitialization strategy).
 *
 * OUTPUT:
 *   Final runtime (ms) and optionally centroid coordinates (printing suppressed if unwanted).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

// File to read data from 
#define FILE_NAME "data_20k.csv"

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
 * Macro to compute squared Euclidean distance between two 2D points.
 * Usage: SQUARED_DISTANCE(a, b) where a and b are Point structs.
 *
 * Complexity: O(1)
 */
#define SQUARED_DISTANCE(a, b) (((a).x - (b).x) * ((a).x - (b).x) + ((a).y - (b).y) * ((a).y - (b).y))

// Struct representing a 2D point with integer coordinates.
typedef struct {
    int x;
    int y;
} Point;

int number_of_threads = 1; // configurable via argv; passed to OpenMP runtime

/**
 * read_csv
 * ----------------------------------------
 *
 * PARAMETERS:
 * filename  - path to CSV file.
 * points    - out parameter; on success will point to a heap-allocated array of Point.
 *
 * RETURNS:
 * Number of points read (>=0), or -1 on file open failure.
 *
 * COMPLEXITY:
 * O(N) where N is number of lines / points.
 *
 * ERROR HANDLING:
 * - On fopen failure, prints perror and returns -1.
 */
int read_csv(const char* filename, Point** points)
{
    FILE* file = fopen(filename, "r");

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

/**
 * k_means
 * -----------------------------------------------------------------------------
 * Parallel Lloyd's k-means with early stopping (average centroid movement).
 *
 * PARAMETERS:
 *   points      : array of N points.
 *   num_points  : number of points (N).
 *   centroids   : array[K] (input initial positions, output final positions).
 *   assignments : array[N] (written each iteration with cluster id per point).
 *
 * PARALLEL PHASES:
 *   1. Assignment      (#pragma omp parallel for) per point nearest centroid.
 *   2. Accumulation    (#pragma omp parallel) local arrays + critical merge.
 *   3. Movement reduce (#pragma omp parallel for reduction(+)) total movement.
 *
 * COMPLEXITY:
 *   Time: O(T * N * K / P + T * K * P) approx (P = threads, T iterations).
 *   Space: O(N + K + K * P) due to per-thread local_sums/local_counts.
 *
 * CONVERGENCE:
 *   Stops early if average movement < EPSILON; else runs to MAX_ITERATIONS.
 *
 * RETURNS:
 *   void (centroids & assignments mutated in place).
 */
void k_means(Point* points, int num_points, Point* centroids, int* assignments)
{
    // Need a place to store old centroids to check for movement
    Point* old_centroids = (Point*)malloc(K * sizeof(Point));
    if (!old_centroids) {
        fprintf(stderr, "Failed to allocate old_centroids\n");
        return; // Abort clustering
    }

    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        // Store the current centroid positions before they are updated
        memcpy(old_centroids, centroids, K * sizeof(Point));

        // ===================== ASSIGNMENT STEP (PARALLEL) =====================
        // For each point, scan all centroids (K small relative to N). Writes to assignments[i].
        // No synchronization required (each i unique). Static scheduling for cache locality.
        #pragma omp parallel for schedule(dynamic) nowait
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

        // ===================== UPDATE STEP (ALLOC GLOBAL BUFFERS) =====================
        // Allocate global accumulation arrays (zeroed). Each thread will build local copies and merge.
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

        // ===================== ACCUMULATION STEP (PARALLEL LOCAL MERGE) =====================
        // Each thread creates stack-local arrays to avoid false sharing & contention.
        // Critical region merges O(K) per thread; acceptable for moderate K.
        // Possible improvement: custom OpenMP reduction or atomic with padding.
        #pragma omp parallel
        {
            Point local_sums[K];
            int local_counts[K];

            for (int c = 0; c < K; ++c)
            {
                local_sums[c].x = 0;
                local_sums[c].y = 0;
                local_counts[c] = 0;
            }

            #pragma omp for schedule(dynamic) nowait
            for (int i = 0; i < num_points; i++) {
                int cid = assignments[i];
                local_sums[cid].x += points[i].x;
                local_sums[cid].y += points[i].y;
                local_counts[cid]++;
            }

            // Merge local into global (critical section per thread) to ensure correctness
            #pragma omp critical
            {
                for (int c = 0; c < K; ++c) {
                    new_centroids[c].x += local_sums[c].x;
                    new_centroids[c].y += local_sums[c].y;
                    counts[c] += local_counts[c];
                }
            }
        }

        // ===================== CENTROID UPDATE =====================
        // Compute arithmetic mean; if empty cluster, leave centroid unchanged (no reinit strategy).
        for (int j = 0; j < K; j++)
        {
            if (counts[j] > 0) {
                centroids[j].x = new_centroids[j].x / counts[j];
                centroids[j].y = new_centroids[j].y / counts[j];
            }
        }
        
        // ===================== MOVEMENT CALC (PARALLEL REDUCTION) =====================
        // Sum Euclidean movement across centroids to gauge convergence.
        double total_movement = 0.0;
        #pragma omp parallel for reduction(+:total_movement) schedule(dynamic)
        for (int j = 0; j < K; j++) {
            total_movement += sqrt(SQUARED_DISTANCE(centroids[j], old_centroids[j]));
        }
        // Divide by K after the sum is complete
        double avg_movement = total_movement / K;
        
        // Early stopping: break when average movement below EPSILON threshold.
        if (avg_movement < EPSILON) {
            free(new_centroids);
            free(counts);
            break; // Exit the loop early
        }

        free(new_centroids);
        free(counts);
    }
    
    // Free the helper array allocated at the start
    free(old_centroids);
}

/**
 * main
 * -----------------------------------------------------------------------------
 *
 * STEPS:
 *   1. Parse argv for thread count; set via omp_set_num_threads.
 *   2. Read points from CSV.
 *   3. Initialize centroids from first K points (for simplicity).
 *   4. Time k_means via omp_get_wtime (wall clock).
 *   5. Print runtime summary.
 *
 * RETURNS:
 *   0 on success; 1 on I/O or allocation failure.
 */
int main(int argc, char** argv)
{
    if (argc > 1)
    {
        number_of_threads = atoi(argv[1]);
        if (number_of_threads < 1) number_of_threads = 1;
    }
    
    omp_set_num_threads(number_of_threads);

    double start, end; // wall-clock timing using OpenMP

    Point* points = NULL;
    int num_points = read_csv(FILE_NAME, &points);

    if (num_points < 0) {
        fprintf(stderr, "Error reading points from CSV\n");
        return 1;
    }

    // Check if K is too large (simple safeguard)
    if (num_points < K) {
        fprintf(stderr, "Error: K (%d) is larger than number of points (%d).\n", K, num_points);
        free(points);
        return 1;
    }

    // Initialize centroids (for simplicity, use first K points; better: random sample or k-means++).
    Point centroids[K];

    for (int i = 0; i < K; i++)
    {
        centroids[i] = points[i];
    }

    int* assignments = (int*) malloc(num_points * sizeof(int));
    if (!assignments)
    {
        fprintf(stderr, "Failed to allocate assignments array\n");
        free(points);
        return 1;
    }

    //Start clock
    start = omp_get_wtime();

    // Run clustering (MAX_ITERATIONS internally controls loop count).
    k_means(points, num_points, centroids, assignments);
    
    //Stop clock
    end = omp_get_wtime();

    printf("K-means clustering completed in %.3f ms using %d thread(s).\n", (end - start) * 1000.0, number_of_threads);
    
    // Cleanup
    free(assignments);
    free(points);
    return 0;
}