/**
 * k_means_sequential.c
 * -----------------------------------------------------------------------------
 * Sequential implementation of Lloyd's k-means for 2D integer points.
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
 *        a. ASSIGNMENT: assign each point to nearest centroid.
 *        b. ACCUMULATION: sum coordinates per cluster.
 *        c. UPDATE: compute new centroid means.
 *        d. MOVEMENT CHECK: compute average Euclidean movement;
 *           stop early if movement < EPSILON.
 *
 * COMPLEXITY:
 *   Time: O(T * N * K) distance ops, T <= MAX_ITERATIONS.
 *   Memory: O(N) points + O(N) assignments + O(K) centroids.
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
#include <time.h>

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

// Sequential version: no threading controls

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
 * Sequential Lloyd's k-means with early stopping (average centroid movement).
 *
 * PARAMETERS:
 *   points      : array of N points.
 *   num_points  : number of points (N).
 *   centroids   : array[K] (input initial positions, output final positions).
 *   assignments : array[N] (written each iteration with cluster id per point).
 *
 * COMPLEXITY:
 *   Time: O(T * N * K), Space: O(N + K).
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

        // ===================== ASSIGNMENT STEP =====================
        // For each point, scan all centroids (K small relative to N). Writes to assignments[i].
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

        // ===================== UPDATE STEP (ALLOC BUFFERS) =====================
        // Allocate accumulation arrays (zeroed).
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

        // ===================== ACCUMULATION STEP =====================
        // Sum points per assigned cluster.
        for (int i = 0; i < num_points; i++) {
            int cid = assignments[i];
            new_centroids[cid].x += points[i].x;
            new_centroids[cid].y += points[i].y;
            counts[cid]++;
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
        
        // ===================== MOVEMENT CALC =====================
        // Sum Euclidean movement across centroids to gauge convergence.
        double total_movement = 0.0;
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
int main(void)
{
    clock_t start, end; // CPU time

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

    // Start clock
    start = clock();

    // Run clustering (MAX_ITERATIONS internally controls loop count).
    k_means(points, num_points, centroids, assignments);
    
    // Stop clock
    end = clock();
    printf("K-means clustering completed in %f ms using 1 thread (sequential).\n", (end - start) * 1000.0/CLOCKS_PER_SEC);
    
    // Cleanup
    free(assignments);
    free(points);
    return 0;
}