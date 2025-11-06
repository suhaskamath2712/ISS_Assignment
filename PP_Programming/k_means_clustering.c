/**
 * k_means_clustering.c
 * ----------------------------------------
 * Minimal C implementation of k-means clustering for 2D integer points.
 *
 * INPUT DATA SOURCE:
 * A CSV file (no header) where each line contains two comma-separated integers:
 * x,y\n
 * Example:
 * 10,42
 * -3,17
 *
 * HIGH-LEVEL ALGORITHM (STANDARD LLOYD'S):
 * 1. Initialize K centroids (here: first K points; can be improved via k-means++).
 * 2. Repeat for up to MAX_ITERATIONS (early stop enabled):
 * a. ASSIGNMENT: For each point, choose the nearest centroid by squared Euclidean distance.
 * b. UPDATE: Recompute each centroid as the arithmetic mean of points assigned to it.
 * 3. Stop early if average centroid movement is about 0 (convergence threshold).
 *
 * PERFORMANCE NOTES:
 * - Time per iteration: O(N * K) distance computations.
 * - Memory: O(N) for points + O(K) for centroids + O(N) for assignments.
 * - Uses squared distance to avoid unnecessary sqrt calls (monotonic w.r.t actual distance).
 *
 * LIMITATIONS / POSSIBLE IMPROVEMENTS:
 * - Early stopping uses an average-movement threshold.
 * - No k-means++ initialization (can reduce iterations / improve clustering quality).
 * - Centroids chosen from first K points may cause degenerate clusters if early points are close.
 * - No handling for empty input or when num_points < K beyond using those points as initial centroids.
 * - Integer division during centroid update may introduce rounding bias (consider storing doubles).
 *
 * EDGE CASES CONSIDERED:
 * - File open failure.
 * - Dynamic growth of point array with realloc.
 * - Empty cluster in update step: centroid left unchanged.
 *
 * BUILDING:
 * gcc -O2 -std=c11 k_means_clustering.c -o kmeans -lm
 * (Note: Added -lm to link the math library for sqrt)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> // Added for memcpy

// File to read data from (relative path). Adjust as needed.
#define FILE_NAME "data_20k.csv"

// Number of clusters (K). Ensure K <= number of points for robust initialization.
#define K 20

// Maximum number of iterations for Lloyd's refinement phase.
#define MAX_ITERATIONS 100

// Convergence threshold (in average centroid movement, Euclidean units).
// When avg movement < EPSILON, iterations stop early.
#define EPSILON 0.01

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

/**
 * read_csv
 * ----------------------------------------
 * Dynamically reads points from a CSV file of the form: "x,y" per line.
 * Allocates (or reallocates) an array of Point and returns its size.
 *
 * PARAMETERS:
 * filename  - path to CSV file.
 * points    - out parameter; on success will point to a heap-allocated array of Point.
 *
 * RETURNS:
 * Number of points read (>=0), or -1 on file open failure.
 *
 * MEMORY MANAGEMENT:
 * Caller is responsible for freeing *points.
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
    if (!*points) {
        fclose(file);
        fprintf(stderr, "Memory allocation failed for points array\n");
        return -1;
    }

    // Read until fscanf fails (EOF or malformed line). Lines not matching pattern are skipped.
    while (fscanf(file, "%d,%d", &(*points)[size].x, &(*points)[size].y) == 2) {
        size++;
        if (size >= capacity) {
            capacity *= 2; // Doubling strategy keeps amortized reallocation cost low.
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
 * ----------------------------------------
 * Performs Lloyd's k-means clustering with early stopping (avg centroid movement).
 *
 * PARAMETERS:
 * points       - array of N points.
 * num_points   - number of points (N).
 * centroids    - array of size k; initialized on entry, updated in-place.
 * assignments  - output array of size N; after each iteration holds centroid index for each point.
 *
 * SIDE EFFECTS:
 * - Updates centroids in-place.
 * - Writes cluster assignment for each point into assignments.
 *
 * COMPLEXITY:
 * O(ITERS * N * K) distance computations, where ITERS <= MAX_ITERATIONS.
 *
 * EDGE CASES:
 * - Empty cluster: centroid remains as previous iteration value (no movement).
 */
void k_means(Point* points, int num_points, Point* centroids, int* assignments) {
    
    // Need a place to store old centroids to check for movement
    Point* old_centroids = (Point*)malloc(K * sizeof(Point));
    if (!old_centroids) {
        fprintf(stderr, "Failed to allocate old_centroids\n");
        return; // Abort clustering
    }

    for (int iter = 0; iter < MAX_ITERATIONS; iter++) {
        
        // Store the current centroid positions *before* they are updated
        memcpy(old_centroids, centroids, K * sizeof(Point));

        // ===================== ASSIGNMENT STEP =====================
        // For each point, find the nearest centroid.
        for (int i = 0; i < num_points; i++) {
            int min_distance = SQUARED_DISTANCE(points[i], centroids[0]);
            assignments[i] = 0;
            for (int j = 1; j < K; j++) {
                int distance = SQUARED_DISTANCE(points[i], centroids[j]);
                if (distance < min_distance) {
                    min_distance = distance;
                    assignments[i] = j;
                }
            }
        }

        // ===================== UPDATE STEP =====================
        // Accumulate sums for each cluster to compute new centroids.
        Point* new_centroids = (Point*)calloc(K, sizeof(Point));
        int* counts = (int*)calloc(K, sizeof(int));
        if (!new_centroids || !counts) {
            fprintf(stderr, "Allocation failed during k-means update step\n");
            free(new_centroids);
            free(counts);
            // We must also free old_centroids before returning
            free(old_centroids); // Ensure all memory is freed on error
            return; 
        }

        for (int i = 0; i < num_points; i++) {
            new_centroids[assignments[i]].x += points[i].x;
            new_centroids[assignments[i]].y += points[i].y;
            counts[assignments[i]]++;
        }

        for (int j = 0; j < K; j++) {
            if (counts[j] > 0) {
                // Integer mean; consider casting to double for more precision if needed.
                centroids[j].x = new_centroids[j].x / counts[j];
                centroids[j].y = new_centroids[j].y / counts[j];
            }
            // else: centroid remains in its 'old_centroid' position (no change)
        }
        
        // Corrected calculation for average centroid movement
        double total_movement = 0.0;
        for (int j = 0; j < K; j++)
        {
            // Compare the NEW position (centroids[j]) with the OLD (old_centroids[j])
            // Use sqrt to get actual Euclidean distance
            total_movement += sqrt(SQUARED_DISTANCE(centroids[j], old_centroids[j]));
        }
        // Divide by K *after* the sum is complete
        double avg_movement = total_movement / K;

        printf("Iteration %d, Average Centroid Movement: %.2f\n", iter + 1, avg_movement);
        
        // Early stopping: stop when average movement is near zero
        if (avg_movement < EPSILON) {
            printf("Converged at iteration %d\n", iter + 1);
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
 * ----------------------------------------
 * Driver program:
 * - Reads points from FILE_NAME.
 * - Initializes centroids with first K points.
 * - Allocates assignments array and runs k-means.
 * - Frees allocated memory.
 *
 * EXIT CODES:
 * 0 on success, 1 on failure to read data.
 *
 * NOTE:
 * Assumes num_points >= K; if not, some centroids may duplicate.
 */
int main(void) {
    Point* points = NULL;
    int num_points = read_csv(FILE_NAME, &points);

    if (num_points < 0) {
        fprintf(stderr, "Error reading points from CSV\n");
        return 1;
    }

    printf("Read %d points from %s\n", num_points, FILE_NAME);

    // Check if K is too large (simple safeguard)
    if (num_points < K) {
        fprintf(stderr, "Error: K (%d) is larger than number of points (%d).\n", K, num_points);
        free(points);
        return 1;
    }

    // Initialize centroids (for simplicity, use first K points; better: random sample or k-means++).
    Point centroids[K];

    for (int i = 0; i < K; i++) {
        centroids[i] = points[i];
    }

    int* assignments = (int*)malloc(num_points * sizeof(int));
    if (!assignments) {
        fprintf(stderr, "Failed to allocate assignments array\n");
        free(points);
        return 1;
    }

    // Run clustering (MAX_ITERATIONS internally controls loop count).
    k_means(points, num_points, centroids, assignments);

    // Optional: Print final centroid locations
    printf("\nFinal Centroid Locations:\n");
    for(int i = 0; i < K; i++) {
        printf("Cluster %d: (%d, %d)\n", i, centroids[i].x, centroids[i].y);
    }

    free(assignments);
    free(points);
    return 0;
}