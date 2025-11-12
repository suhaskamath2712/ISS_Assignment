"""
Simple script to compute average Time_ms and Speedup per Threads 
from a CSV file and plot line graphs for both.

Usage: python analyze_kmeans.py
Outputs: 
- Prints average times and speedups.
- Saves kmeans_avg_time.png (Runtime Plot)
- Saves kmeans_avg_speedup.png (Speedup Plot)
"""

import matplotlib.pyplot as plt
import os # Import os to check if file exists

# --- Hardcoded File Paths ---
INPUT_CSV_PATH = r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\OpenMP\analyse_kmeans\output_parallel.csv"
TIME_PLOT_PATH = r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\OpenMP\analyse_kmeans\kmeans_avg_time.png"
SPEEDUP_PLOT_PATH = r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\OpenMP\analyse_kmeans\kmeans_avg_speedup.png"
# --- End Hardcoded File Paths ---


def read_data():
    """
    Reads data from CSV output (output_parallel.csv) using the hardcoded path.
    Returns two dictionaries: dynamic_data and static_data
    Format: {threads: [time1, time2, ...], ...}
    """
    static_data = {}
    dynamic_data = {}
    
    if not os.path.exists(INPUT_CSV_PATH):
        print(f"Error: Input file not found at '{INPUT_CSV_PATH}'")
        print("Please check the hardcoded INPUT_CSV_PATH variable.")
        return None, None

    with open(INPUT_CSV_PATH, "r") as f:
        for line in f:
            try:
                is_dynamic, threads, time_ms = line.strip().split(",")
                
                threads = int(threads)
                time_ms = float(time_ms)

                if is_dynamic == "1":
                    if threads not in dynamic_data:
                        dynamic_data[threads] = []
                    dynamic_data[threads].append(time_ms)
                else:
                    if threads not in static_data:
                        static_data[threads] = []
                    static_data[threads].append(time_ms)
            except ValueError:
                print(f"Skipping malformed line: {line.strip()}")

    return dynamic_data, static_data


def main():
    dynamic_data, static_data = read_data()
    
    if dynamic_data is None or static_data is None:
        return # Exit if file wasn't read

    # --- Process Dynamic Data ---
    dynamic_threads = sorted(dynamic_data.keys())
    dynamic_avgs = [sum(dynamic_data[t]) / len(dynamic_data[t]) for t in dynamic_threads]

    print("--- Dynamic Schedule Averages ---")
    print("Threads, AvgTime_ms, NumRuns")
    for t, avg in zip(dynamic_threads, dynamic_avgs):
        print(f"{t:7d}, {avg:11.3f}, {len(dynamic_data[t]):7d}")

    # --- Process Static Data ---
    static_threads = sorted(static_data.keys())
    static_avgs = [sum(static_data[t]) / len(static_data[t]) for t in static_threads]

    print("\n--- Static Schedule Averages ---")
    print("Threads, AvgTime_ms, NumRuns")
    for t, avg in zip(static_threads, static_avgs):
        print(f"{t:7d}, {avg:11.3f}, {len(static_data[t]):7d}")

    # --- Calculate Speedup ---
    
    # Find the sequential baseline (1-thread static time)
    sequential_baseline_time = None
    if 1 in static_threads:
        sequential_baseline_time = static_avgs[static_threads.index(1)]
        print(f"\nUsing sequential baseline (1-thread static): {sequential_baseline_time:.3f} ms")
    else:
        print("\nError: Cannot calculate speedup.")
        print("No data found for 1-thread static runs.")
        print("Please run './kmeans_parallel 1 0' at least once.")
        return

    # Calculate speedups
    dynamic_speedups = [sequential_baseline_time / avg for avg in dynamic_avgs]
    static_speedups = [sequential_baseline_time / avg for avg in static_avgs]

    print("\n--- Dynamic Schedule Speedup ---")
    print("Threads, Speedup")
    for t, speedup in zip(dynamic_threads, dynamic_speedups):
        print(f"{t:7d}, {speedup:7.2f}x")

    print("\n--- Static Schedule Speedup ---")
    print("Threads, Speedup")
    for t, speedup in zip(static_threads, static_speedups):
        print(f"{t:7d}, {speedup:7.2f}x")

    # --- Plot 1: Average Runtimes ---
    plt.figure(figsize=(12, 7), dpi=150)
    plt.plot(dynamic_threads, dynamic_avgs, "-o", label="Dynamic Schedule")
    plt.plot(static_threads, static_avgs, "-o", label="Static Schedule")
    plt.title("K-Means Average Runtime vs Threads")
    plt.xlabel("Number of Threads")
    plt.ylabel("Average Time (ms)")
    plt.xticks(sorted(list(set(dynamic_threads + static_threads)))) # Show all thread counts
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.savefig(TIME_PLOT_PATH)
    print(f"\nSaved runtime plot to '{TIME_PLOT_PATH}'")

    # --- Plot 2: Speedup ---
    plt.figure(figsize=(12, 7), dpi=150)
    
    # Plot ideal speedup line
    all_threads = sorted(list(set(dynamic_threads + static_threads)))
    plt.plot(all_threads, all_threads, "--", label="Ideal Speedup", color="gray")
    
    # Plot actual speedups
    plt.plot(dynamic_threads, dynamic_speedups, "-o", label="Dynamic Schedule")
    plt.plot(static_threads, static_speedups, "-o", label="Static Schedule")
    
    plt.title(f"K-Means Speedup vs Threads (Baseline: {sequential_baseline_time:.2f} ms)")
    plt.xlabel("Number of Threads")
    plt.ylabel("Speedup Factor (X)")
    plt.xticks(all_threads)
    plt.yticks(all_threads) # Set y-ticks to match x-ticks for easy "ideal" comparison
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.savefig(SPEEDUP_PLOT_PATH)
    print(f"Saved speedup plot to '{SPEEDUP_PLOT_PATH}'")


if __name__ == "__main__":
    main()