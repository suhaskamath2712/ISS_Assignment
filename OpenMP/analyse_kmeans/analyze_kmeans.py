"""
Simple script to compute average Time_ms per Threads from hardcoded data
and plot a line graph.

Usage: python analyze_kmeans.py
Outputs: prints averages and saves kmeans_avg_time.png
"""

import matplotlib.pyplot as plt

#Read data from CSV output of output_parallel.csv
#Returns two 2-D arrays: dynamic_data and static_data
def read_data():
    static_data = {}
    dynamic_data = {}
    file_path = r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\OpenMP\analyse_kmeans\output_parallel.csv"
    with open(file_path, "r") as f:
        for line in f:
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
    return dynamic_data, static_data


def main():
    dynamic_data, static_data = read_data()
    dynamic_threads = sorted(dynamic_data.keys())
    dynamic_avgs = [sum(dynamic_data[t]) / len(dynamic_data[t]) for t in dynamic_threads]

    # Print averages
    print("Dynamic Schedule Averages:")
    print("Threads, AvgTime_ms, NumRuns")
    for t, avg in zip(dynamic_threads, dynamic_avgs):
        print(f"{t}, {avg:.3f}, {len(dynamic_data[t])}")

    static_threads = sorted(static_data.keys())
    static_avgs = [sum(static_data[t]) / len(static_data[t]) for t in static_threads]

    # Print averages
    print("Static Schedule Averages:")
    print("Threads, AvgTime_ms, NumRuns")
    for t, avg in zip(static_threads, static_avgs):
        print(f"{t}, {avg:.3f}, {len(static_data[t])}")

    # Plot
    plt.figure(figsize=(16,9), dpi=300)
    plt.plot(dynamic_threads, dynamic_avgs, "-o", label="Dynamic")
    plt.plot(static_threads, static_avgs, "-o", label="Static")
    plt.title("K-Means Average Runtime vs Threads")
    plt.xlabel("Threads")
    plt.ylabel("Average Time (ms)")
    plt.xticks(dynamic_threads)
    plt.grid(True)
    plt.legend()
    plt.savefig("kmeans_avg_time.png")

if __name__ == "__main__":
    main()
