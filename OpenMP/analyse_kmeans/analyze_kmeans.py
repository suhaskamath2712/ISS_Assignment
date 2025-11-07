"""
Simple script to compute average Time_ms per Threads from hardcoded data
and plot a line graph.

Usage: python analyze_kmeans.py
Outputs: prints averages and saves kmeans_avg_time.png
"""

import matplotlib.pyplot as plt

DATA = {
    # Timings (ms) taken from OpenMP/cluster/output.log (five runs each)
    1:  [27.744, 27.768, 27.860, 27.598, 27.593],
    4:  [8.852, 9.075, 8.878, 8.961, 8.843],
    8:  [5.512, 5.526, 5.650, 5.783, 6.081],
    16: [4.548, 4.843, 4.680, 4.597, 4.665],
    32: [3.776, 3.705, 3.668, 3.650, 3.655],
}


def main():
    threads = sorted(DATA.keys())
    avgs = [sum(DATA[t]) / len(DATA[t]) for t in threads]

    # Print averages
    print("Threads, AvgTime_ms, NumRuns")
    for t, avg in zip(threads, avgs):
        print(f"{t}, {avg:.3f}, {len(DATA[t])}")

    # Plot
    plt.figure(figsize=(7, 4.5), dpi=140)
    plt.plot(threads, avgs, "-o")
    plt.title("K-Means Average Runtime vs Threads")
    plt.xlabel("Threads")
    plt.ylabel("Average Time (ms)")
    plt.xticks(threads)
    plt.grid(True, linestyle=":", linewidth=0.7, alpha=0.6)
    plt.tight_layout()
    plt.savefig("kmeans_avg_time.png")


if __name__ == "__main__":
    main()
