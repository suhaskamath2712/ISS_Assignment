import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Configuration
FILE_NAME = r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\MPI\computation_times_oversubscribed.csv"
EXPECTED_RUNS = 5

def analyze_and_plot():
    # 1. Load Data
    try:
        df = pd.read_csv(FILE_NAME)
        df.columns = df.columns.str.strip()
        if '#processes' in df.columns:
            df.rename(columns={'#processes': 'processes'}, inplace=True)
            
    except FileNotFoundError:
        print(f"Error: {FILE_NAME} not found.")
        return

    # Ensure data is numeric
    df['time_ms'] = pd.to_numeric(df['time_ms'], errors='coerce')
    df.dropna(subset=['time_ms'], inplace=True)

    process_counts = sorted(df['processes'].unique())
    
    summary_data = []
    average_times_dict = {} # To store avg times for speedup calculation

    print("--- Performance Analysis ---")
    print("\n=== Average Computation Time & Load Imbalance (across 5 runs) ===")

    # 2. Calculate Parameters
    for p in process_counts:
        subset = df[df['processes'] == p]
        times = subset['time_ms'].values
        
        try:
            if len(times) % p != 0:
                print(f"\nWarning: Data count mismatch for {p} processes. Skipping.")
                continue
            runs_matrix = times.reshape(-1, p)
            
        except ValueError:
            print(f"\nError reshaping data for {p} processes. Skipping.")
            continue

        # --- CALCULATION 1: Time taken by each run ---
        run_times = np.max(runs_matrix, axis=1)
        
        # --- CALCULATION 2: Average time ---
        avg_time = np.mean(run_times)
        average_times_dict[p] = avg_time # Store for speedup
        
        # --- CALCULATION 3: Load Imbalance Factor ---
        max_times = np.max(runs_matrix, axis=1)
        min_times = np.min(runs_matrix, axis=1)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            imbalance_per_run = (max_times - min_times) / min_times
            imbalance_per_run = np.nan_to_num(imbalance_per_run)
            
        avg_imbalance = np.mean(imbalance_per_run)

        # --- Add data for plotting ---
        summary_data.append({
            'processes': p,
            'avg_time': avg_time,
            'avg_imbalance': avg_imbalance
        })
        
        # --- RESTORED CONSOLE OUTPUT ---
        print(f"\nProcesses: {p:<3} | Avg. Time: {avg_time:>8.2f} ms")
        
        if p > 1: # Imbalance is meaningless for 1 process
            print(f"  Load Imbalance Factor (per run):")
            for i, imbalance_val in enumerate(imbalance_per_run):
                print(f"    Run {i+1}: {imbalance_val:7.4f} (Max: {max_times[i]:.2f} ms, Min: {min_times[i]:.2f} ms)")

    # --- RESTORED SPEEDUP SUMMARY ---
    print("\n=== Speedup vs. Serial (Based on Avg. Time) ===")
    if 1 in average_times_dict:
        serial_avg = average_times_dict[1]
        print(f"Serial (1 proc): {serial_avg:>8.2f} ms (1.00x)")
        for procs, avg_time in average_times_dict.items():
            if procs > 1:
                speedup = serial_avg / avg_time
                print(f"Parallel ({procs:<3} procs): {avg_time:>8.2f} ms ({speedup:5.2f}x)")
    else:
        print("Serial (1 process) data not found. Cannot calculate speedup.")


    # --- 3. Plotting ---
    plot_df = pd.DataFrame(summary_data)
    sns.set_style("whitegrid")
    
    # --- Plot 1: Avg Time vs Number of Processes ---
    plt.figure(figsize=(16,9))
    plt.plot(plot_df['processes'], plot_df['avg_time'], linewidth=2, marker='o')
    plt.title('Average Execution Time vs. Number of Processes')
    plt.xlabel('Number of Processes')
    plt.ylabel('Time (ms)')
    plt.xticks(process_counts)
    plt.grid(True)
    plt.savefig(r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\MPI\plot_execution_time.png", dpi=300)
    print("\nSaved plot: plot_execution_time.png")
    plt.close()

    # --- Plot 2: Avg Load Imbalance vs Number of Processes ---
    plt.figure(figsize=(16,9))
    plt.plot(plot_df['processes'], plot_df['avg_imbalance'], linewidth=2, marker='o')
    plt.title('Average Load Imbalance Factor vs. Number of Processes')
    plt.xlabel('Number of Processes')
    plt.ylabel('Load Imbalance Factor ((Max-Min)/Min)')
    plt.xticks(process_counts)
    plt.grid(True)
    plt.savefig(r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\MPI\plot_load_imbalance.png", dpi=300)
    print("Saved plot: plot_load_imbalance.png")
    plt.close()

    # --- Plot 3: Computation Time per Process for 64-Process Runs ---
    RANK_PLOT_PROCESS_COUNT = 64
    df_target = df[df['processes'] == RANK_PLOT_PROCESS_COUNT].copy()
    
    if df_target.empty:
        print(f"Info: No data found for {RANK_PLOT_PROCESS_COUNT} processes. Skipping rank plot.")
    else:
        # Group by 'rank' and calculate the mean 'time_ms' across the 5 runs
        avg_rank_time = df_target.groupby('rank')['time_ms'].mean().reset_index()
        overall_avg_time = avg_rank_time['time_ms'].mean()

        plt.figure(figsize=(16,9))
        plt.plot(avg_rank_time['rank'], avg_rank_time['time_ms'], linewidth=2, marker='o')
        plt.axhline(y=overall_avg_time, color='r', linestyle='--', label=f'Overall Avg. ({overall_avg_time:.2f} ms)')
        plt.title(f'Average Computation Time per Process (All {RANK_PLOT_PROCESS_COUNT}-Process Runs)', fontsize=16)
        plt.xlabel('Process Rank')
        plt.ylabel('Average Computation Time (ms)')
        
        plt.xticks(range(0, RANK_PLOT_PROCESS_COUNT, 8)) # Show a tick every 8 ranks
        plt.xlim(-1, RANK_PLOT_PROCESS_COUNT)
        
        plt.legend()
        plt.grid(True, which='both', linestyle='--', alpha=0.7)
        
        plot_filename = rf"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\MPI\plot_{RANK_PLOT_PROCESS_COUNT}_process_rank_times.png"
        plt.savefig(plot_filename, dpi=300)
        print(f"Saved plot: plot_{RANK_PLOT_PROCESS_COUNT}_process_rank_times.png")
        plt.close()

if __name__ == "__main__":
    analyze_and_plot()