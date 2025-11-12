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
        # We assume the user removed '#' symbols as stated.
        # If the header is 'processes,rank,time_ms', this works automatically.
        df = pd.read_csv(FILE_NAME)
        
        # robustness: strip whitespace from column names just in case
        df.columns = df.columns.str.strip()
        
        # robustness: Rename columns if they differ slightly (e.g. if '#' was left in header)
        if '#processes' in df.columns:
            df.rename(columns={'#processes': 'processes'}, inplace=True)
            
    except FileNotFoundError:
        print(f"Error: {FILE_NAME} not found.")
        return

    # Ensure data is numeric
    df['time_ms'] = pd.to_numeric(df['time_ms'], errors='coerce')
    df.dropna(subset=['time_ms'], inplace=True)

    # Get unique process counts sorted
    process_counts = sorted(df['processes'].unique())
    
    summary_data = []

    print(f"{'Processes':<10} | {'Avg Time (ms)':<15} | {'Avg Imbalance':<15}")
    print("-" * 45)

    # 2. Calculate Parameters
    for p in process_counts:
        # Filter data for this process count
        subset = df[df['processes'] == p]
        
        # Extract the time column
        times = subset['time_ms'].values
        
        # Reshape: (Number of Runs, Number of Processes per Run)
        # If we have 5 runs of 8 processes, we reshape to (5, 8)
        try:
            # Check if data points match expected count
            if len(times) % p != 0:
                print(f"Warning: Data count mismatch for {p} processes. Skipping.")
                continue
                
            runs_matrix = times.reshape(-1, p)
            
        except ValueError:
            print(f"Error reshaping data for {p} processes.")
            continue

        # --- CALCULATION 1: Time taken by each run ---
        # The execution time of a parallel run is determined by the SLOWEST process (max time)
        run_times = np.max(runs_matrix, axis=1)
        
        # --- CALCULATION 2: Average time ---
        avg_time = np.mean(run_times)
        
        # --- CALCULATION 3: Load Imbalance Factor ---
        # Formula: (Max - Min) / Min
        max_times = np.max(runs_matrix, axis=1)
        min_times = np.min(runs_matrix, axis=1)
        
        # Handle division by zero for very fast/dummy runs
        with np.errstate(divide='ignore', invalid='ignore'):
            imbalance = (max_times - min_times) / min_times
            # For 1 process, min=max, result is 0/0 = NaN or 0. We force 0.
            imbalance = np.nan_to_num(imbalance)
            
        avg_imbalance = np.mean(imbalance)

        summary_data.append({
            'processes': p,
            'avg_time': avg_time,
            'avg_imbalance': avg_imbalance
        })
        
        print(f"{p:<10} | {avg_time:<15.2f} | {avg_imbalance:<15.4f}")

    # Convert summary to DataFrame for easy plotting
    plot_df = pd.DataFrame(summary_data)

    # 3. Plotting
    sns.set_style("whitegrid")
    
    # --- Plot 1: Avg Time vs Number of Processes ---
    plt.figure(figsize=(10, 6))
    plt.plot(plot_df['processes'], plot_df['avg_time'], marker='o', linewidth=2, color='b')
    plt.title('Average Execution Time vs. Number of Processes', fontsize=14)
    plt.xlabel('Number of Processes', fontsize=12)
    plt.ylabel('Time (ms)', fontsize=12)
    plt.xticks(process_counts) # Ensure all process counts are shown on X-axis
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Save Plot 1
    plt.savefig(r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\MPI\plot_execution_time.png")
    print("\nSaved plot: plot_execution_time.png")
    plt.close()

    # --- Plot 2: Avg Load Imbalance vs Number of Processes ---
    plt.figure(figsize=(10, 6))
    plt.plot(plot_df['processes'], plot_df['avg_imbalance'], marker='s', linewidth=2, color='r')
    plt.title('Average Load Imbalance Factor vs. Number of Processes', fontsize=14)
    plt.xlabel('Number of Processes', fontsize=12)
    plt.ylabel('Load Imbalance Factor ((Max-Min)/Min)', fontsize=12)
    plt.xticks(process_counts)
    plt.grid(True, linestyle='--', alpha=0.7)

    # Save Plot 2
    plt.savefig(r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\MPI\plot_load_imbalance.png")
    print("Saved plot: plot_load_imbalance.png")
    plt.close()

if __name__ == "__main__":
    analyze_and_plot()