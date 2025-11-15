import re
import os
import numpy as np
import matplotlib.pyplot as plt

# --- Configuration ---
CSV_FILE = r"H:\My Drive\Introduction to Scalable Systems\ISS_Assignment\MPI\times.csv"
EXPECTED_RUNS = 5

def analyze_data(csv_file):
    """
    Loads the CSV data and calculates:
    1. Average total time for each process count.
    2. Average computation time for each process count.
    3. Load imbalance factor for each run.
    """
    
    comp_times = {} # Stores computation times: {procs: [run1_times, run2_times, ...]}
    total_times = {} # Stores total times: {procs: [run1_total, run2_total, ...]}

    # Base directory for output plots (use same directory as the CSV file)
    base_dir = os.path.dirname(os.path.abspath(csv_file)) if csv_file else os.getcwd()

    try:
        with open(csv_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Could not find the file '{csv_file}'.")
        print("Please make sure the CSV file is in the same directory.")
        return
    except Exception as e:
        print(f"Error loading CSV: {e}")
        return

    current_comp_run = []
    current_comp_procs = 0 # This will hold the procs for the run we are building

    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        # Line is a computation time line
        if line.startswith('#'):
            try:
                _, procs, rank, time_ms = line.split(',')
                procs = int(procs)
                time_ms = float(time_ms)
                
                # Set the process count for the run we are currently building
                if current_comp_procs == 0:
                    current_comp_procs = procs
                
                current_comp_run.append(time_ms)
                
            except ValueError:
                print(f"Warning: Skipping malformed line: {line}")
                
        # Line is a total time line
        elif line.startswith('$'):
            try:
                _, procs, time_ms = line.split(',')
                procs = int(procs)
                time_ms = float(time_ms)
                
                if procs not in total_times:
                    total_times[procs] = []
                total_times[procs].append(time_ms)
                
                # --- THIS IS THE FIX ---
                # The '$' line marks the END of a run.
                # We save the computation run that just finished
                # and then reset for the next run.
                
                if current_comp_run:
                    if current_comp_procs not in comp_times:
                        comp_times[current_comp_procs] = []
                    comp_times[current_comp_procs].append(current_comp_run)
                
                # Reset for the next run
                current_comp_run = []
                current_comp_procs = 0
                # --- END FIX ---
                
            except ValueError:
                print(f"Warning: Skipping malformed line: {line}")

    # --- Start Analysis ---
    print("--- Performance Analysis ---")
    
    process_counts = sorted(list(set(comp_times.keys()) | set(total_times.keys())))
    
    avg_comp_times_dict = {} # For speedup calculation and plotting (avg comp time per procs)
    avg_total_times_dict = {} # For plotting (avg total time per procs)
    avg_imbalance_by_procs = {} # For plotting (avg load factor per procs)
    
    for procs in process_counts:
        print(f"\n=== Analysis for {procs} Processes ===")
        
        # --- 1. Total Time ---
        if procs in total_times:
            runs = total_times[procs]
            if len(runs) != EXPECTED_RUNS:
                print(f"Warning: Found {len(runs)} total time runs, expected {EXPECTED_RUNS}.")
            avg_total_time = np.mean(runs)
            print(f"  Average TOTAL Time:   {avg_total_time:>8.2f} ms")
            avg_total_times_dict[procs] = float(avg_total_time)
        else:
            print("  Average TOTAL Time:   N/A")

        # --- 2. Computation Time & Imbalance ---
        if procs in comp_times:
            runs_matrix = comp_times[procs]
            if len(runs_matrix) != EXPECTED_RUNS:
                print(f"Warning: Found {len(runs_matrix)} computation runs, expected {EXPECTED_RUNS}.")

            # Time for each run is the max (slowest) process
            time_per_run = [np.max(run) for run in runs_matrix]
            avg_comp_time = np.mean(time_per_run)
            avg_comp_times_dict[procs] = avg_comp_time
            print(f"  Average COMP Time:  {avg_comp_time:>8.2f} ms")
            
            if procs > 1:
                print(f"  Load Imbalance Factor (per run):")
                imbalances_for_procs = []
                for i, run in enumerate(runs_matrix):
                    run = np.array(run) # Ensure it's a numpy array for min/max
                    max_time = np.max(run)
                    min_time = np.min(run)
                    imbalance = (max_time - min_time) / min_time if min_time > 0 else 0
                    imbalances_for_procs.append(imbalance)
                    print(f"    Run {i+1}: {imbalance:7.4f} (Max: {max_time:.2f} ms, Min: {min_time:.2f} ms)")
                if imbalances_for_procs:
                    avg_imbalance_by_procs[procs] = float(np.mean(imbalances_for_procs))
        else:
            print("  Average COMP Time:  N/A")

    # --- 3. Speedup (based on Computation Time) ---
    print("\n=== Computation Speedup vs. Serial ===")
    if 1 in avg_comp_times_dict:
        serial_avg = avg_comp_times_dict[1]
        print(f"Serial (1 proc): {serial_avg:>8.2f} ms (1.00x)")
        for procs, avg_time in avg_comp_times_dict.items():
            if procs > 1:
                speedup = serial_avg / avg_time
                print(f"Parallel ({procs:<3} procs): {avg_time:>8.2f} ms ({speedup:5.2f}x)")
    else:
        print("Serial (1 process) computation data not found. Cannot calculate speedup.")

    # --- 4. Average Load Factor by Process Count ---
    print("\n=== Average Load Factor by Process Count ===")
    if comp_times:
        # Recompute averages where missing and print summary
        for procs in sorted(comp_times.keys()):
            if procs not in avg_imbalance_by_procs:
                runs_matrix = comp_times[procs]
                imbalances = []
                for run in runs_matrix:
                    run_arr = np.array(run)
                    max_time = np.max(run_arr)
                    min_time = np.min(run_arr)
                    imbalance = (max_time - min_time) / min_time if min_time > 0 else 0.0
                    imbalances.append(imbalance)
                if imbalances:
                    avg_imbalance_by_procs[procs] = float(np.mean(imbalances))
            avg_imb = avg_imbalance_by_procs.get(procs, None)
            if avg_imb is not None:
                print(f"{procs:>4} procs: average load factor = {avg_imb:.4f}")
            else:
                print(f"{procs:>4} procs: average load factor = N/A")
    else:
        print("No computation time data to compute load factors.")

    # --- 5. Plots ---
    try:
        # 5.1 Average computation time vs process count
        widescreen = (16, 9)
        high_dpi = 600
        if avg_comp_times_dict:
            x = sorted(avg_comp_times_dict.keys())
            y = [avg_comp_times_dict[p] for p in x]
            plt.figure(figsize=widescreen)
            plt.plot(x, y, marker='o', linewidth=2)
            plt.title('Average Computation Time vs Process Count')
            plt.xlabel('Process Count')
            plt.ylabel('Avg Computation Time (ms)')
            plt.grid(True, linestyle='--', alpha=0.4)
            plt.tight_layout()
            plt.savefig(os.path.join(base_dir, 'avg_comp_time.png'), dpi=high_dpi)
            plt.close()

        # 5.2 Average total (execution) time vs process count
        if avg_total_times_dict:
            x = sorted(avg_total_times_dict.keys())
            y = [avg_total_times_dict[p] for p in x]
            plt.figure(figsize=widescreen)
            plt.plot(x, y, marker='o', color='tab:orange', linewidth=2)
            plt.title('Average Total Execution Time vs Process Count')
            plt.xlabel('Process Count')
            plt.ylabel('Avg Total Time (ms)')
            plt.grid(True, linestyle='--', alpha=0.4)
            plt.tight_layout()
            plt.savefig(os.path.join(base_dir, 'avg_total_time.png'), dpi=high_dpi)
            plt.close()

        # 5.3 Average load imbalance vs process count
        if avg_imbalance_by_procs:
            x = sorted(avg_imbalance_by_procs.keys())
            y = [avg_imbalance_by_procs[p] for p in x]
            plt.figure(figsize=widescreen)
            plt.plot(x, y, marker='o', color='tab:green', linewidth=2)
            plt.title('Average Load Imbalance vs Process Count')
            plt.xlabel('Process Count')
            plt.ylabel('Avg Load Factor ((max-min)/min)')
            plt.grid(True, linestyle='--', alpha=0.4)
            plt.tight_layout()
            plt.savefig(os.path.join(base_dir, 'avg_load_imbalance.png'), dpi=high_dpi)
            plt.close()

        # 5.4 Average computation time per process rank for 64 processes (across all runs)
        target_procs = 64
        if target_procs in comp_times:
            runs_matrix = comp_times[target_procs]
            runs_matrix_np = np.array(runs_matrix)  # shape: (num_runs, num_ranks)
            avg_per_rank = np.mean(runs_matrix_np, axis=0)  # shape: (num_ranks,)
            x = np.arange(len(avg_per_rank))
            plt.figure(figsize=widescreen)
            plt.plot(x, avg_per_rank, color='tab:blue', alpha=0.85)
            plt.title('Average Computation Time per Process Rank (64 Processes)')
            plt.xlabel('Rank')
            plt.ylabel('Avg Computation Time (ms)')
            plt.grid(True, linestyle='--', axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(base_dir, 'comp_times_64procs.png'), dpi=high_dpi)
            plt.close()
        else:
            print("Note: No 64-process computation runs found; skipping per-rank plot.")
    except Exception as e:
        print(f"Plotting error: {e}")


if __name__ == "__main__":
    analyze_data(CSV_FILE)