import csv
from typing import List, Tuple
import matplotlib.pyplot as plt
import numpy as np

def read_csv(file_path: str) -> Tuple[List[int], List[int]]:
    print("Reading file.")
    x, y, z = [], [], []
    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header
        for row in reader:
            x.append(int(row[0]))
            y.append(int(row[1]))
            z.append(int(row[2]))
    return x, y, z

def plot_data_time(x: List[int], y: List[int]):
    plt.scatter(x, y, label='Data points')
    plt.scatter(x, x, label='Reference points y = x', color='orange')
    plt.title('Timing Plot')
    plt.xlabel('Number of parcels')
    plt.ylabel('Time (μs)')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_data_memory(x: List[int], z: List[int]):
    plt.scatter(x, z, label='Data points')
    plt.scatter(x, x, label='Reference points y = x', color='orange')
    plt.title('Memory Plot')
    plt.xlabel('Number of parcels')
    plt.ylabel('Memory (KB)')
    plt.legend()
    plt.grid(True)
    plt.show()

print("Starting plot script.")
x, y, z = read_csv("q1_tests/q1_timings.csv")
print("Printing graph.")
plot_data_time(x, y)
plot_data_memory(x, z)