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
            x.append(int(row[1]))
            y.append(int(row[3]))
            z.append(int(row[4]))
    return x, y, z

def plot_data_time(x: List[int], y: List[int]):
    plt.scatter(x, y, label='Data points')
    fx = [i * np.log(i) for i in x]
    plt.scatter(x, fx, color='orange', label='y = x log(x)')
    plt.title('Timing Plot')
    plt.xlabel('Number of Edges')
    plt.ylabel('Time (μs)')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_data_memory(x: List[int], z: List[int]):
    plt.scatter(x, z, label='Data points')
    fx = [i * np.log(i) for i in x]
    plt.scatter(x, fx, color='orange', label='y = x log(x)')
    plt.title('Memory Plot')
    plt.xlabel('Number of Edges')
    plt.ylabel('Memory (KB)')
    plt.legend()
    plt.grid(True)
    plt.show()

print("Starting plot script.")
x, y, z = read_csv("q3_tests/q3_timings.csv")
print("Printing graph.")
plot_data_time(x, y)
plot_data_memory(x, z)