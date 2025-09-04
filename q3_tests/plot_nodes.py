import csv
from typing import List, Tuple
import matplotlib.pyplot as plt
import numpy as np

def read_csv(file_path: str) -> Tuple[List[int], List[int]]:
    print("Reading file.")
    x, y = [], []
    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header
        for row in reader:
            x.append(int(row[0]))
            y.append(int(row[3]))
    return x, y

def plot_data(x: List[int], y: List[int]):
    plt.scatter(x, y, label='Data points')
    fx = [50 * i * np.log(i) for i in x]
    plt.scatter(x, fx, color='orange', label='y = x log(x)')
    plt.title('Timing Plot')
    plt.xlabel('Number of Nodes')
    plt.ylabel('Time (μs)')
    plt.legend()
    plt.show()

print("Starting plot script.")
x, y = read_csv("q3_tests/q3_timings.csv")
print("Printing graph.")
plot_data(x, y)