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
            z.append(int(row[3]))
    return x, y, z

def plot_data_3d(x: List[int], y: List[int], z: List[int]):
    #add 3d plot code here
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    fxy = []
    for i in range(len(x)):
        fxy.append(8 * max(x[i], y[i]) * np.log(x[i]))
    ax.scatter(x, y, z, color = 'blue', label='Data points')
    ax.scatter(x, y, fxy, color = 'orange', label='Reference points')
    ax.set_title('3D Timing Plot')
    ax.set_xlabel('Number of Nodes')
    ax.set_ylabel('Number of Edges')
    ax.set_zlabel('Time (μs)')
    ax.legend()
    
    plt.show()

print("Starting plot script.")
x, y, z = read_csv("q3_tests/q3_timings.csv")
print("Printing graph.")
plot_data_3d(x, y, z)