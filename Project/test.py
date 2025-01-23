import numpy as np
import matplotlib.pyplot as plt
import os

def read_bin_file(filepath, shape):
    data = np.fromfile(filepath, dtype=np.float64)
    return data.reshape(shape)

def visualize_velocity_magnitude(step):
    N = 128
    u_x = read_bin_file(f'./bin_results/u_x_{step}.bin', (N, N))
    u_y = read_bin_file(f'./bin_results/u_y_{step}.bin', (N, N))

    # Calculate the magnitude of the velocity field
    magnitude = np.sqrt(u_x**2 + u_y**2)

    # Plot the heatmap
    plt.imshow(magnitude, cmap='hot', interpolation='nearest')
    plt.colorbar(label='Velocity Magnitude')
    plt.title(f'Velocity Magnitude Heatmap at Step {step}')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()

# Example usage
visualize_velocity_magnitude(0)