import numpy as np
import matplotlib.pyplot as plt
import os
from PIL import Image
import csv
import pandas as pd

# Load simulation parameters
params = pd.read_csv('simulation_parameters.csv')
NX, NY, NSTEPS, SAVE_EVERY = params.iloc[0][['NX', 'NY', 'NSTEPS', 'NSAVE']].astype(int).values
UMAX = params.iloc[0]['UMAX']

# Path for read and write
data_directory = './bin_results/'
output_directory = './gifs/'
os.makedirs(output_directory, exist_ok=True)

def read_binary_file(filename, shape):
    """
    Legge un file binario e restituisce un array numpy con la forma specificata.
    """
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float64)
    return data.reshape(shape)

def create_gifs(nx, ny, nsteps, save_every):
    """
    Crea delle GIF per i campi di velocità e densità a partire dai file binari.
    """
    density_images = []
    velocity_images = []
    
    # Iterate on every timestep
    for n in range(0, nsteps + 1, save_every):
        # Input
        #rho_file = os.path.join(data_directory, f'rho{n}.bin')
        ux_file = os.path.join(data_directory, f'u_x_{n}.bin')
        uy_file = os.path.join(data_directory, f'u_y_{n}.bin')
        
        # Read data
        #rho = read_binary_file(rho_file, (ny, nx))
        ux = read_binary_file(ux_file, (ny, nx))
        uy = read_binary_file(uy_file, (ny, nx))
        
        # Compute velocity module
        velocity_magnitude = np.sqrt(ux**2 + uy**2) / UMAX
        
        # Create velocity image
        fig, ax = plt.subplots()
        cax = ax.imshow(velocity_magnitude, origin='lower', cmap='plasma', vmin=0, vmax=1)
        fig.colorbar(cax)
        ax.set_title(f'Modulo Velocità - Timestep {n}')
        plt.axis('off')
        
        # Save temp images
        temp_velocity_path = os.path.join(output_directory, f'temp_vel_{n:04d}.png')
        plt.savefig(temp_velocity_path, bbox_inches='tight')
        plt.close(fig)
        velocity_images.append(Image.open(temp_velocity_path))
    
    
    # Create velocity GIF
    velocity_gif_path = os.path.join(output_directory, 'velocity_evolution.gif')
    velocity_images[0].save(velocity_gif_path, save_all=True, append_images=velocity_images[1:], duration=200, loop=0)
    
    # Remove temp files
    for img in velocity_images:
        os.remove(img.filename)

    print(f'Velocity GIF saved in: {velocity_gif_path}')

# Create GIFs
create_gifs(NX, NY, NSTEPS, SAVE_EVERY)
