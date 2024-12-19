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

# Directory dove sono salvati i file binari
data_directory = '../bin_results/'  # Modifica questo percorso con quello corretto
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
    
    # Itera su ogni timestep salvato
    for n in range(0, nsteps + 1, save_every):
        # File di input
        rho_file = os.path.join(data_directory, f'rho{n:04d}.bin')
        ux_file = os.path.join(data_directory, f'ux{n:04d}.bin')
        uy_file = os.path.join(data_directory, f'uy{n:04d}.bin')
        
        # Leggi i dati
        rho = read_binary_file(rho_file, (ny, nx))
        ux = read_binary_file(ux_file, (ny, nx))
        uy = read_binary_file(uy_file, (ny, nx))
        
        # Calcola il modulo della velocità
        velocity_magnitude = np.sqrt(ux**2 + uy**2)
        
        # Crea la figura per la densità
        fig, ax = plt.subplots()
        cax = ax.imshow(rho, origin='lower', cmap='viridis')
        fig.colorbar(cax)
        ax.set_title(f'Densità - Timestep {n}')
        plt.axis('off')
        
        # Salva l'immagine come temporanea
        temp_density_path = os.path.join(output_directory, f'temp_rho_{n:04d}.png')
        plt.savefig(temp_density_path, bbox_inches='tight')
        plt.close(fig)
        density_images.append(Image.open(temp_density_path))
        
        # Crea la figura per la velocità
        fig, ax = plt.subplots()
        cax = ax.imshow(velocity_magnitude, origin='lower', cmap='plasma')
        fig.colorbar(cax)
        ax.set_title(f'Modulo Velocità - Timestep {n}')
        plt.axis('off')
        
        # Salva l'immagine come temporanea
        temp_velocity_path = os.path.join(output_directory, f'temp_vel_{n:04d}.png')
        plt.savefig(temp_velocity_path, bbox_inches='tight')
        plt.close(fig)
        velocity_images.append(Image.open(temp_velocity_path))
    
    # Crea la GIF per la densità
    density_gif_path = os.path.join(output_directory, 'density_evolution.gif')
    density_images[0].save(density_gif_path, save_all=True, append_images=density_images[1:], duration=200, loop=0)
    
    # Crea la GIF per la velocità
    velocity_gif_path = os.path.join(output_directory, 'velocity_evolution.gif')
    velocity_images[0].save(velocity_gif_path, save_all=True, append_images=velocity_images[1:], duration=200, loop=0)
    
    # Rimuovi i file temporanei
    for img in density_images + velocity_images:
        os.remove(img.filename)
    
    print(f'GIF di densità salvata in: {density_gif_path}')
    print(f'GIF di velocità salvata in: {velocity_gif_path}')

# Crea le GIF
create_gifs(NX, NY, NSTEPS, SAVE_EVERY)
