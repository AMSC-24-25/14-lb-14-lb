import numpy as np
import matplotlib.pyplot as plt
import os
from PIL import Image
import pandas as pd
import vtk
from vtk.util import numpy_support

# Directory dove salvare i file di output
output_directory = './gifs/'
vtk_output_directory = './vtk/'
os.makedirs(output_directory, exist_ok=True)
os.makedirs(vtk_output_directory, exist_ok=True)

def export_to_vtk(rho, ux, uy, nx, ny, timestep):
    """
    Esporta i dati della simulazione in formato VTK per la visualizzazione in ParaView.
    """
    # Crea un oggetto vtkImageData
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(nx, ny, 1)
    imageData.SetSpacing(1.0, 1.0, 1.0)

    # Converti i dati numpy in array VTK
    vtk_rho = numpy_support.numpy_to_vtk(rho.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    vtk_ux = numpy_support.numpy_to_vtk(ux.ravel(), deep=True, array_type=vtk.VTK_FLOAT)
    vtk_uy = numpy_support.numpy_to_vtk(uy.ravel(), deep=True, array_type=vtk.VTK_FLOAT)

    # Aggiungi i dati come scalari al vtkImageData
    vtk_rho.SetName('Density')
    vtk_ux.SetName('Velocity_U')
    vtk_uy.SetName('Velocity_V')
    imageData.GetPointData().AddArray(vtk_rho)
    imageData.GetPointData().AddArray(vtk_ux)
    imageData.GetPointData().AddArray(vtk_uy)

    # Salva i dati in un file VTK
    writer = vtk.vtkXMLImageDataWriter()
    vtk_filename = os.path.join(vtk_output_directory, f'simulation_timestep_{timestep:04d}.vti')
    writer.SetFileName(vtk_filename)
    writer.SetInputData(imageData)
    writer.Write()

def create_gifs_from_csv(csv_filename):
    """
    Crea delle GIF per i campi di velocità e densità a partire dai dati salvati in un file CSV.
    """
    # Leggi il file CSV
    df = pd.read_csv(csv_filename)
    
    # Ottieni il numero di timestep unici
    timesteps = df['Timestep'].unique()
    
    density_images = []
    velocity_images = []
    
    # Itera su ogni timestep
    for t in timesteps:
        # Filtra i dati per il timestep corrente
        df_t = df[df['Timestep'] == t]
        
        # Crea una matrice per la densità e per le componenti della velocità
        nx = df_t['X'].max() + 1
        ny = df_t['Y'].max() + 1
        rho = df_t.pivot(index='Y', columns='X', values='Density').values
        ux = df_t.pivot(index='Y', columns='X', values='Velocity_U').values
        uy = df_t.pivot(index='Y', columns='X', values='Velocity_V').values
        
        # Esporta i dati in formato VTK
        export_to_vtk(rho, ux, uy, nx, ny, t)
        
        # Calcola il modulo della velocità
        velocity_magnitude = np.sqrt(ux**2 + uy**2)
        
        # Crea la figura per la densità
        fig, ax = plt.subplots()
        cax = ax.imshow(rho, origin='lower', cmap='viridis')
        fig.colorbar(cax)
        ax.set_title(f'Densità - Timestep {t}')
        plt.axis('off')
        
        # Salva l'immagine come temporanea
        temp_density_path = os.path.join(output_directory, f'temp_rho_{t:04d}.png')
        plt.savefig(temp_density_path, bbox_inches='tight')
        plt.close(fig)
        density_images.append(Image.open(temp_density_path))
        
        # Crea la figura per la velocità
        fig, ax = plt.subplots()
        cax = ax.imshow(velocity_magnitude, origin='lower', cmap='plasma')
        fig.colorbar(cax)
        ax.set_title(f'Modulo Velocità - Timestep {t}')
        plt.axis('off')
        
        # Salva l'immagine come temporanea
        temp_velocity_path = os.path.join(output_directory, f'temp_vel_{t:04d}.png')
        plt.savefig(temp_velocity_path, bbox_inches='tight')
        plt.close(fig)
        velocity_images.append(Image.open(temp_velocity_path))
    
    # Crea la GIF per la densità
    density_gif_path = os.path.join(output_directory, 'density_evolution_from_csv.gif')
    density_images[0].save(density_gif_path, save_all=True, append_images=density_images[1:], duration=200, loop=0)
    
    # Crea la GIF per la velocità
    velocity_gif_path = os.path.join(output_directory, 'velocity_evolution_from_csv.gif')
    velocity_images[0].save(velocity_gif_path, save_all=True, append_images=velocity_images[1:], duration=200, loop=0)
    
    # Rimuovi i file temporanei
    for img in density_images + velocity_images:
        os.remove(img.filename)
    
    print(f'GIF di densità salvata in: {density_gif_path}')
    print(f'GIF di velocità salvata in: {velocity_gif_path}')

# Nome del file CSV
csv_filename = "simulation_data.csv"

# Crea le GIF dai dati del CSV
create_gifs_from_csv(csv_filename)
