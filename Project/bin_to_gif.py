import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import imageio
import os

# Load simulation parameters
params = pd.read_csv('simulation_parameters.csv')
NX, NY, NZ, NSTEPS, SAVE_EVERY, Re = params.iloc[0][['NX', 'NY', 'NZ', 'NSTEPS', 'NSAVE', 'RE']].astype(int).values
UMAX = params.iloc[0]['U_MAX']

def read_binary_file(filename, shape):
    """Legge un file binario con double e lo converte in un array numpy di forma 'shape'."""
    data = np.fromfile(filename, dtype=np.float64)
    return data.reshape(shape)

def create_velocity_gif(input_folder, output_gif, nx, ny, nz, steps, slice_z=0, fps=5):
    """
    Crea una GIF che mostra la mappa 2D del modulo di velocità.
    
    Parametri:
    - input_folder: cartella in cui si trovano i file binari u_x_*.bin, u_y_*.bin, u_z_*.bin
    - output_gif: nome (o path) del file GIF in uscita
    - nx, ny, nz: dimensioni del dominio
    - steps: lista o range dei 'time step' da visualizzare
    - slice_z: se >0 e il dominio è 3D, indica l'indice z in cui estrarre la sezione
    - fps: numero di frame al secondo nella GIF
    """
    frames = []  # raccogliamo qui le immagini (in memoria) per comporre la GIF
    
    # Creiamo una figura fuori dal loop così da riusare la stessa
    fig, ax = plt.subplots()
    
    for step in steps:
        # Costruiamo i nomi dei file binari per le 3 componenti
        ux_file = os.path.join(input_folder, f"u_x_{step}.bin")
        uy_file = os.path.join(input_folder, f"u_y_{step}.bin")
        uz_file = os.path.join(input_folder, f"u_z_{step}.bin")
        
        # Leggiamo i file (come array 3D, forma (nz, ny, nx))
        u_x = read_binary_file(ux_file, (nz, ny, nx))
        u_y = read_binary_file(uy_file, (nz, ny, nx))
        u_z = read_binary_file(uz_file, (nz, ny, nx))
        
        # Calcoliamo il modulo della velocità
        # shape: (nz, ny, nx)
        velocity_magnitude = np.sqrt(u_x**2 + u_y**2 + u_z**2)
        
        # Se il dominio è 2D (nz=1) oppure vogliamo una slice specifica in 3D
        # estraiamo la fetta z = slice_z
        # (attenzione che 'slice_z' sia un indice valido se nz>1)
        velocity_2d = velocity_magnitude[:, slice_z, :]
        
        # Puliamo l'axes per disegnare un nuovo frame
        ax.clear()
        
        # Mostriamo la mappa 2D
        # origin='lower' per avere y dal basso in alto, dipende dal tuo convenzione
        im = ax.imshow(velocity_2d, cmap='jet', origin='lower', 
                       extent=(0, nx, 0, ny), vmin=0, vmax=velocity_2d.max())
        
        ax.set_title(f"Step = {step}, z = {slice_z}")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        
        # Aggiorniamo la colorbar
        # Per evitare di creare 1000 colorbar, possiamo farlo una volta fuori dal ciclo
        # o possiamo rimuovere e ricreare. Qui per semplicità la creiamo la prima volta e poi
        # la riusiamo.
        
        # Convertiamo la figura in un'immagine RGBA
        fig.canvas.draw()
        # Leggiamo i pixel dalla figura
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        
        frames.append(image)
    
    plt.close(fig)  # chiudiamo la figura quando abbiamo finito

    # Creiamo la GIF
    # duration=1/fps => durata di ogni frame
    imageio.mimsave(output_gif, frames, fps=fps)
    print(f"GIF salvata in: {output_gif}")

# ESEMPIO di utilizzo
if __name__ == "__main__":
    input_folder = "./bin_results"
    output_gif = "vel_field.gif"
    steps = range(0, NSTEPS, SAVE_EVERY)   # ad esempio generiamo un frame ogni 40 step
    slice_z = 64               # se vogliamo la "metà" del dominio in z
    fps = 5
    
    create_velocity_gif(input_folder, output_gif, NX, NY, NZ, steps, slice_z, fps)
