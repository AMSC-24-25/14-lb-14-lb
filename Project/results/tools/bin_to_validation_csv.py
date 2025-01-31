import numpy as np
import os
import csv
import pandas as pd

# Load simulation parameters
params = pd.read_csv('results/simulation_parameters.csv')
NX, NY, NZ, NSTEPS, SAVE_EVERY, RE= params.iloc[0][['NX', 'NY', 'NZ', 'NSTEPS', 'NSAVE', 'RE']].astype(int).values
UMAX = params.iloc[0]['U_MAX']

def read_binary_file(filename, shape):
    """
    Legge un file binario in double precision
    e lo rimappa in un array NumPy di forma `shape`,
    assumendo che l'ordine in memoria sia C-order.
    shape = (nz, ny, nx).
    """
    data = np.fromfile(filename, dtype=np.float64)
    return data.reshape(shape)


def extract_z_velocity_line(input_folder, step, nx, ny, nz, output_csv):
    """
    - Legge la componente z della velocità da 'u_z_<step>.bin'
    - Estrae i valori di u_z lungo la linea parallela all'asse x
      che passa per (y = ny//2, z = nz//2).
    - Salva in un CSV con colonne: step, x, y, z, u_z
    """
    # File binario della componente z
    uz_file = os.path.join(input_folder, f"u_z_{step}.bin")
    
    # Carichiamo i dati in un array (nz, ny, nx)
    u_z = read_binary_file(uz_file, (nz, ny, nx))
    
    # Calcoliamo l'indice di y e z al centro
    yc = ny // 2
    zc = nz // 2
    
    # Apriamo un file CSV in scrittura
    with open(output_csv, 'w') as out:
        # Scriviamo un'intestazione (facoltativa)
        out.write("step,x,y,z,u_z\n")
        
        # Cicliamo su x da 0 a nx-1
        for x in range(nx):
            # Leggiamo la velocità u_z in (z=zc, y=yc, x)
            val_uz = u_z[zc, yc, x]
            
            # Scriviamo la riga CSV
            # Esempio di formato: 1000, 12, 32, 32, 0.02345
            out.write(f"{step},{x},{yc},{zc},{val_uz}\n")

def extract_x_velocity_line(input_folder, step, nx, ny, nz, output_csv):
    """
    - Legge la componente z della velocità da 'u_z_<step>.bin'
    - Estrae i valori di u_z lungo la linea parallela all'asse x
      che passa per (y = ny//2, z = nz//2).
    - Salva in un CSV con colonne: step, x, y, z, u_z
    """
    # File binario della componente x
    ux_file = os.path.join(input_folder, f"u_x_{step}.bin")
    
    # Carichiamo i dati in un array (nz, ny, nx)
    u_x = read_binary_file(ux_file, (nz, ny, nx))
    
    # Calcoliamo l'indice di y e z al centro
    xc = nx // 2
    yc = ny // 2
    
    # Apriamo un file CSV in scrittura
    with open(output_csv, 'w') as out:
        # Scriviamo un'intestazione (facoltativa)
        out.write("step,x,y,z,u_x\n")
        
        # Cicliamo su x da 0 a nx-1
        for z in range(nz):
            # Leggiamo la velocità u_x in (x=xc, y=yc, z)
            val_ux = u_x[z, yc, xc]
            
            # Scriviamo la riga CSV
            # Esempio di formato: 1000, 12, 32, 32, 0.02345
            out.write(f"{step},{xc},{yc},{z},{val_ux}\n")


# Esempio d'uso
if __name__ == "__main__":
    step = (NSTEPS // SAVE_EVERY) * SAVE_EVERY
    
    # Cartella dove trovi i file binari
    input_folder = "results/bin_results"
    
    # File CSV di uscita
    if not os.path.exists(f"./validation/{RE}/"):
        os.makedirs(f"./validation/{RE}/")

    output_csv_z = f"./validation/{RE}/u_z_line_center.csv"
    output_csv_x = f"./validation/{RE}/u_x_line_center.csv"
    
    extract_z_velocity_line(input_folder, step, NX, NY, NZ, output_csv_z)
    extract_x_velocity_line(input_folder, step, NX, NY, NZ, output_csv_x)
    print(f"CSV generati! Guarda '{output_csv_z}' e '{output_csv_x}'.")
