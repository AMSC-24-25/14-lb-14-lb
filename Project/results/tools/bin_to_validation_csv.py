import numpy as np
import os
import csv
import pandas as pd

# Carichiamo i parametri di simulazione
params = pd.read_csv('./results/simulation_parameters.csv')
NX, NY, NZ, NSTEPS, SAVE_EVERY, RE= params.iloc[0][['NX', 'NY', 'NZ', 'NSTEPS', 'NSAVE', 'RE']].astype(int).values
UMAX = params.iloc[0]['U_MAX']
world_size = params.iloc[0]['WORLD_SIZE']

# Esempio di definizione partizionamento lungo x.
# Se NON sai a priori i confini, puoi impostarli a mano, oppure
# eseguire uno script per leggere le subdirectory 'partition_*-*' e dedurre i x_splits.
# Qui assumiamo un esempio con 2 partizioni: [0, NX//2, NX].
def read_binary_file(filename, shape):
    """
    Legge un file binario in double precision
    e lo rimappa in un array NumPy di forma `shape`,
    assumendo che l'ordine in memoria sia C-order.
    shape = (nz, ny, local_nx).
    """
    data = np.fromfile(filename, dtype=np.float64)
    return data.reshape(shape)

def get_partition_points(domain_size, world_size):
    partition_x_size = float(domain_size) / world_size
    partition_points = [round(partition_x_size * i) for i in range(int(world_size) + 1)]
    return partition_points

def read_partition_variable(
    input_folder, x_start, x_end, step, NZ, NY, varName
):
    """
    Legge il file binario 'u_<varName>_<step>.bin' per la partizione [x_start, x_end),
    e lo restituisce come array di shape (NZ, NY, local_nx).
    Esempio: varName = 'z' => 'u_z_<step>.bin'
    """
    local_nx = x_end - x_start
    part_dir = os.path.join(input_folder, f"partition_{x_start}-{x_end}")
    filename = os.path.join(part_dir, f"u_{varName}_{step}.bin")
    
    # Leggiamo i dati
    arr_local = read_binary_file(filename, (NZ, NY, local_nx))
    return arr_local

def assemble_full_variable(
    input_folder, x_splits, NX, NY, NZ, step, varName
):
    """
    Ricostruisce la variabile 'u_<varName>' sull'intero dominio (NZ, NY, NX),
    unendo i dati dalle varie partizioni (lungo x).
    Ritorna un array di shape (NZ, NY, NX).
    """
    full_arr = np.zeros((NZ, NY, NX), dtype=np.float64)
    
    for i in range(len(x_splits) - 1):
        x_start = x_splits[i]
        x_end   = x_splits[i+1]
        
        part_local = read_partition_variable(input_folder, x_start, x_end, step, NZ, NY, varName)
        # Inseriamo nella porzione x_start:x_end
        full_arr[:, :, x_start:x_end] = part_local
    
    return full_arr

def extract_z_velocity_line(input_folder, step, nx, ny, nz, output_csv):
    """
    - Ricostruisce la componente z della velocità unendo le partizioni
    - Estrae i valori di u_z lungo la linea parallela all'asse x,
      che passa per (y=ny//2, z=nz//2).
    - Salva in un CSV con colonne: step,x,y,z,u_z
    """
    # Ricostruisci u_z (forma (nz, ny, nx))
    u_z = assemble_full_variable(input_folder, x_splits, nx, ny, nz, step, varName='z')
    
    # Calcoliamo l'indice di y e z al centro
    yc = ny // 2
    zc = nz // 2
    
    # Apriamo un file CSV in scrittura
    with open(output_csv, 'w') as out:
        out.write("step,x,y,z,u_z\n")
        # Cicliamo su x
        for x in range(nx):
            val_uz = u_z[zc, yc, x]
            out.write(f"{step},{x},{yc},{zc},{val_uz}\n")

def extract_x_velocity_line(input_folder, step, nx, ny, nz, output_csv):
    """
    - Ricostruisce la componente x della velocità unendo le partizioni
    - Estrae i valori di u_x lungo la linea parallela all'asse z,
      che passa per (y=ny//2, x=nx//2). (Nota: sto interpretando "parallel=asse x"
      in modo coerente con la funzione org, ma sembra più un "parallel=asse z" data la descrizione)
    - Salva in un CSV con colonne: step,x,y,z,u_x
    """
    # Ricostruisci u_x (forma (nz, ny, nx))
    u_x = assemble_full_variable(input_folder, x_splits, nx, ny, nz, step, varName='x')
    
    # Indici centrali
    xc = nx // 2
    yc = ny // 2
    
    # Salviamo i valori scorrendo z in [0..nz-1], e x=xc, y=yc fissi
    with open(output_csv, 'w') as out:
        out.write("step,x,y,z,u_x\n")
        for z in range(nz):
            val_ux = u_x[z, yc, xc]
            out.write(f"{step},{xc},{yc},{z},{val_ux}\n")

# Esempio d'uso
if __name__ == "__main__":
    # Scegli lo step da cui estrarre i dati
    step = (NSTEPS // SAVE_EVERY) * SAVE_EVERY
    x_splits = get_partition_points(NX, world_size)
    
    # Cartella dove trovi i file binari MPI
    input_folder = "./results/bin_results"
    
    # Cartella di output per i CSV
    os.makedirs(f"./validation/{RE}/", exist_ok=True)
    
    output_csv_z = f"./validation/{RE}/u_z_line_center.csv"
    output_csv_x = f"./validation/{RE}/u_x_line_center.csv"
    
    extract_z_velocity_line(input_folder, step, NX, NY, NZ, output_csv_z)
    extract_x_velocity_line(input_folder, step, NX, NY, NZ, output_csv_x)
    
    print(f"CSV generati! Guarda '{output_csv_z}' e '{output_csv_x}'.")
