#!/u/sw/toolchains/gcc-glibc/11.2.0/base/bin/python
import numpy as np
import os

# Funzione per leggere un file binario
def read_binary_file(filename, shape):
    """Legge un file binario e lo converte in un array numpy."""
    data = np.fromfile(filename, dtype=np.float64)
    return data.reshape(shape)

# Funzione per scrivere i file VTK

def write_vtk(output_file, velocity_magnitude, nx, ny, nz):
    """Scrive i dati del modulo della velocità in formato VTK per Paraview."""
    with open(output_file, 'w') as vtk_file:
        # Header VTK
        vtk_file.write("# vtk DataFile Version 3.0\n")
        vtk_file.write("Velocity magnitude field\n")
        vtk_file.write("ASCII\n")
        vtk_file.write("DATASET STRUCTURED_POINTS\n")
        vtk_file.write(f"DIMENSIONS {nx} {ny} {nz}\n")
        vtk_file.write("ORIGIN 0 0 0\n")
        vtk_file.write("SPACING 1 1 1\n")
        vtk_file.write(f"POINT_DATA {nx * ny * nz}\n")

        # Scrivi il modulo della velocità
        vtk_file.write("SCALARS velocity_magnitude float\n")
        vtk_file.write("LOOKUP_TABLE default\n")
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    vtk_file.write(f"{velocity_magnitude[z, y, x]}\n")

# Funzione principale
def convert_bin_to_vtk(input_folder, output_folder, nx, ny, nz, steps):
    """Converte i file binari in file VTK per Paraview."""
    for step in steps:
        # File binari per ogni componente della velocità
        ux_file = os.path.join(input_folder, f"u_x_{step}.bin")
        uy_file = os.path.join(input_folder, f"u_y_{step}.bin")
        uz_file = os.path.join(input_folder, f"u_z_{step}.bin")

        # Leggi i file binari
        u_x = read_binary_file(ux_file, (nz, ny, nx))
        u_y = read_binary_file(uy_file, (nz, ny, nx))
        u_z = read_binary_file(uz_file, (nz, ny, nx))

        # Calcola il modulo della velocità
        velocity_magnitude = np.sqrt(u_x**2 + u_y**2 + u_z**2)

        # Scrivi il file VTK
        vtk_filename = os.path.join(output_folder, f"velocity_magnitude_{step}.vtk")
        write_vtk(vtk_filename, velocity_magnitude, nx, ny, nz)
        print(f"File VTK creato: {vtk_filename}")

# Parametri del dominio e del simulatore
input_folder = "./bin_results"  # Cartella dei risultati binari
output_folder = "./vtk_results"  # Cartella dei risultati binari
nx, ny, nz = 128, 128, 128  # Dimensioni del dominio
steps = list(range(18000, 20000, 40))  # Passi temporali salvati

# Converte i file binari in formato VTK
convert_bin_to_vtk(input_folder, output_folder, nx, ny, nz, steps)
