import numpy as np
import os

def load_field(filename, nx, ny):
    """Load binary field file."""
    if os.path.exists(filename):
        data = np.fromfile(filename, dtype=np.float64)
        if data.size == nx * ny:
            print(f"Loaded {filename}: min={data.min()}, max={data.max()}")
            return data.reshape((ny, nx))
        else:
            print(f"Error: File {filename} has incorrect size {data.size}, expected {nx * ny}")
    else:
        print(f"File {filename} not found.")
    return np.zeros((ny, nx))

def save_to_vtk(filename, data, field_name, nx, ny):
    """Save a 2D field to a VTK file."""
    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"{field_name} field\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx} {ny} 1\n")
        f.write("ORIGIN 0 0 0\n")
        f.write("SPACING 1 1 1\n")
        f.write(f"POINT_DATA {nx * ny}\n")
        f.write(f"SCALARS {field_name} double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(ny):
            for i in range(nx):
                f.write(f"{data[j, i]}\n")
        print(f"Saved {filename}")

def convert_bin_to_vtk(bin_dir, vtk_dir, field_names, nx, ny):
    """Convert binary field files to VTK files."""
    os.makedirs(vtk_dir, exist_ok=True)
    for field in field_names:
        for timestep in sorted(os.listdir(bin_dir)):
            if timestep.startswith(field) and timestep.endswith('.bin'):
                bin_file = os.path.join(bin_dir, timestep)
                vtk_file = os.path.join(vtk_dir, timestep.replace('.bin', '.vtk'))
                data = load_field(bin_file, nx, ny)
                save_to_vtk(vtk_file, data, field, nx, ny)

# Example usage
if __name__ == "__main__":
    BIN_DIR = "../bin_results"
    VTK_DIR = "../vtk_results"
    FIELD_NAMES = ["rho", "ux", "uy"]
    import pandas as pd

# Load simulation parameters
params = pd.read_csv('simulation_parameters.csv')
NX, NY = params.iloc[0][['NX', 'NY']].astype(int).values  # Example dimensions

convert_bin_to_vtk(BIN_DIR, VTK_DIR, FIELD_NAMES, NX, NY)
