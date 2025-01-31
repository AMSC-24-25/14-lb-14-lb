import numpy as np
import csv
import pandas as pd
import os

# Load simulation parameters
params = pd.read_csv('results/simulation_parameters.csv')
NX, NY, NZ, NSTEPS, SAVE_EVERY, Re = params.iloc[0][['NX', 'NY', 'NZ', 'NSTEPS', 'NSAVE', 'RE']].astype(int).values
UMAX = params.iloc[0]['U_MAX']
world_size = params.iloc[0]['WORLD_SIZE']

def read_binary_file(filename, shape):
    """
    Reads a binary file (double) and maps it into a Numpy array
    with shape `shape`, assuming 'C-order' ordering.
    """
    data = np.fromfile(filename, dtype=np.float64)
    return data.reshape(shape)

def read_partition_3d(input_folder, x_start, x_end, step, NZ, NY):
    """
    Reads the files u_x, u_y, u_z for the partition [x_start, x_end) at step `step`.
    Returns (u_x_local, u_y_local, u_z_local), with shape (NZ, NY, local_nx).
    """
    local_nx = x_end - x_start
    
    part_dir = os.path.join(input_folder, f"partition_{x_start}-{x_end}")
    
    ux_file = os.path.join(part_dir, f"u_x_{step}.bin")
    uy_file = os.path.join(part_dir, f"u_y_{step}.bin")
    uz_file = os.path.join(part_dir, f"u_z_{step}.bin")

    ux_local = read_binary_file(ux_file, (NZ, NY, local_nx))
    uy_local = read_binary_file(uy_file, (NZ, NY, local_nx))
    uz_local = read_binary_file(uz_file, (NZ, NY, local_nx))

    return ux_local, uy_local, uz_local

def assemble_full_domain_3d(input_folder, x_splits, NX, NY, NZ, step):
    """
    Reconstructs the fields (u_x, u_y, u_z) of the entire domain (NZ, NY, NX)
    by merging the partition data.
    """
    ux_full = np.zeros((NZ, NY, NX), dtype=np.float64)
    uy_full = np.zeros((NZ, NY, NX), dtype=np.float64)
    uz_full = np.zeros((NZ, NY, NX), dtype=np.float64)
    
    for i in range(len(x_splits) - 1):
        x_start = x_splits[i]
        x_end   = x_splits[i+1]
        
        ux_loc, uy_loc, uz_loc = read_partition_3d(input_folder, x_start, x_end, step, NZ, NY)
        ux_full[:, :, x_start:x_end] = ux_loc
        uy_full[:, :, x_start:x_end] = uy_loc
        uz_full[:, :, x_start:x_end] = uz_loc
    
    return ux_full, uy_full, uz_full

def get_partition_points(domain_size, world_size):
    partition_x_size = float(domain_size) / world_size
    partition_points = [round(partition_x_size * i) for i in range(int(world_size) + 1)]
    return partition_points

def write_vtk_scalar_speed(filename, ux, uy, uz, NX, NY, NZ):
    """
    Writes a VTK file (legacy format, ASCII) with a SCALAR field
    'velocity_magnitude', calculated as sqrt(u_x^2 + u_y^2 + u_z^2).

    The grid is 'STRUCTURED_POINTS', dimensions NX, NY, NZ.
    We will have POINT_DATA = NX*NY*NZ and a "SCALARS velocity_magnitude float".
    """
    npoints = NX * NY * NZ
    
    with open(filename, 'w') as vtk_file:
        # VTK Header
        vtk_file.write("# vtk DataFile Version 3.0\n")
        vtk_file.write("LBM speed field\n")
        vtk_file.write("ASCII\n")
        vtk_file.write("DATASET STRUCTURED_POINTS\n")
        vtk_file.write(f"DIMENSIONS {NX} {NY} {NZ}\n")
        vtk_file.write("ORIGIN 0 0 0\n")
        vtk_file.write("SPACING 1 1 1\n")
        vtk_file.write(f"POINT_DATA {npoints}\n")
        
        # Declaration of a scalar field
        vtk_file.write("SCALARS velocity_magnitude float 1\n")
        vtk_file.write("LOOKUP_TABLE default\n")

        # Iterate over nodes (z,y,x) for compatibility with shape (NZ, NY, NX)
        for z in range(NZ):
            for y in range(NY):
                for x in range(NX):
                    vx = ux[z, y, x]
                    vy = uy[z, y, x]
                    vz = uz[z, y, x]
                    speed = np.sqrt(vx*vx + vy*vy + vz*vz)
                    
                    # Write one scalar value per line
                    vtk_file.write(f"{speed}\n")

def export_speed_to_vtk_3d(
    input_folder,
    output_folder,
    x_splits,
    NX, NY, NZ,
    steps
):
    """
    For each step, reconstructs the entire domain and saves a .vtk file
    """
    os.makedirs(output_folder, exist_ok=True)
    
    for step in steps:
        # Reconstruct the velocity fields
        ux, uy, uz = assemble_full_domain_3d(input_folder, x_splits, NX, NY, NZ, step)
        
        # Generate the .vtk file name
        out_vtk = os.path.join(output_folder, f"speed_{step}.vtk")
        
        # Write the file with SCALARS
        write_vtk_scalar_speed(out_vtk, ux, uy, uz, NX, NY, NZ)
        print(f"Created VTK file: {out_vtk}")

# Example usage
if __name__ == "__main__":
    steps = range(0, NSTEPS, SAVE_EVERY)
    x_splits = get_partition_points(NX, world_size)
    
    input_folder = "results/bin_results"
    output_folder = "results/vtk_results"
    
    export_speed_to_vtk_3d(
        input_folder=input_folder,
        output_folder=output_folder,
        x_splits=x_splits,
        NX=NX, NY=NY, NZ=NZ,
        steps=steps
    )
