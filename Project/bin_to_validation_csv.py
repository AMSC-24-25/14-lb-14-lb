import numpy as np
import os
import csv
import pandas as pd

# Load simulation parameters
params = pd.read_csv('simulation_parameters.csv')
NX, NY, NZ, NSTEPS, SAVE_EVERY, RE= params.iloc[0][['NX', 'NY', 'NZ', 'NSTEPS', 'NSAVE', 'RE']].astype(int).values
UMAX = params.iloc[0]['U_MAX']

def read_binary_file(filename, shape):
    """
    Reads a binary file in double precision
    and maps it into a NumPy array of shape `shape`,
    assuming that the memory order is C-order.
    shape = (nz, ny, nx).
    """
    data = np.fromfile(filename, dtype=np.float64)
    return data.reshape(shape)


def extract_z_velocity_line(input_folder, step, nx, ny, nz, output_csv):
    """
    - Reads the z component of the velocity from 'u_z_<step>.bin'
    - Extracts the values of u_z along the line parallel to the x-axis
      that passes through (y = ny//2, z = nz//2).
    - Saves in a CSV with columns: step, x, y, z, u_z
    """
    # Binary file of the z component
    uz_file = os.path.join(input_folder, f"u_z_{step}.bin")
    
    # Load the data into an array (nz, ny, nx)
    u_z = read_binary_file(uz_file, (nz, ny, nx))
    
    # Calculate the index of y and z at the center
    yc = ny // 2
    zc = nz // 2
    
    # Open a CSV file for writing
    with open(output_csv, 'w') as out:
        # Write a header (optional)
        out.write("step,x,y,z,u_z\n")
        
        # Loop over x from 0 to nx-1
        for x in range(nx):
            # Read the velocity u_z at (z=zc, y=yc, x)
            val_uz = u_z[zc, yc, x]
            
            # Write the CSV row
            # Example format: 1000, 12, 32, 32, 0.02345
            out.write(f"{step},{x},{yc},{zc},{val_uz}\n")

def extract_x_velocity_line(input_folder, step, nx, ny, nz, output_csv):
    """
    - Reads the x component of the velocity from 'u_x_<step>.bin'
    - Extracts the values of u_x along the line parallel to the x-axis
      that passes through (y = ny//2, z = nz//2).
    - Saves in a CSV with columns: step, x, y, z, u_x
    """
    # Binary file of the x component
    ux_file = os.path.join(input_folder, f"u_x_{step}.bin")
    
    # Load the data into an array (nz, ny, nx)
    u_x = read_binary_file(ux_file, (nz, ny, nx))
    
    # Calculate the index of y and z at the center
    xc = nx // 2
    yc = ny // 2
    
    # Open a CSV file for writing
    with open(output_csv, 'w') as out:
        # Write a header (optional)
        out.write("step,x,y,z,u_x\n")
        
        # Loop over z from 0 to nz-1
        for z in range(nz):
            # Read the velocity u_x at (x=xc, y=yc, z)
            val_ux = u_x[z, yc, xc]
            
            # Write the CSV row
            # Example format: 1000, 12, 32, 32, 0.02345
            out.write(f"{step},{xc},{yc},{z},{val_ux}\n")


# Example usage
if __name__ == "__main__":
    step = (NSTEPS // SAVE_EVERY) * SAVE_EVERY
    
    # Folder where the binary files are located
    input_folder = "./bin_results"
    
    # Output CSV file
    if not os.path.exists(f"./validation/{RE}/"):
        os.makedirs(f"./validation/{RE}/")

    output_csv_z = f"./validation/{RE}/u_z_line_center.csv"
    output_csv_x = f"./validation/{RE}/u_x_line_center.csv"
    
    extract_z_velocity_line(input_folder, step, NX, NY, NZ, output_csv_z)
    extract_x_velocity_line(input_folder, step, NX, NY, NZ, output_csv_x)
    print(f"CSV files generated! Check '{output_csv_z}' and '{output_csv_x}'.")
