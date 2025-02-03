import numpy as np
import matplotlib.pyplot as plt
import imageio
import csv
import pandas as pd
import os

# Load simulation parameters
params = pd.read_csv('results/simulation_parameters.csv')
NX, NY, NZ, NSTEPS, SAVE_EVERY, Re = params.iloc[0][['NX', 'NY', 'NZ', 'NSTEPS', 'NSAVE', 'RE']].astype(int).values
UMAX = params.iloc[0]['U_MAX']
world_size = params.iloc[0]['WORLD_SIZE']

def read_binary_file(filename, shape):
    """Reads a binary file of doubles and maps it into a Numpy array of shape `shape` (C-order)."""
    data = np.fromfile(filename, dtype=np.float64)
    return data.reshape(shape)

def read_partition_u(input_folder, x_start, x_end, step, NZ, NY):
    """
    Reads the u_x, u_y, u_z files for the partition [x_start, x_end) and the requested step.
    Returns three arrays (ux, uy, uz) of shape (NZ, NY, local_nx).
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

def assemble_full_domain(input_folder, x_splits, NX, NY, NZ, step):
    """
    Merges the data of all partitions (divided along x) into three arrays
    (u_x, u_y, u_z) of shape (NZ, NY, NX).
    """
    ux_full = np.zeros((NZ, NY, NX), dtype=np.float64)
    uy_full = np.zeros((NZ, NY, NX), dtype=np.float64)
    uz_full = np.zeros((NZ, NY, NX), dtype=np.float64)
    
    for i in range(len(x_splits) - 1):
        x_start = x_splits[i]
        x_end   = x_splits[i+1]
        
        ux_loc, uy_loc, uz_loc = read_partition_u(input_folder, x_start, x_end, step, NZ, NY)
        ux_full[:, :, x_start:x_end] = ux_loc
        uy_full[:, :, x_start:x_end] = uy_loc
        uz_full[:, :, x_start:x_end] = uz_loc
    
    return ux_full, uy_full, uz_full

def get_partition_points(domain_size, world_size):
    partition_x_size = float(domain_size) / world_size
    partition_points = [round(partition_x_size * i) for i in range(int(world_size) + 1)]
    return partition_points

def create_velocity_gifs_mpi(
    input_folder,
    output_gif,
    NX, NY, NZ,
    x_splits,
    steps,
    slice_y=0,      # fixed y => XZ plane
    fps=5,
    arrow_skip=4,
    arrow_gif_suffix="_arrows"
):
    """
    Creates three GIFs from the simulation on the XZ plane (with fixed y):
     1) normalized velocity map (without arrows)
     2) the same map with arrows (quiver)
     3) the same map with iso-velocity lines (contours)
    """
    frames_no_arrows = []
    frames_arrows    = []
    frames_iso       = []  # For the third GIF (contours)

    # FIGURE WITHOUT ARROWS
    fig_no, ax_no = plt.subplots(figsize=(6,5))
    # FIGURE WITH ARROWS
    fig_ar, ax_ar = plt.subplots(figsize=(6,5))
    # FIGURE WITH ISO-VELOCITY CONTOURS
    fig_iso, ax_iso = plt.subplots(figsize=(6,5))
    
    for step in steps:
        # Reconstruction of the entire domain
        ux_full, uy_full, uz_full = assemble_full_domain(input_folder, x_splits, NX, NY, NZ, step)
        
        # Calculation of normalized velocity
        velocity_magnitude = np.sqrt(ux_full**2 + uy_full**2 + uz_full**2) / UMAX
        
        # Extraction of the XZ plane (shape: (NZ, NX))
        velocity_2d = velocity_magnitude[:, slice_y, :]
        ux_2d = ux_full[:, slice_y, :]
        uz_2d = uz_full[:, slice_y, :]
        
        # --- FIGURE 1: Without arrows ---
        ax_no.clear()
        im_no = ax_no.imshow(
            velocity_2d,
            cmap='coolwarm',
            origin='lower',
            extent=(0, NX, 0, NZ),
            vmin=0,
            vmax=1
        )
        ax_no.set_title(f"Step = {step}, y={slice_y}")
        ax_no.set_xlabel("x")
        ax_no.set_ylabel("z")
        
        if step == steps[0]:
            fig_no.colorbar(im_no, ax=ax_no)
            fig_no.suptitle("Velocity field normalized to u_max", fontsize=14)
        
        fig_no.canvas.draw()
        buf_no = fig_no.canvas.buffer_rgba()
        image_no = np.frombuffer(buf_no, dtype=np.uint8)
        image_no = image_no.reshape(fig_no.canvas.get_width_height()[::-1] + (4,))
        image_no = image_no[..., :3]
        image_no = image_no.copy()  # Copy to separate memory
        
        frames_no_arrows.append(image_no)
        
        # --- FIGURE 2: With arrows (quiver) ---
        ax_ar.clear()
        im_ar = ax_ar.imshow(
            velocity_2d,
            cmap='coolwarm',
            origin='lower',
            extent=(0, NX, 0, NZ),
            vmin=0,
            vmax=1
        )
        ax_ar.set_title(f"Step = {step}, y={slice_y}")
        ax_ar.set_xlabel("x")
        ax_ar.set_ylabel("z")

        if step == steps[0]:
            fig_ar.colorbar(im_ar, ax=ax_ar)
            fig_ar.suptitle("Velocity field normalized to u_max", fontsize=14)
        
        # Creation of arrows (quiver) => (u_x, u_z) in (x,z)
        arrow_x = np.arange(0, NX, arrow_skip)
        arrow_z = np.arange(0, NZ, arrow_skip)
        Xarrow, Zarrow = np.meshgrid(arrow_x, arrow_z)
        
        Uarrow = ux_2d[Zarrow, Xarrow]
        Varrow = uz_2d[Zarrow, Xarrow]
        
        ax_ar.quiver(
            Xarrow, Zarrow,
            Uarrow, Varrow,
            color='white',
            scale=2.0
        )
        
        fig_ar.canvas.draw()
        buf_ar = fig_ar.canvas.buffer_rgba()
        image_ar = np.frombuffer(buf_ar, dtype=np.uint8)
        image_ar = image_ar.reshape(fig_ar.canvas.get_width_height()[::-1] + (4,))
        image_ar = image_ar[..., :3]
        image_ar = image_ar.copy()  # Copy to separate memory
        
        frames_arrows.append(image_ar)
        
        # --- FIGURE 3: With iso-velocity lines (contours) ---
        ax_iso.clear()
        im_iso = ax_iso.imshow(
            velocity_2d,
            cmap='Greys',
            origin='lower',
            extent=(0, NX, 0, NZ),
            vmin=0,
            vmax=1
        )
        # Construct x and z coordinates consistent with extent
        x_coords = np.linspace(0, NX, velocity_2d.shape[1])
        z_coords = np.linspace(0, NZ, velocity_2d.shape[0])
        # Define contour levels (iso-velocity)
        levels = np.linspace(0, 1, 10)
        CS = ax_iso.contour(
            x_coords,
            z_coords,
            velocity_2d,
            levels=levels,
            colors='red',
            linewidths=1.5
        )
        # Add labels to contours
        ax_iso.clabel(CS, inline=True, fontsize=8)
        ax_iso.set_title(f"Step = {step}, y={slice_y} (Iso-velocity)")
        ax_iso.set_xlabel("x")
        ax_iso.set_ylabel("z")
        
        if step == steps[0]:
            fig_iso.colorbar(im_iso, ax=ax_iso)
            fig_iso.suptitle("Velocity field normalized to u_max (Iso-velocity lines)", fontsize=14)
        
        fig_iso.canvas.draw()
        buf_iso = fig_iso.canvas.buffer_rgba()
        image_iso = np.frombuffer(buf_iso, dtype=np.uint8)
        image_iso = image_iso.reshape(fig_iso.canvas.get_width_height()[::-1] + (4,))
        image_iso = image_iso[..., :3]
        image_iso = image_iso.copy()  # Copy to separate memory
        
        frames_iso.append(image_iso)
    
    plt.close(fig_no)
    plt.close(fig_ar)
    plt.close(fig_iso)

    n = len(frames_arrows)
    frame_duration = 1000.0 / fps
    durations = [frame_duration] * (n - 1)
    durations.append(2000.0)  # Duration of the last frame = 2s

    # GIF 1 (without arrows)
    imageio.mimsave(output_gif, frames_no_arrows, duration=durations, loop=0)
    print(f"GIF saved (without arrows): {output_gif}")
    
    # GIF 2 (with arrows)
    base, ext = os.path.splitext(output_gif)
    output_gif_arrows = base + arrow_gif_suffix + ext
    imageio.mimsave(output_gif_arrows, frames_arrows, duration=durations, loop=0)
    print(f"GIF saved (with arrows): {output_gif_arrows}")
    
    # GIF 3 (with iso-velocity lines)
    output_gif_iso = base + "_iso" + ext
    imageio.mimsave(output_gif_iso, frames_iso, duration=durations, loop=0)
    print(f"GIF saved (with iso-velocity lines): {output_gif_iso}")

# Example usage
if __name__ == "__main__":
    steps = range(0, NSTEPS, SAVE_EVERY)
    x_splits = get_partition_points(NX, world_size)
    
    input_folder = "results/bin_results"
    output_gif   = "./vel_field.gif"
    
    slice_y    = 32
    fps        = 10
    arrow_skip = 4
    
    create_velocity_gifs_mpi(
        input_folder=input_folder,
        output_gif=output_gif,
        NX=NX, NY=NY, NZ=NZ,
        x_splits=x_splits,
        steps=steps,
        slice_y=slice_y,
        fps=fps,
        arrow_skip=arrow_skip
    )
