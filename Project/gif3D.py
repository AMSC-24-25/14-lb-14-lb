import numpy as np
import matplotlib.pyplot as plt
import imageio
import os

def read_binary_file(filename, shape):
    """Legge un file binario di double e lo rimappa in un array Numpy di forma `shape` (C-order)."""
    data = np.fromfile(filename, dtype=np.float64)
    return data.reshape(shape)

def read_partition_u(input_folder, x_start, x_end, step, NZ, NY):
    """
    Legge i file u_x, u_y, u_z per la partizione [x_start, x_end) e lo step richiesto.
    Ritorna tre array (ux, uy, uz) di shape (NZ, NY, local_nx).
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
    Unisce i dati di tutte le partizioni (divise lungo x) in tre array
    (u_x, u_y, u_z) di shape (NZ, NY, NX).
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

def create_velocity_gifs_mpi(
    input_folder,
    output_gif,
    NX, NY, NZ,
    x_splits,
    steps,
    slice_y=0,      # fisso y => piano XZ
    fps=5,
    arrow_skip=4,
    arrow_gif_suffix="_arrows"
):
    """
    Crea due GIF (3D → estrazione piano XZ a y=slice_y):
     1) sola mappa 2D del modulo velocità
     2) la stessa mappa con frecce (quiver).
    """
    frames_no_arrows = []
    frames_arrows = []
    
    fig_no, ax_no = plt.subplots(figsize=(6,5))
    fig_ar, ax_ar = plt.subplots(figsize=(6,5))
    
    for step in steps:
        # Ricostruisci l'intero dominio
        ux_full, uy_full, uz_full = assemble_full_domain(input_folder, x_splits, NX, NY, NZ, step)
        
        # Modulo velocità
        velocity_magnitude = np.sqrt(ux_full**2 + uy_full**2 + uz_full**2)
        
        # Estrazione piano XZ (shape: (NZ, NX))
        velocity_2d = velocity_magnitude[:, slice_y, :]
        ux_2d = ux_full[:, slice_y, :]
        uz_2d = uz_full[:, slice_y, :]
        
        # --- FIGURA 1: no arrows ---
        ax_no.clear()
        im_no = ax_no.imshow(
            velocity_2d,
            cmap='jet',
            origin='lower',
            extent=(0, NX, 0, NZ),
            vmin=0,
            vmax=velocity_2d.max()
        )
        ax_no.set_title(f"Step = {step}, y={slice_y} (no arrows)")
        ax_no.set_xlabel("x")
        ax_no.set_ylabel("z")
        
        fig_no.canvas.draw()
        # Prima usavi fig.canvas.tostring_rgb(), ora usi buffer_rgba()
        buf_no = fig_no.canvas.buffer_rgba()
        # Converti in array NumPy
        image_no = np.frombuffer(buf_no, dtype=np.uint8)
        # buffer_rgba produce un array con shape (h, w, 4)
        image_no = image_no.reshape(fig_no.canvas.get_width_height()[::-1] + (4,))
        # Se vuoi rimuovere l'alpha (ultime dimensioni) e tenere solo RGB:
        image_no = image_no[..., :3]
        image_no = image_no.copy()  # Copia in memoria separata
        
        frames_no_arrows.append(image_no)
        
        # --- FIGURA 2: con arrows ---
        ax_ar.clear()
        im_ar = ax_ar.imshow(
            velocity_2d,
            cmap='jet',
            origin='lower',
            extent=(0, NX, 0, NZ),
            vmin=0,
            vmax=velocity_2d.max()
        )
        ax_ar.set_title(f"Step = {step}, y={slice_y} (+ arrows)")
        ax_ar.set_xlabel("x")
        ax_ar.set_ylabel("z")
        
        # Quiver => (u_x, u_z) in (x,z)
        arrow_x = np.arange(0, NX, arrow_skip)
        arrow_z = np.arange(0, NZ, arrow_skip)
        Xarrow, Zarrow = np.meshgrid(arrow_x, arrow_z)
        
        Uarrow = ux_2d[Zarrow, Xarrow]
        Varrow = uz_2d[Zarrow, Xarrow]
        
        ax_ar.quiver(
            Xarrow, Zarrow,
            Uarrow, Varrow,
            color='white',
            scale=20
        )
        
        fig_ar.canvas.draw()
        buf_ar = fig_ar.canvas.buffer_rgba()
        image_ar = np.frombuffer(buf_ar, dtype=np.uint8)
        image_ar = image_ar.reshape(fig_ar.canvas.get_width_height()[::-1] + (4,))
        image_ar = image_ar[..., :3]
        image_ar = image_ar.copy()  # Copia in memoria separata
        
        frames_arrows.append(image_ar)
    
    plt.close(fig_no)
    plt.close(fig_ar)
    
    # GIF 1 (senza frecce)
    imageio.mimsave(output_gif, frames_no_arrows, fps=fps)
    print(f"GIF salvata (senza frecce): {output_gif}")
    
    # GIF 2 (con frecce)
    base, ext = os.path.splitext(output_gif)
    output_gif_arrows = base + arrow_gif_suffix + ext
    imageio.mimsave(output_gif_arrows, frames_arrows, fps=fps)
    print(f"GIF salvata (con frecce): {output_gif_arrows}")

# Esempio di utilizzo
if __name__ == "__main__":
    # Parametri esempio
    #x_splits = [0, 128]
    x_splits = [0, 16, 32, 48, 64, 80, 96, 112, 128]
    NX, NY, NZ = 128, 128, 128
    steps = range(0, 121, 40)
    
    input_folder = "./bin_results"
    output_gif   = "vel_field.gif"
    
    slice_y    = 32
    fps        = 5
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
