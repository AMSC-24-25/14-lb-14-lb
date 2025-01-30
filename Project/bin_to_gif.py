import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import imageio
import os

# Carichiamo i parametri da CSV (come facevi già)
params = pd.read_csv('simulation_parameters.csv')
NX, NY, NZ, NSTEPS, SAVE_EVERY, Re = params.iloc[0][['NX', 'NY', 'NZ', 'NSTEPS', 'NSAVE', 'RE']].astype(int).values
UMAX = params.iloc[0]['U_MAX']

def read_binary_file(filename, shape):
    """Legge un file binario con double e lo converte in un array numpy di forma 'shape'."""
    data = np.fromfile(filename, dtype=np.float64)
    return data.reshape(shape)

def create_velocity_gif_with_arrows(
    input_folder,
    output_gif,
    nx, ny, nz,
    steps,
    slice_y=0,
    fps=5,
    arrow_skip=4,
    arrow_gif_suffix="_arrows"
):
    """
    Crea DUE GIF:
      1) La mappa 2D del modulo della velocità
      2) La stessa mappa 2D con sovrapposto un quiver (frecce) ogni 'arrow_skip' celle.

    Parametri aggiuntivi rispetto alla versione base:
    - arrow_skip: ogni quanti punti tracciare le frecce
    - arrow_gif_suffix: suffisso da aggiungere al nome della seconda GIF
    """
    # Memorizziamo i frame di entrambe le GIF
    frames_no_arrows = []
    frames_arrows = []

    # Creiamo due figure separate: una per la mappa liscia, una per la mappa con frecce
    fig_no, ax_no = plt.subplots(figsize=(6,5))
    fig_ar, ax_ar = plt.subplots(figsize=(6,5))

    # Loop su tutti gli step da salvare
    for step in steps:
        # Costruiamo i nomi dei file binari per le 3 componenti
        ux_file = os.path.join(input_folder, f"u_x_{step}.bin")
        uy_file = os.path.join(input_folder, f"u_y_{step}.bin")
        uz_file = os.path.join(input_folder, f"u_z_{step}.bin")
        
        # Leggiamo i file (array di shape (nz, ny, nx))
        u_x = read_binary_file(ux_file, (nz, ny, nx))
        u_y = read_binary_file(uy_file, (nz, ny, nx))
        u_z = read_binary_file(uz_file, (nz, ny, nx))
        
        # Calcoliamo il modulo della velocità: shape (nz, ny, nx)
        velocity_magnitude = np.sqrt(u_x**2 + u_y**2 + u_z**2)
        
        # Estraiamo la "fetta" 2D con dimensioni (nz, nx), 
        # facendo: velocity_2d = velocity_magnitude[:, slice_y, :] 
        # (coerente con il tuo codice originale)
        velocity_2d = velocity_magnitude[:, slice_y, :]  # shape (nz, nx)

        # Prepariamo anche le 2D delle singole componenti (x e z) se vogliamo il quiver
        # supponendo che nel piano (z,x) ci interessino (u_z, u_x) come coordinate delle frecce
        # Oppure vuoi (u_x, u_y)? Dipende dal tuo scopo. Qui ipotizziamo (u_x, u_z).
        ux_2d = u_x[:, slice_y, :]  # shape (nz, nx)
        uz_2d = u_z[:, slice_y, :]  # shape (nz, nx)

        ## ----------------- PRIMA FIGURA: senza frecce ----------------------
        ax_no.clear()
        # Mostriamo la mappa 2D - in origin='lower' stiamo interpretando l'indice 0 come "in basso"
        im_no = ax_no.imshow(
            velocity_2d,
            cmap='jet',
            origin='lower',
            extent=(0, nx, 0, nz),  # asse orizzontale: x in [0..nx], verticale: z in [0..nz]
            vmin=0,
            vmax=velocity_2d.max()
        )
        ax_no.set_title(f"Step = {step}, slice_y = {slice_y} (no arrows)")
        ax_no.set_xlabel("x")
        ax_no.set_ylabel("z")

        # Convertiamo la figura in un'immagine e la salviamo nei frames_no_arrows
        fig_no.canvas.draw()
        image_no = np.frombuffer(fig_no.canvas.tostring_rgb(), dtype='uint8')
        image_no = image_no.reshape(fig_no.canvas.get_width_height()[::-1] + (3,))
        frames_no_arrows.append(image_no)

        ## ----------------- SECONDA FIGURA: con frecce ----------------------
        ax_ar.clear()
        # Mostriamo la stessa mappa
        im_ar = ax_ar.imshow(
            velocity_2d,
            cmap='jet',
            origin='lower',
            extent=(0, nx, 0, nz),
            vmin=0,
            vmax=velocity_2d.max()
        )
        ax_ar.set_title(f"Step = {step}, slice_y = {slice_y} (+ arrows)")
        ax_ar.set_xlabel("x")
        ax_ar.set_ylabel("z")

        # Costruiamo la griglia di campionamento per le frecce
        # arrow_skip controlla ogni quanti punti prendiamo in x e z
        arrow_x = np.arange(0, nx, arrow_skip)
        arrow_z = np.arange(0, nz, arrow_skip)
        Xarrow, Zarrow = np.meshgrid(arrow_x, arrow_z)  # shape (len(arrow_z), len(arrow_x))

        # Estraggo i valori corrispondenti di (u_x, u_z)
        # occhio che l'indice riga corrisponde a z, la colonna a x
        Uarrow = ux_2d[Zarrow, Xarrow]
        Varrow = uz_2d[Zarrow, Xarrow]

        # Disegniamo le frecce. Scale/regole di scaling da adattare a piacere.
        # color='white' per contrastare con la colormap 'jet'.
        ax_ar.quiver(
            Xarrow, Zarrow,
            Uarrow, Varrow,
            color='white',
            scale=0.8  # puoi regolare scale a tuo piacimento
        )

        # Convertiamo la figura con frecce in un'immagine e la salviamo
        fig_ar.canvas.draw()
        image_ar = np.frombuffer(fig_ar.canvas.tostring_rgb(), dtype='uint8')
        image_ar = image_ar.reshape(fig_ar.canvas.get_width_height()[::-1] + (3,))
        frames_arrows.append(image_ar)

    # Fine del loop su steps: chiudiamo le figure
    plt.close(fig_no)
    plt.close(fig_ar)

    # Salviamo la prima GIF (senza frecce)
    imageio.mimsave(output_gif, frames_no_arrows, fps=fps)
    print(f"GIF salvata (senza frecce): {output_gif}")

    # Salviamo la seconda GIF (con frecce)
    base, ext = os.path.splitext(output_gif)
    output_gif_arrows = base + arrow_gif_suffix + ext
    imageio.mimsave(output_gif_arrows, frames_arrows, fps=fps)
    print(f"GIF salvata (con frecce): {output_gif_arrows}")


# ESEMPIO di utilizzo
if __name__ == "__main__":
    input_folder = "./bin_results"
    output_gif = "vel_field.gif"
    steps = range(0, NSTEPS, SAVE_EVERY)
    slice_y = 64
    fps = 5
    arrow_skip = 6  # ogni 4 celle disegniamo una freccia

    create_velocity_gif_with_arrows(
        input_folder=input_folder,
        output_gif=output_gif,
        nx=NX, ny=NY, nz=NZ,
        steps=steps,
        slice_y=slice_y,
        fps=fps,
        arrow_skip=arrow_skip
    )
