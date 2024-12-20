# Lattice Boltzmann Method (LBM) Simulation

## Overview
This project implements a Lattice Boltzmann Method (LBM) simulation, a numerical method widely used for fluid dynamics and related computational physics applications. The project is written in C++ and includes core functionality to initialize the simulation, handle collisions and streaming, apply boundary conditions, and visualize results.

## Features
- Flexible grid size and simulation parameters.
- Implements core LBM operations (collision, streaming, boundary conditions).
- Parallel impementation using OPEN_MP.

## File Structure
- **`LBM.h`**: Header file defining the `LBM` class, constants and its methods.
- **`LBM.cpp`**: Implementation of the `LBM` class, including the main computational methods.
- **`main.cpp`**: Entry point for the simulation. Sets up the `LBM` object, initializes parameters, runs the simulation, and outputs results.

## Requirements
- g++ compiler with C++17 support or higher.
- OPEN_MP dev package.

## Compilation and Execution
1. Compile the project using the included bash file:
   ```bash
   ./compile.sh
   ```

2. Run the executable:
   ```bash
   ./sim
   ```

## How to Use
1. Modify the simulation parameters in `main.cpp` to suit your needs (e.g., grid size, time steps).
2. Compile and run the simulation.
3. Run the bin_to_gif.py to create the visualization representation of the simulation.

## Example Output
The simulation outputs results such as velocity fields, density distributions, or other relevant physical quantities. 
Each of them is saved in a bin file for each timestep and can be converted to a GIF.

The final result is a GIF like this:

<p align="center">
  <img width="400" alt="" src="extra/velocity_evolution.gif"> <img width="400" alt="" src="extra/velocity_evolution_arrows.gif">
</p>

## Code Validation
In order to verify the correctness and the physical behaviour of the model, a comparison with results from <a href="https://www.sciencedirect.com/science/article/pii/0021999182900584">"U. Ghia, K. N. Ghia, C. T. Shin, *High-Re solutions for incompressible flow using Navier-Stokes equations and multigrid method*"</a>. Data are refered to the geometrical center of the cavity, with Re=100.

<p align="center">
  <img width="400" alt="" src="extra/u_plot.png"> <img width="400" alt="" src="extra/v_plot.png">
</p>

## Strong and Weak Scalability Test
The model has been parallelized using OpenMP. In order to understand how performance scales with the dimension of the domain and with the number of threads used during the computation, some scalability tests have been performed.

<p align="center">
  <img width="400" alt="" src="extra/strong_scalability.png"> <img width="400" alt="" src="extra/weak_scalability.png">
</p>

## Authors
Luca Donato, \
Leonardo Arnaboldi, \
Fabio Ghattas,\
 Luis Felipe Epia Realpe,\
 Tommaso Angelaccio

