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

## Example Output
The simulation outputs results such as velocity fields, density distributions, or other relevant physical quantities. 
Check simulation_data.csv as example.


## Authors
Luca Donato, \
Leonardo Arnaboldi, \
Fabio Ghattas,\
 Luis Felipe Epia Realpe,\
 Tommaso Angelaccio

