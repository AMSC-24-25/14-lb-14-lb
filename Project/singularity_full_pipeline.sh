#!/bin/bash

CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp"

# for mac:
# g++ LBM.o seconds.o main.o -o sim

#classes at the moment
cmake .
make

singularity  build --fakeroot --force lbm.sif conf.def

singularity exec --bind ./results/bin_results:$HOME/results/bin_results lbm.sif mpirun -np 8 ./build/lbm_sim
