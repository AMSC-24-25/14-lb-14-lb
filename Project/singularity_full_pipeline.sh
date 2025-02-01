#!/bin/bash -l

module load mpi/openmpi-x86_64

CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp"

# for mac:
# g++ LBM.o seconds.o main.o -o sim

#classes at the moment
cmake . && make

singularity  build --fakeroot --force lbm.sif conf.def

singularity exec --bind bin_results:$HOME/bin_results lbm.sif mpirun -np 8 ./lbm_sim
