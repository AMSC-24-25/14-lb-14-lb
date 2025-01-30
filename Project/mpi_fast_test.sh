#!/bin/bash
module load mpi/openmpi-x86_64

CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp"

# for mac:
# g++ LBM.o seconds.o main.o -o sim

rm ldc_sim_mpi
mpic++ ${CXXFLAGS} -I/usr/include LBM.hpp LBM.cpp VelocitySet.cpp LidDrivenCavity3D.hpp LidDrivenCavity3D.cpp ObstacleLiftDrag.hpp ObstacleLiftDrag.cpp  main.cpp -o ldc_sim_mpi 
mpirun -n 14 ./ldc_sim_mpi