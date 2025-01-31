#!/bin/bash -l
module load mpi/openmpi-x86_64

CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp -I  /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3 -I /usr/include/ -I /usr/include/x86_64-linux-gnu -fopenmp -O3"

# for mac:
# g++ LBM.o seconds.o main.o -o sim

rm ldc_sim_mpi
mpic++ ${CXXFLAGS} -I/usr/include LBM.hpp LBM.cpp VelocitySet.cpp LidDrivenCavity3D.hpp LidDrivenCavity3D.cpp ObstacleLiftDrag.hpp ObstacleLiftDrag.cpp  main.cpp -o ldc_sim_mpi 
mpirun -n 8 ./ldc_sim_mpi
python bin_to_gif.py