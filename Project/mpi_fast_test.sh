#!/bin/bash -l
module load mpi/openmpi-x86_64

# for mac:
# g++ LBM.o seconds.o main.o -o sim

rm ldc_sim_mpi
cmake . && make
mpirun -n 8 ./ldc_sim_mpi
python bin_to_gif.py