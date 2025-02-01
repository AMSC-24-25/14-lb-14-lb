#!/bin/bash -l
module load mpi/openmpi-x86_64

rm -rf ./build/*
cmake . && make
mpirun -n 8 ./ldc_sim_mpi
python bin_to_gif.py