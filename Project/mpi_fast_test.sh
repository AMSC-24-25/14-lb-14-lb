#!/bin/bash -l
module load mpi/openmpi-x86_64

rm -rf ./build/*
cmake . && make
mpirun -n 8 ./build/lbm_sim
python ./results/tools/bin_to_gif.py