#!/bin/bash -l
module load mpi/openmpi-x86_64

rm -rf ./build/*
cmake . && make

singularity  build --fakeroot --force lbm.sif conf.def

singularity exec --bind ./results/bin_results:$HOME/results/bin_results lbm.sif mpirun -np 8 ./build/lbm_sim
