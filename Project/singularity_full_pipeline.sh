#!/bin/bash -l

cd build
cmake ..
make
cd ..
singularity  build --fakeroot --force lbm.sif conf.def

singularity exec --bind bin_results:$HOME/bin_results lbm.sif mpirun -np 8 ./ldc_sim_mpi
