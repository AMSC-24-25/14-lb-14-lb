#!/bin/bash -l

cd build
cmake ..
make
cd ..
./build/lbm_sim
python results/tools/bin_to_vtk.py
#python results/tools/bin_to_validation_csv.py
python results/tools/bin_to_gif.py