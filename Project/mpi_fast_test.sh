#!/bin/bash -l

cd build
cmake ..
make
cd ..
./build/lbm_sim
python results/tools/bin_to_gif.py