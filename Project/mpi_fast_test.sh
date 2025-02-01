#!/bin/bash

cmake .
make
./build/lbm_sim
python results/tools/bin_to_gif.py