#!/bin/bash

rm -f ldc_sim
rm -f ldc_sim

CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp -I /usr/include -fopenmp -O3"

# for mac:
# g++ LBM.o seconds.o main.o -o sim

#classes at the moment
g++ ${CXXFLAGS} LBMobj.hpp LBMobj.cpp VelocitySet.cpp LidDrivenCavity.hpp LidDrivenCavity.cpp  main.cpp -o ldc_sim