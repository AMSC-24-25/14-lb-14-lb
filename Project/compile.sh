#!/bin/bash

rm -f ldc_sim

CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp -I /usr/include/ -fopenmp -O3"


# for mac:
# g++ LBM.o seconds.o main.o -o sim

#classes at the moment
g++ ${CXXFLAGS} LBM.hpp LBM.cpp VelocitySet.cpp LidDrivenCavity3D.hpp LidDrivenCavity3D.cpp  main.cpp -o ldc_sim
