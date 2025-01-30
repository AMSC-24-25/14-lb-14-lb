#!/bin/bash

rm -f ldc_sim


CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp -I  /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3 -I /usr/include/ -I /usr/include/x86_64-linux-gnu"


# for mac:
# g++ LBM.o seconds.o main.o -o sim

#classes at the moment
g++ ${CXXFLAGS} LBM.hpp LBM.cpp VelocitySet.cpp LidDrivenCavity3D.hpp LidDrivenCavity3D.cpp ObstacleLiftDrag.hpp ObstacleLiftDrag.cpp  main.cpp -o ldc_sim
