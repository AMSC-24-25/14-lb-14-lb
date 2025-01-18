#!/bin/bash

rm -f LBM.o seconds.o main.o sim

CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp"

g++ ${CXXFLAGS} -c LBM.cpp -o LBM.o
g++ ${CXXFLAGS} -c seconds.cpp -o seconds.o
g++ ${CXXFLAGS} -c main.cpp -o main.o
 
g++ -fopenmp LBM.o seconds.o main.o -o sim -lrt

rm -f LBM.o seconds.o main.o

# for mac:
# g++ LBM.o seconds.o main.o -o sim

#classes at the moment
#g++ -Wall -I ${mkEigenInc} -fopenmp LBMobj.hpp LBMobj.cpp VelocitySet.cpp LidDrivenCavity.hpp LidDrivenCavity.cpp  main.cpp -o ldc_sim
