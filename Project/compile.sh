#!/bin/bash

rm -f LBM.o seconds.o main.o sim

CXXFLAGS="-std=c++17"

nvcc ${CXXFLAGS} -c LBM.cu -o LBM.o
nvcc ${CXXFLAGS} -c seconds.cpp -o seconds.o
nvcc ${CXXFLAGS} -c main.cu -o main.o
 
nvcc LBM.o seconds.o main.o -o sim -lrt

rm -f LBM.o seconds.o main.o
./sim
# for mac:
# g++ LBM.o seconds.o main.o -o sim

