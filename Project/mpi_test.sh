CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp -I ${mkEigenInc} -fopenmp -O3"

# for mac:
# g++ LBM.o seconds.o main.o -o sim

#classes at the moment
/usr/lib64/openmpi/bin/mpic++ ${CXXFLAGS} -I/usr/include LBMobj.hpp LBMobj.cpp VelocitySet.cpp LidDrivenCavity.hpp LidDrivenCavity.cpp  main.cpp -o ldc_sim_mpi && \
/usr/lib64/openmpi/bin/mpirun -n 2 ./ldc_sim_mpi