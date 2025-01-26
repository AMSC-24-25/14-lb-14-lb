module load mpi/openmpi-x86_64

CXXFLAGS="-std=c++17 -pedantic -O3 -Wall -fopenmp -I ${mkEigenInc} -fopenmp -O3"

# for mac:
# g++ LBM.o seconds.o main.o -o sim

#classes at the moment
mpic++ ${CXXFLAGS} -I/usr/include LBMobj.hpp LBMobj.cpp VelocitySet.cpp LidDrivenCavity.hpp LidDrivenCavity.cpp ObstacleLiftDrag.hpp ObstacleLiftDrag.cpp main.cpp -o ldc_sim_mpi &&\
mpirun -n 2 ./ldc_sim_mpi