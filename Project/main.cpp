#include "seconds.h"
#include "LBMobj.hpp"
#include "LidDrivenCavity.hpp"
#include "ObstacleLiftDrag.hpp"
#include <iostream>
#include <mpi.h>

const int scale = 2;
const int NX = 64 * scale;
const int NY = NX;
const double nu = 1.0 / 6.0;
const int NSTEPS = 7; // Added + 1 just to test the code
const int NSAVE = 50;
const int NMSG = 50 * scale * scale;

int world_size;
int world_rank;

int main(int argc, char *argv[])
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    printf("Hello world from processor %s, rank %d out of %d processors\n", processor_name, world_rank, world_size);

    // lbm object initialization
    LBM::dimensions d = {NX, NY, 1};
    LBM::partition_config dims = {world_rank, world_size};
    LBM lbm = LBM(LBM::VelocitySet::D2Q9, dims, d, nu);
    lbm.setInitialCondition(LidDrivenCavityInitial);
    std::cout << "Initial condition set for Lid Driven Cavity simulation." << std::endl;
    lbm.addBoundaryCondition(BounceBackEast);
    std::cout << "Set East Boundary for bounce back." << std::endl;
    lbm.addBoundaryCondition(BounceBackWest);
    std::cout << "Set West Boundary for bounce back." << std::endl;
    lbm.addBoundaryCondition(BounceBackSouth);
    std::cout << "Set South Boundary for bounce back." << std::endl;
    lbm.addBoundaryCondition(MovingWallNorth);
    std::cout << "Set North Boundary as moving wall." << std::endl;
    lbm.addBoundaryCondition(ObstacleEast);
    std::cout << "Set East Boundary for obstacle." << std::endl;
    lbm.addBoundaryCondition(ObstacleWest);
    std::cout << "Set West Boundary for obstacle." << std::endl;
    lbm.addBoundaryCondition(ObstacleSouth);
    std::cout << "Set South Boundary for obstacle." << std::endl;
    lbm.addBoundaryCondition(ObstacleNorth);
    std::cout << "Set North Boundary for obstacle." << std::endl;

    lbm.init_equilibrium();
    lbm.init_obstacle();

    lbm.init_equilibrium();
    // main simulation loop; take NSTEPS time steps
    for (int n = 0; n < NSTEPS; ++n)
    {
        lbm.stream_collide_save();
        if (n % NSAVE == 0)
        {
            lbm.saveToBin(n, world_rank);
        }
    }

    MPI_Finalize();
}
