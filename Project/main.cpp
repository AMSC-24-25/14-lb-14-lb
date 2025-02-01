#include <LatticeBoltzmannMethod/LBM.hpp>
#include <LatticeBoltzmannMethod/LidDrivenCavity3D.hpp>
#include <LatticeBoltzmannMethod/ObstacleLiftDrag.hpp>
#include <iostream>
#include <fstream>
#include <mpi.h>

const unsigned int scale = 2;
const unsigned int NX = 64*scale;
const unsigned int NY = NX;
const unsigned int NZ = NX;
const double nu = 1.0 / 6.0;
const unsigned int NSTEPS = 500*scale*scale + 1; 
const unsigned int NSAVE  =  10*scale*scale;
const unsigned int NMSG   =  50*scale*scale;

namespace lb = LatticeBoltzmannMethod;

int main(int argc, char* argv[])
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

    //lbm object initialization
    lb::LBM::dimensions d = {NX, NY, NZ};
    lb::LBM lbm = lb::LBM(world_rank, world_size, lb::LBM::VelocitySet::D3Q27, d, nu);

    if (world_rank == 0){
        std::cout << "\n------------------------------------------------------" << std::endl;        
    }
    lbm.setInitialCondition(lb::LidDrivenCavityInitial);
    if (world_rank == 0) {
        std::cout << "Initial condition set for Lid Driven Cavity simulation." << std::endl;
    }
    lbm.addBoundaryCondition(lb::BounceBackAllBut001);
    if (world_rank == 0) {
        std::cout << "Set any other wall for bounce-back." << std::endl;
    }
    lbm.addBoundaryCondition(lb::MovingWall001);
    if (world_rank == 0) {
        std::cout << "Set upper (considering z axis) boundary as moving wall." << std::endl;
    }
    auto ob = std::make_shared<lb::ObstacleLiftDrag>(3.0/7.0 * NX, 5.5/7.0 * NY, 3.0/7.0 * NZ , 1.0/7.0 * NX, 1.0/7.0 * NY, 1.0/7.0 * NZ);
    lbm.addObstacle(ob);
    if (world_rank == 0) {
        std::cout << "Set boundary for obstacle" << std::endl;
    }


    // Write simulation parameters to CSV file
    if (world_rank ==0){
        const double u_max = lb::Re/(6*NX);
        std::ofstream csvFile("results/simulation_parameters.csv");
        if (csvFile.is_open())
        {
            csvFile << "NX,NY,NZ,NSTEPS,NSAVE,U_MAX,RE,WORLD_SIZE\n";
            csvFile << NX << "," << NY << "," << NZ << "," << NSTEPS << "," << NSAVE << "," << u_max << "," << lb::Re << "," << world_size << "\n";
            csvFile.close();
            std::cout << "Simulation parameters saved to simulation_parameters.csv." << std::endl;
        }
        else
        {
            std::cerr << "Unable to open file simulation_parameters.csv" << std::endl;
        }
    }

    lbm.init_equilibrium();

    if (world_rank == 0){
        std::cout << "------------------------------------------------------\n" << std::endl;        
    }
    
    // main simulation loop; take NSTEPS time steps
    for(unsigned int n = 0; n < NSTEPS; ++n)
    {

        lbm.stream_collide_save();

        if(n % NSAVE == 0)
        {
            lbm.saveToBin(n);
        }

    }

 
    MPI_Finalize();
}

