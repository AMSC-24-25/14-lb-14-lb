#include "seconds.h"
#include "LBM.hpp"
#include "LidDrivenCavity3D.hpp"
#include "ObstacleLiftDrag.hpp"
#include <iostream>
#include <omp.h>
#include <fstream>

const unsigned int scale = 2;
const unsigned int NX = 64*scale;
const unsigned int NY = NX;
const unsigned int NZ = NX;
const double nu = 1.0 / 6.0;
const unsigned int NSTEPS = 50*scale*scale + 1; //Added + 1 just to test the code
const unsigned int NSAVE  =  10*scale*scale;
const unsigned int NMSG   =  50*scale*scale;

int main(int argc, char* argv[])
{
    omp_set_num_threads(omp_get_max_threads());

    //lbm object initialization
    LBM::dimensions d = {NX, NY, NZ};
    LBM lbm = LBM(LBM::VelocitySet::D3Q27, d, nu);
    lbm.setInitialCondition(LidDrivenCavityInitial);
    std::cout << "Initial condition set for Lid Driven Cavity simulation." << std::endl;
    lbm.addBoundaryCondition(BounceBackAllBut001);
    std::cout << "Set any other wall for bounce-back." << std::endl; 
    lbm.addBoundaryCondition(MovingWall001);
    std::cout << "Set upper (consifdering z axis) boundary as moving wall." << std::endl; 
    lbm.addBoundaryCondition(BounceBackObstacle);
    std::cout << "Set boundary for obstacle" << std::endl; 

    // Write simulation parameters to CSV file
    const double u_max = Re/(6*NX);
    std::ofstream csvFile("simulation_parameters.csv");
    if (csvFile.is_open())
    {
        csvFile << "NX,NY,NZ,NSTEPS,NSAVE,U_MAX,RE\n";
        csvFile << NX << "," << NY << "," << NZ << "," << NSTEPS << "," << NSAVE << "," << u_max << "," << Re << "\n";
        csvFile.close();
        std::cout << "Simulation parameters saved to simulation_parameters.csv." << std::endl;
    }
    else
    {
        std::cerr << "Unable to open file simulation_parameters.csv" << std::endl;
    }

    //double start = seconds();

    lbm.init_equilibrium();
    lbm.init_obstacle();
    // main simulation loop; take NSTEPS time steps
    for(unsigned int n = 0; n < NSTEPS; ++n)
    {

        lbm.stream_collide_save();

        if(n % NSAVE == 0)
        {
            lbm.saveToBin(n);
        }

    }

    
    //double end = seconds();
    //double runtime = end-start;

 
    return 0;
}

