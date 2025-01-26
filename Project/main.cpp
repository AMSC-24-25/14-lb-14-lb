#include "seconds.h"
#include "LBMobj.hpp"
#include "LidDrivenCavity.hpp"
#include "ObstacleLiftDrag.hpp"
#include <iostream>

const unsigned int scale = 2;
const unsigned int NX = 64*scale;
const unsigned int NY = NX;
const double nu = 1.0 / 6.0;
const unsigned int NSTEPS = 50*scale*scale + 1; //Added + 1 just to test the code
const unsigned int NSAVE  =  10*scale*scale;
const unsigned int NMSG   =  50*scale*scale;

int main(int argc, char* argv[])
{

    //lbm object initialization
    LBM::dimensions d = {NX, NY, 1};
    LBM lbm = LBM(LBM::VelocitySet::D2Q9, d, nu);
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

    //double start = seconds();
    #pragma omp parallel master 
    {
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
    }
    
    //double end = seconds();
    //double runtime = end-start;

 
    return 0;
}

