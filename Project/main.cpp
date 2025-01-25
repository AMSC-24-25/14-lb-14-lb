#include "seconds.h"
#include "LBM.hpp"
#include "LidDrivenCavity3D.hpp"
#include <iostream>

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

    //lbm object initialization
    LBM::dimensions d = {NX, NY, NZ};
    LBM lbm = LBM(LBM::VelocitySet::D3Q19, d, nu);
    lbm.setInitialCondition(LidDrivenCavityInitial);
    std::cout << "Initial condition set for Lid Driven Cavity simulation." << std::endl;
    lbm.addBoundaryCondition(MovingWall001);
    std::cout << "Set upper (consifdering z axis) boundary as moving wall." << std::endl; 
    lbm.addBoundaryCondition(BounceBackAllBut001);
    std::cout << "Set any other wall for bounce-back." << std::endl; 

    //double start = seconds();
    #pragma omp parallel master 
    {
        lbm.init_equilibrium();
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

