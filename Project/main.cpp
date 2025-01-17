#include <stdio.h>
#include <stdlib.h>

#include "seconds.h"
#include "LBMobj.hpp"
#include "LidDrivenCavity.hpp"
#include <utility>

const unsigned int scale = 2;
const unsigned int NX = 128*scale;
const unsigned int NY = NX;
const double nu = 1.0 / 6.0;
const unsigned int NSTEPS = 3000*scale*scale;
const unsigned int NSAVE  =  100*scale*scale;
const unsigned int NMSG   =  50*scale*scale;

int main(int argc, char* argv[])
{

    //lbm object initialization
    LBM::dimensions d = {NX, NY, 1};
    LBM lbm = LBM(LBM::VelocitySet::D2Q9, d, nu);
    lbm.setInitialCondition(LidDrivenCavityInitial);
    lbm.addBoundaryCondition(BounceBackEast);
    lbm.addBoundaryCondition(BounceBackWest);
    lbm.addBoundaryCondition(BounceBackSouth);
    lbm.addBoundaryCondition(MovingWallNorth);

    double start = seconds();
    #pragma omp parallel master 
    {
// main simulation loop; take NSTEPS time steps
    for(unsigned int n = 0; n < NSTEPS; ++n)
    {

        lbm.stream_collide_save();

    }
    }
    
    double end = seconds();
    double runtime = end-start;

 
    return 0;
}

