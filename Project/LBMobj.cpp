#include "LBMobj.hpp"
#include <iostream>

using Eigen::VectorXd;
using Eigen::ArrayXd;

void LBM::init_equilibrium()
{

    #pragma omp parallel for default(none) shared(this) schedule(static)
    for(unsigned int x = 0; x < this->N.x; ++x)
    {
        for(unsigned int y = 0; y < this->N.y; ++y)
        {
            for(unsigned int z = 0; z < this->N.z; ++z)
            {
                applyInitial(x,y,z);

                double rho = get_rho(x,y,z);
                VectorXd u = get_u(x,y,z);
                double unorm = u.norm();
                double omusq = 1.0 - 1.5*(unorm * unorm);
                ArrayXd c3u = this->v->get_c() * (u * 3.0);
                VectorXd f = this->v->get_w().array() * (1.0 * omusq + c3u * (c3u*0.5 + 1.0));

                savePopulation(x,y,z, f);
            }
        }
    }
}

void LBM::stream_collide_save()
{
        // useful constants
    const double tauinv = 2.0/(6.0*this->nu+1.0); // 1/tau
    const double omtauinv = 1.0-tauinv;     // 1 - 1/tau

    #pragma omp parallel for default(none) shared(this, tauinv, omtauinv) schedule(static)
    for(unsigned int x = 0; x < this->N.x; ++x)
    {
        for(unsigned int y = 0; y < this->N.y; ++y)
        {
            for(unsigned int z = 0; z < this->N.z; ++z)
            {
                //classical stream
                VectorXd stream = populationAdjacent(x,y,z);
                //modify stream applying boundaries
                applyBoundary(x,y,z, stream);

                //compute moments
                double rho = stream.sum();
                set_rho(x,y,z, rho );

                double rhoinv = 1.0/rho;

                VectorXd cf_rhoinv = (this->v->get_c().transpose() * stream) * rhoinv;
                VectorXd& u = cf_rhoinv;
                set_u(x,y,z, u );

                //collision
                double unorm = u.norm();
                double omusq = 1.0 - 1.5*(unorm * unorm);
                u *= 3.0;
                ArrayXd c3u = this->v->get_c() * (u * 3.0);
                VectorXd twr = this->v->get_w() * tauinv * rho;

                stream *= omtauinv;
                stream += (twr.array() * ( omusq + c3u * (1.0 + 0.5*c3u) )).matrix(); 

                this->savePopulation(x,y,z, stream);
            }
        }
    }
}


LBM::LBM(LBM::VelocitySet::StandardSet vSet, LBM::dimensions d,  double nu) : N(d), nu(nu)
{
    v = std::make_unique<VelocitySet>(vSet);
    if(N.x == 0 || N.y == 0 || N.z == 0)
    {
        std::cerr << "Setting problem size to zero!!\n";
        exit(1);
    }
    if((v->getD() == 1 && (N.y != 1 || N.z != 1)) || (v->getD() == 2 && N.z != 1))
    {
        std::cerr << "Incompatible size and velocity set dimension!! The problem is "  << v->getD() << "D\n";
        exit(1);
    }    
    population = (double*) malloc(sizeof(double) * N.x * N.y * N.z * v->getQ() * 2);
    rho = (double*) malloc(sizeof(double) * N.x * N.y * N.z);
    u = (double*) malloc(sizeof(double) * N.x * N.y * N.z * v->getD());
}
