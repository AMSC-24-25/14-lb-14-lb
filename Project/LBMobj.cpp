#include "LBMobj.hpp"

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
                double rho = get_rho(x,y,z);
                VectorXd u = get_u(x,y,z);
                double unorm = u.norm();
                double omusq = 1.0 - 1.5*(unorm * unorm);
                ArrayXd c3u = this->v.get_c() * (u * 3.0);
                VectorXd f = this->v.get_w().array().cwiseproduct(1.0 * omusq + c3u.cwiseproduct(c3u*0.5 + 1.0));

                savePopulation(x,y,z, &f);
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

                VectorXd cf_rhoinv = (this->v.get_c().transpose() * stream) * rhoinv;
                VectorXd& u = cf_rhoinv;
                set_u(x,y,z, u );

                //collision
                double unorm = u.norm();
                double omusq = 1.0 - 1.5*(unorm * unorm);
                u *= 3.0;
                ArrayXd c3u = this->v.get_c() * (u * 3.0);
                VectorXd twr = this->v.get_w() * tauinv * rho;

                stream *= omtauinv;
                stream += twr.array().cwiseproduct( omusq + c3u.cwiseproduct(1.0 + 0.5*c3u) ); 

                this->savePopulation(x,y,z, stream);
            }
        }
    }
}


