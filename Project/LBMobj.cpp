#include "LBMobj.hpp"


void init_equilibrium()
{

    #pragma omp parallel for default(none) shared(this) schedule(static)
    for(unsigned int x = 0; x < this.Nx(); ++x)
    {
        for(unsigned int y = 0; y < this.Ny(); ++y)
        {
            for(unsigned int z = 0; z < this.Nz(), ++z)
            {
                double rho = this.r(x,y,z);
                auto u  = this.u(x,y,z);
                
                double w0 = this.v.w(0);
                                 
                this.setRest(x,y,z, w0 * rho * omusq);

                for(unsigned int i = 0; i < this.v.getQ(); ++i)
                {
                    auto ci = this.v.c(i);
                    double wi = this.v.w(i);
                    double cidot3u = 0;
                    for(unsigned int k = 0; k < this.v.getD(); ++k) cidot3u = ci.get(k) * 3.0 * u.get(k);

                    this.setPopulation(x,y,z,i, wi * (omusq + cidot3u*(1.0+0.5*cidot3u)) );
                }
            }
        }
    }
}