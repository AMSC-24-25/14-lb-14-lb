#include <LatticeBoltzmannMethod/LidDrivenCavity3D.hpp>

using Eigen::VectorXd;

namespace LatticeBoltzmannMethod{
    void LidDrivenCavityInitial(unsigned int x, unsigned int y, unsigned int z, LBM& l)
    {
        l.set_rho(x,y,z, rho0);
        
        VectorXd u(3);
        if(z == l.N.z - 1)
        {
            double u_max = Re/(6*l.N.x);
            u << u_max, 0.0, 0.0;
        }
        else u << 0.0, 0.0, 0.0;

        l.set_u(x,y,z, u);
    }



    void MovingWall001(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l)
    {
        if(z == (l.N.z - 1))
        {
            double u_lid = Re/(6*l.N.x);
            double rho_top = l.get_rho(x,y,z);
            

            

            const Eigen::MatrixXd& c = l.v->get_c();
            VectorXd p = l.getPopulation(x,y,z);
            for(unsigned int i = 0; i < l.v->getQ(); i++)
            {
                if(c(i, 2) == 1.0)
                {
                    //page 180 of the book
                    double coeff = 2.0 * l.v->get_w()(i) * (1.0 / l.cs) * (1.0 / l.cs) * rho_top;
                    unsigned int j = l.v->directionIndexInvert(i);
                    f(j) = p(i) - coeff * (c(i, 0) * u_lid);
                }
            }
        }
    }

    void BounceBackAllBut001(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l)
    {
        double x_bounce = 0.0, y_bounce = 0.0, z_bounce = 0.0;
        if(x == 0) x_bounce = -1.0;
        else if(x == l.N.x - 1) x_bounce = 1.0;
        if(y == 0) y_bounce = -1.0;
        else if(y == l.N.y - 1) y_bounce = 1.0;
        if(z == 0) z_bounce = -1.0;

        if(x_bounce == 0.0 && y_bounce == 0.0 && z_bounce == 0.0) return;

        const Eigen::MatrixXd& c = l.v->get_c();
        VectorXd p = l.getPopulation(x,y,z);
        for(unsigned int i = 0; i < l.v->getQ(); i++)
        {
            if(c(i, 0) == x_bounce || c(i, 1) == y_bounce || c(i, 2) == z_bounce)
            {
                unsigned int j = l.v->directionIndexInvert(i);
                f(j) = p(i); 
            }
        }
    }
}