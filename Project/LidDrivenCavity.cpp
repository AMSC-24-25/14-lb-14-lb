#include "LBMobj.hpp"

using Eigen::VectorXd;

const double Re = 100.0;
const double rho0 = 1.0;

void LidDrivenCavityInitial(int x, int y, int z, LBM& l);
//void BounceBackNorth(int x, int y, int z, Eigen::VectorXd& f, LBM& l);
void BounceBackEast(int x, int y, int z, Eigen::VectorXd& f, LBM& l);
void BounceBackWest(int x, int y, int z, Eigen::VectorXd& f, LBM& l);
void BounceBackSouth(int x, int y, int z, Eigen::VectorXd& f, LBM& l);
void MovingWallNorth(int x, int y, int z, Eigen::VectorXd& f, LBM& l);

void LidDrivenCavityInitial(int x, int y, int z, LBM& l)
{
    l.set_rho(x,y,z, rho0);
    
    VectorXd u(2);
    if(y == l.N.y - 1)
    {
        double u_max = Re/(6*l.N.x);
        u << u_max, 0.0;
    }
    else u << 0.0, 0.0;

    l.set_u(x,y,z, u);
}

void BounceBackEast(int x, int y, int z, Eigen::VectorXd& f, LBM& l)
{
    if(x == (l.N.x - 1))
    {
        VectorXd p = l.getPopulation(x,y,z);
        f(3) = p(1);
        f(7) = p(5);
        f(6) = p(8);
    }
}

void BounceBackWest(int x, int y, int z, Eigen::VectorXd& f, LBM& l)
{
    if(x == 0)
    {
        VectorXd p = l.getPopulation(x,y,z);
        f(1) = p(3);
        f(5) = p(7);
        f(8) = p(6);
    }
}

void BounceBackSouth(int x, int y, int z, Eigen::VectorXd& f, LBM& l)
{
    if(y == 0)
    {
        VectorXd p = l.getPopulation(x,y,z);
        f(2) = p(4);
        f(5) = p(7);
        f(6) = p(8);
    }
}

void MovingWallNorth(int x, int y, int z, Eigen::VectorXd& f, LBM& l)
{
    if(y == (l.N.y - 1))
    {
        double u_lid = Re/(6*l.N.x);
        double rho_top = l.get_rho(x,y,z);
        VectorXd p = l.getPopulation(x,y,z);

        double coeff = 2.0 * (1.0/36.0) * (1.0 / l.cs) * (1.0 / l.cs) * rho_top;

        f(4) = p(2);
        f(7) = p(5) - coeff * u_lid;
        f(8) = p(6) + coeff * u_lid;
    }
}