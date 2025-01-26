#include "LidDrivenCavity.hpp"

using Eigen::VectorXd;

const double Re = 100.0;
const double rho0 = 1.0;

void LidDrivenCavityInitial(unsigned int x, unsigned int y, unsigned int z, LBM &l)
{
    l.set_rho(x, y, z, rho0);

    VectorXd u(2);
    if (y == l.N.y - 1)
    {
        double u_max = Re / (6 * l.N.x);
        u << u_max, 0.0;
    }
    else
        u << 0.0, 0.0;

    l.set_u(x, y, z, u);
}

// CHANGED ALL BOUNCEBACK FUNCTIONS TO REFLECT ONLY THE NORMAL COMPONENT TO THE BOUNDARY
void BounceBackEast(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd &f, LBM &l)
{
    if (x == (l.N.x - 1))
    {
        VectorXd p = l.getPopulation(x, y, z);
        f(3) = p(1);
        f(7) = p(8);
        f(6) = p(5);
    }
}

void BounceBackWest(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd &f, LBM &l)
{
    if (x == 0)
    {
        VectorXd p = l.getPopulation(x, y, z);
        f(1) = p(3);
        f(5) = p(6);
        f(8) = p(7);
    }
}

void BounceBackSouth(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd &f, LBM &l)
{
    if (y == 0)
    {
        VectorXd p = l.getPopulation(x, y, z);
        f(2) = p(4);
        f(5) = p(8);
        f(6) = p(7);
    }
}

void MovingWallNorth(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd &f, LBM &l)
{
    if (y == (l.N.y - 1))
    {
        double u_lid = Re / (6 * l.N.x);
        double rho_top = l.get_rho(x, y, z);
        VectorXd p = l.getPopulation(x, y, z);

        double coeff = 2.0 * (1.0 / 36.0) * (1.0 / l.cs) * (1.0 / l.cs) * rho_top;

        f(4) = p(2);
        f(7) = p(5) - coeff * u_lid;
        f(8) = p(6) + coeff * u_lid;
    }
}