#ifndef __LDC_H
#define __LDC_H
#include <LatticeBoltzmannMethod/LBM.hpp>

namespace LatticeBoltzmannMethod{

const double Re = 100.0;
const double rho0 = 1.0;

void LidDrivenCavityInitial(unsigned int x, unsigned int y, unsigned int z, LBM& l);
void MovingWall001(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);
void BounceBackAllBut001(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);
}
#endif