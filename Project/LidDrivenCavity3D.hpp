#ifndef __LDC_H
#define __LDC_H
#include "LBM.hpp"
void LidDrivenCavityInitial(unsigned int x, unsigned int y, unsigned int z, LBM& l);
void MovingWall001(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);
void BounceBackAllBut001(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);
#endif