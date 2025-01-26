

#ifndef __LDC_H
#define __LDC_H
#include "LBMobj.hpp"
void LidDrivenCavityInitial(unsigned int x, unsigned int y, unsigned int z, LBM &l);
// void BounceBackNorth(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);
void BounceBackEast(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd &f, LBM &l);
void BounceBackWest(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd &f, LBM &l);
void BounceBackSouth(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd &f, LBM &l);
void MovingWallNorth(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd &f, LBM &l);

#endif