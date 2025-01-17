#include "LBMobj.hpp"

void LidDrivenCavityInitial(int x, int y, int z, LBM& l);
//void BounceBackNorth(int x, int y, int z, Eigen::VectorXd& f, LBM& l);
void BounceBackEast(int x, int y, int z, Eigen::VectorXd& f, LBM& l);
void BounceBackWest(int x, int y, int z, Eigen::VectorXd& f, LBM& l);
void BounceBackSouth(int x, int y, int z, Eigen::VectorXd& f, LBM& l);
void MovingWallNorth(int x, int y, int z, Eigen::VectorXd& f, LBM& l);
