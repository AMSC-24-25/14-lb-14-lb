#ifndef __OLD_H
#define __OLD_H
#include "LBM.hpp"

void BounceBackObstacle(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);
/*void ObstacleEast(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);
void ObstacleWest(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);
void ObstacleSouth(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);
void ObstacleNorth(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l);*/

#endif