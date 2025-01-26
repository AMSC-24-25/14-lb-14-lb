#include "ObstacleLiftDrag.hpp"
#include <eigen3/Eigen/Dense>
#include "LBMobj.hpp"
using Eigen::VectorXd;





/*void init_obstacle(int x, int y, int z, int lenght, int height, int depth) {    
    obstacle.x = x;
    obstacle.y = y;
    obstacle.z = z;
    obstacle.length = lenght;
    obstacle.height = height;
    obstacle.depth = depth;
}*/


void ObstacleEast(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l)
{
    //int right = l.obstacle.x + l.obstacle.length;
    unsigned int top = l.obstacle.y + l.obstacle.height;
    //int dep = l.obstacle.z + l.obstacle.depth;
    //if(x >= obstacle.x && x < right && y >= obstacle.y && y < top && z >= obstacle.z && z < dep)
    if(x == l.obstacle.x && y >= l.obstacle.y && y < top /*&& z >= obstacle.z && z < dep*/)
    {
        VectorXd p = l.getPopulation(x,y,z);
        f(3) = p(1);
        f(7) = p(8);
        f(6) = p(5);
        //f(1) = 0;
        //f(5) = 0;
        //f(8) = 0;
    }
}

void ObstacleWest(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l)
{
    unsigned int right = l.obstacle.x + l.obstacle.length;
    unsigned int top = l.obstacle.y + l.obstacle.height;
    //int dep = l.obstacle.z + l.obstacle.depth;
    if(x == right && y >= l.obstacle.y && y < top /*&& z >= obstacle.z && z < dep*/)
    {
        VectorXd p = l.getPopulation(x,y,z);
        f(1) = p(3);
        f(5) = p(6);
        f(8) = p(7);
        //f(3) = 0;
        //f(7) = 0;
        //f(6) = 0;
    }
}

void ObstacleSouth(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l)
{
    unsigned int right = l.obstacle.x + l.obstacle.length;
    unsigned int top = l.obstacle.y + l.obstacle.height;
    //int dep = l.obstacle.z + l.obstacle.depth;
    if(x >= l.obstacle.x && x < right && y == top /*&& z >= obstacle.z && z < dep*/)
    {
        VectorXd p = l.getPopulation(x,y,z);
        f(2) = p(4);
        f(5) = p(8);
        f(6) = p(7);
        //f(4) = 0;
        //f(8) = 0;
        //f(7) = 0;
    }
}

void ObstacleNorth(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l)
{
    unsigned int right = l.obstacle.x + l.obstacle.length;
    //int top = l.obstacle.y + l.obstacle.height;
    //int dep = l.obstacle.z + l.obstacle.depth;
    if(x >= l.obstacle.x && x < right && y == l.obstacle.y /*&& z >= obstacle.z && z < dep*/)
    {
        VectorXd p = l.getPopulation(x,y,z);
        f(4) = p(2);
        f(8) = p(5);
        f(7) = p(6);
        //f(2) = 0;
        //f(5) = 0;
        //f(6) = 0;
    }
}