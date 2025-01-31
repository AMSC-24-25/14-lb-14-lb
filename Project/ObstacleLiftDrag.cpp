#include "ObstacleLiftDrag.hpp"

using Eigen::VectorXd;

namespace LatticeBoltzmannMethod{

    ObstacleLiftDrag::ObstacleLiftDrag(unsigned int x_start, unsigned int y_start, unsigned int z_start, unsigned int lenght, unsigned int height, unsigned int depth)
    {
        this->x0 = x_start;
        this->y0 = y_start;
        this->z0 = z_start;
        this->length = lenght;
        this->height = height;
        this->depth = depth;
    }

    inline bool ObstacleLiftDrag::is_inside_point(unsigned int x, unsigned int y, unsigned int z) const
    {
        bool x_inside = (x >= this->x0) && (x <= (this-> x0 + this->length));
        bool y_inside = (y >= this->y0) && (y <= (this-> y0 + this->height));
        bool z_inside = (z >= this->z0) && (z <= (this-> z0 + this->depth));
        return x_inside && y_inside && z_inside;
    }


    void ObstacleLiftDrag::operator()(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l) const
    {
        bool inside = is_inside_point(x, y, z);
        if(!inside) return;

        

        unsigned int obx = this->x0;
        unsigned int oby = this->y0;
        unsigned int obz = this->z0;
        unsigned int right = this->x0 + this->length;
        unsigned int top = this->y0 + this->height;
        unsigned int dep = this->z0 + this->depth;
        double x_bounce = 0.0, y_bounce = 0.0, z_bounce = 0.0;
        if(x == right) x_bounce = -1.0;
        else if(x == obx) x_bounce = 1.0;
        if(y == top) y_bounce = -1.0;
        else if(y == oby) y_bounce = 1.0;
        if(z == dep) z_bounce = -1.0;
        else if (z == obz) z_bounce = 1.0;

        

        if(x_bounce == 0.0 && y_bounce == 0.0 && z_bounce == 0.0)
        {
            f = l.get_rho(x, y, z) * l.v->get_w();
            return;
        }

        const Eigen::MatrixXd& c = l.v->get_c();
        VectorXd p = l.getPopulation(x,y,z);
        for(unsigned int i = 0; i < l.v->getQ(); i++)
        {
            if(c(i, 0) == x_bounce || c(i, 1) == y_bounce || c(i, 2) == z_bounce)
            {
                unsigned int j = l.v->directionIndexInvert(i);
                f(j) = p(i); 
                //std::cout << "yessir" << std::endl;
            }
        }
        //std::cout << "(" << x <<  "," << y << ", " << z <<  ")         ";
    }
}
/*
void ObstacleEast(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l)
{
    //int right = l.obstacle.x + l.obstacle.length;
    unsigned int top = l.obstacle.y + l.obstacle.height;
    int dep = l.obstacle.z + l.obstacle.depth;
    //if(x >= obstacle.x && x < right && y >= obstacle.y && y < top && z >= obstacle.z && z < dep)
    if(x == l.obstacle.x && y >= l.obstacle.y && y < top && z >= obstacle.z && z < dep)
    {
        VectorXd p = l.getPopulation(x,y,z);
        f(3) = p(l.directionIndexInvert(3));
        f(7) = p(l.directionIndexInvert(7));
        f(6) = p(l.directionIndexInvert(6));
        //f(1) = 0;
        //f(5) = 0;
        //f(8) = 0;
    }
}

void ObstacleWest(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l)
{
    unsigned int right = l.obstacle.x + l.obstacle.length;
    unsigned int top = l.obstacle.y + l.obstacle.height;
    int dep = l.obstacle.z + l.obstacle.depth;
    if(x == right && y >= l.obstacle.y && y < top && z >= obstacle.z && z < dep)
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
    int dep = l.obstacle.z + l.obstacle.depth;
    if(x >= l.obstacle.x && x < right && y == top && z >= obstacle.z && z < dep)
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
    int dep = l.obstacle.z + l.obstacle.depth;
    if(x >= l.obstacle.x && x < right && y == l.obstacle.y && z >= obstacle.z && z < dep)
    {
        VectorXd p = l.getPopulation(x,y,z);
        f(4) = p(2);
        f(8) = p(5);
        f(7) = p(6);
        //f(2) = 0;
        //f(5) = 0;
        //f(6) = 0;
    }
}*/