#include <LatticeBoltzmannMethod/ObstacleLiftDrag.hpp>

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
            }
        }

    }
}