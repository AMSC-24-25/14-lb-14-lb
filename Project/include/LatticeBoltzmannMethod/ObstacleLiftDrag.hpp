#ifndef __OLD_H
#define __OLD_H
#include <Obstacle.hpp>
#include <LBM.hpp>

namespace LatticeBoltzmannMethod{

    class ObstacleLiftDrag : public Obstacle  //a parallelepipedal obstacle
    {
        public:
            unsigned int x0;   //starting vertex of the obstacle is (x0, y0, z0)
            unsigned int y0;
            unsigned int z0;
            unsigned int length;
            unsigned int height;
            unsigned int depth;

            ObstacleLiftDrag(unsigned int x_start, unsigned int y_start, unsigned int z_start, unsigned int lenght, unsigned int height, unsigned int depth);

            void operator()(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l) const override;

            inline bool is_inside_point(unsigned int x, unsigned int y, unsigned int z) const;
    };
}
#endif