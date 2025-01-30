#ifndef __OLD_H
#define __OLD_H
#include "Obstacle.hpp"
#include "LBM.hpp"

namespace LatticeBoltzmannMethod{

    class ObstacleLiftDrag : public Obstacle
    {
        public:
            unsigned int x;
            unsigned int y;
            unsigned int z;
            unsigned int length;
            unsigned int height;
            unsigned int depth;

            ObstacleLiftDrag(unsigned int x0, unsigned int y0, unsigned int0, unsigned int lenght, unsigned int height, unsigned int depth);

            void operator()(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l) const override;

            bool is_inside_point(unsigned int x, unsigned int y, unsigned int z) const;
    };
}
#endif