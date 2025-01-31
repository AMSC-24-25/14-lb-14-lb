#ifndef _LBM_OBSTALCE
#define _LBM_OBSTALCE
#include <eigen3/Eigen/Dense>

namespace LatticeBoltzmannMethod{

    class LBM;

    class Obstacle
    {
        public:
            virtual void operator()(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l) const = 0;
            virtual ~Obstacle() = default;
    };
}
#endif