#ifndef _LBM_OBSTALCE
#define _LBM_OBSTALCE
#include <Eigen/Dense>

namespace LatticeBoltzmannMethod{

    class LBM;
    //an abstract class to add customized obstacles inside an LBM domain
    //(derived classe's instantitation)
    class Obstacle
    {
        public:
            virtual void operator()(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f, LBM& l) const = 0;
            virtual ~Obstacle() = default;
    };
}
#endif