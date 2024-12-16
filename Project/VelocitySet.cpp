#include "LBMobj.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

LBM::VelocitySet::VelocitySet(StandardSet set) : D(fromStdD(set)), Q(fromStdQ(set))
{
    c = MatrixXd(D, Q);
    w = VectorXd(Q);

    if( Q == 9)
    {
        c << 0.0, 0.0,
             1.0, 0.0,
             0.0, 1.0,
            -1.0, 0.0,
            0.0, -1.0,
             1.0, 1.0,
            -1.0, 1.0,
           -1.0, -1.0,
            1.0, -1.0;

        constexpr double w0 = 4.0/9.0;  // zero weight
        constexpr double ws = 1.0/9.0;  // adjacent weight
        constexpr double wd = 1.0/36.0;
        w << w0, ws, ws, ws, ws, wd, wd, wd, wd;
    }
} 

const int LBM::VelocitySet::getD() { return this->D; }

const int LBM::VelocitySet::getQ() { return this->Q; }

const MatrixXd& LBM::VelocitySet::get_c() { return this->c; }

const VectorXd& LBM::VelocitySet::get_w() { return this->w; }

const int LBM::VelocitySet::fromStdD(StandardSet std)
{
    if(std == D2Q9) return 2;
    //if(std == ...) ...
}

const int LBM::VelocitySet::fromStdQ(StandardSet std)
{
    if(std == D2Q9) return 9;
    //if(std == ...) ...
}