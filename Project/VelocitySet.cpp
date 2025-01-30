#include "LBM.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

LBM::VelocitySet::VelocitySet(StandardSet set) : D(fromStdD(set)), Q(fromStdQ(set))
{
    c = MatrixXd(Q, D);
    w = VectorXd(Q);

    /*
    6 2 5
    3 0 1
    7 4 8
    */

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
    if(Q == 15){
        c << 0.0, 0.0, 0.0,
             1.0, 0.0, 0.0,
            -1.0, 0.0, 0.0, 
             0.0, 1.0, 0.0,
            0.0, -1.0, 0.0,
             0.0, 0.0, 1.0,
            0.0, 0.0, -1.0,
             1.0, 1.0, 1.0,
          -1.0, -1.0, -1.0,
            1.0, 1.0, -1.0,
           -1.0, -1.0, 1.0,
            1.0, -1.0, 1.0,
           -1.0, 1.0, -1.0,
            -1.0, 1.0, 1.0,
           1.0, -1.0, -1.0;

        constexpr double w0 = 2.0/9.0;  // zero weight
        constexpr double ws = 1.0/9.0;  // adjacent weight
        constexpr double wd = 1.0/72.0;
        w << w0, ws, ws, ws, ws, ws, ws, wd, wd, wd, wd, wd, wd, wd, wd;
    }
    if( Q == 19)
    {
        c << 0.0, 0.0, 0.0,
             1.0, 0.0, 0.0,
            -1.0, 0.0, 0.0, 
             0.0, 1.0, 0.0,
            0.0, -1.0, 0.0,
             0.0, 0.0, 1.0,
            0.0, 0.0, -1.0,
             1.0, 1.0, 0.0,
           -1.0, -1.0, 0.0,
             1.0, 0.0, 1.0,
           -1.0, 0.0, -1.0,
             0.0, 1.0, 1.0,
           0.0, -1.0, -1.0,
            1.0, -1.0, 0.0,
            -1.0, 1.0, 0.0,
            1.0, 0.0, -1.0,
            -1.0, 0.0, 1.0,
            0.0, 1.0, -1.0,
            0.0, -1.0, 1.0; 
             
        constexpr double w0 = 1.0/3.0;  // zero weight
        constexpr double ws = 1.0/18.0;  // adjacent weight
        constexpr double wd = 1.0/36.0;
        w << w0, ws, ws, ws, ws, ws, ws, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd;
    }
    if(Q == 27)
    {
        c << 0.0, 0.0, 0.0,
             1.0, 0.0, 0.0,
            -1.0, 0.0, 0.0, 
             0.0, 1.0, 0.0,
            0.0, -1.0, 0.0,
             0.0, 0.0, 1.0,
            0.0, 0.0, -1.0,
             1.0, 1.0, 0.0,
           -1.0, -1.0, 0.0,
             1.0, 0.0, 1.0,
           -1.0, 0.0, -1.0,
             0.0, 1.0, 1.0,
           0.0, -1.0, -1.0,
            1.0, -1.0, 0.0,
            -1.0, 1.0, 0.0,
            1.0, 0.0, -1.0,
            -1.0, 0.0, 1.0,
            0.0, 1.0, -1.0,
            0.0, -1.0, 1.0,
             1.0, 1.0, 1.0,
          -1.0, -1.0, -1.0,
            1.0, 1.0, -1.0,
           -1.0, -1.0, 1.0,
            1.0, -1.0, 1.0,
           -1.0, 1.0, -1.0,
            -1.0, 1.0, 1.0,
           1.0, -1.0, -1.0;

        constexpr double w0 = 8.0/27.0;  // zero weight
        constexpr double ws = 2.0/27.0;  // adjacent weight
        constexpr double wd = 1.0/54.0;
        constexpr double wv = 1.0/216.0;
        w << w0, ws, ws, ws, ws, ws, ws, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wv, wv, wv, wv, wv, wv, wv, wv;
    }
    

    for(unsigned int i = 0; i < Q; ++i)
    {
        for(unsigned int j = 0; j < Q; ++j)
        {
            if(c.row(i) == -c.row(j))
            {
                inverse[i] = j;
            }
        }
    }
    std::vector<bool> visited(Q, false);

    for(unsigned int i = 0; i < Q; ++i) {
        if (!visited[i]) {
            int inv = inverse[i];
            
            // Mark the current element and its inverse as visited
            visited[i] = true;
            visited[inv] = true;
            
            // Add to basis
            base.push_back(i);
        }
    }
} 

const unsigned int LBM::VelocitySet::getD() { return this->D; }

const unsigned int LBM::VelocitySet::getQ() { return this->Q; }

const MatrixXd& LBM::VelocitySet::get_c() { return this->c; }

const VectorXd& LBM::VelocitySet::get_w() { return this->w; }

const unsigned int LBM::VelocitySet::fromStdD(StandardSet std)
{
    if(std == D1Q3) return 1;
    else if(std == D2Q9) return 2;
    else if(std == D3Q15) return 3;
    else if(std == D3Q19) return 3;
    else if(std == D3Q27) return 3;
    else return 0;
}

const unsigned int LBM::VelocitySet::fromStdQ(StandardSet std)
{
    if(std == D1Q3) return 3;
    else if(std == D2Q9) return 9;
    else if(std == D3Q15) return 15;
    else if(std == D3Q19) return 19;
    else if(std == D3Q27) return 27;
    else return 0;
}

const unsigned int LBM::VelocitySet::directionIndexInvert(unsigned int i)
{
    return inverse[i];
}
const unsigned int LBM::VelocitySet::directionIndexBase(unsigned int i)
{
    return base[i];
}
