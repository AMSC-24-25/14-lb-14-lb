/* This code accompanies
 *   The Lattice Boltzmann Method: Principles and Practice
 *   T. Krüger, H. Kusumaatmaja, A. Kuzmin, O. Shardt, G. Silva, E.M. Viggen
 *   ISBN 978-3-319-44649-3 (Electronic) 
 *        978-3-319-44647-9 (Print)
 *   http://www.springer.com/978-3-319-44647-9
 *
 * This code is provided under the MIT license. See LICENSE.txt.
 *
 * Author: Orest Shardt
 *
 */
#include <tuple>
#include <Eigen/Dense>

#ifndef __LBM_H
#define __LBM_H

class LBM 
{

    private:
        Velocity v;
        unsigned int scale;
        struct N {
            int x;
            int y;
            int z;
        } N;
        const double nu = 1.0/6.0;
        const double tau = 3.0*nu+0.5;
        const double cs = 1.0/1.732;

        bool computeFlowProperties;
        bool quiet;
        int step;

        double *population;
        double *rho;
        double *u;

    void init_equilibrium();
    void stream_collide_save();

};

#endif /* __LBM_H */

