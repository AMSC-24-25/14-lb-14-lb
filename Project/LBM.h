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

#include <memory>

#ifndef __LBM_H
#define __LBM_H

const bool bTest = false;

const unsigned int scale = 2;
const unsigned int NX = 32 * scale;
const unsigned int NY = NX;

// The number of directions in the lattice
const unsigned int ndir = 9;

const size_t size_0dir = NX * NY;
const size_t size_n0dir = NX * NY * (ndir - 1);
const size_t size_scalar = NX * NY;

const double w0 = 4.0 / 9.0;  // zero weight
const double ws = 1.0 / 9.0;  // adjacent weight
const double wd = 1.0 / 36.0; // diagonal weight

// Re = u*N/(cs^2 * (tau - 1/2))
// cs^2 = 1/3
const double Re = 100.0;

const double nu = 1.0 / 6.0;
const double tau = 3.0 * nu + 0.5;
const double cs = 1.0 / 1.732;

// Having large u_max destabilizes the simulation, higher NX is better
const double u_max = Re / (6 * NX);
const double rho0 = 1.0;

const unsigned int NSTEPS = 1000 * scale * scale;
const unsigned int NSAVE = 50 * scale * scale;
const unsigned int NMSG = 50 * scale * scale;

// compute L2 error and energy?
// disable for speed testing
const bool computeFlowProperties = false;

// suppress verbose output
const bool quiet = true;

void lid_driven_cavity(std::unique_ptr<double[]> &r, std::unique_ptr<double[]> &u, std::unique_ptr<double[]> &v);
void stream_collide_save(std::unique_ptr<double[]> &f0, std::unique_ptr<double[]> &f1, std::unique_ptr<double[]> &f2, std::unique_ptr<double[]> &r, std::unique_ptr<double[]> &u, std::unique_ptr<double[]> &v, bool save);
void init_equilibrium(std::unique_ptr<double[]> &f0, std::unique_ptr<double[]> &f1, std::unique_ptr<double[]> &r, std::unique_ptr<double[]> &u, std::unique_ptr<double[]> &v);
void compute_flow_properties(unsigned int t, std::unique_ptr<double[]> &r, std::unique_ptr<double[]> &u, std::unique_ptr<double[]> &v, std::unique_ptr<double[]> &prop);
void report_flow_properties(unsigned int t, std::unique_ptr<double[]> &rho, std::unique_ptr<double[]> &ux, std::unique_ptr<double[]> &uy);
void apply_lid_boundary(std::unique_ptr<double[]> &f1, std::unique_ptr<double[]> &rho, double u_lid);
void apply_bounce_back(std::unique_ptr<double[]> &f1);
void save_to_csv(const char *filename, unsigned int t, std::unique_ptr<double[]> &rho, std::unique_ptr<double[]> &ux, std::unique_ptr<double[]> &uy);

inline size_t field0_index(unsigned int x, unsigned int y)
{
    return NX * y + x;
}

inline size_t scalar_index(unsigned int x, unsigned int y)
{
    return NX * y + x;
}

inline size_t fieldn_index(unsigned int x, unsigned int y, unsigned int d)
{
    return (ndir - 1) * (NX * y + x) + (d - 1);
}

#endif /* __LBM_H */
