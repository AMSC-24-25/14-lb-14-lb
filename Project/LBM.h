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
#ifndef __LBM_H
#define __LBM_H

const unsigned int scale = 2;
const unsigned int NX = 2048*scale;
const unsigned int NY = NX;

const unsigned int ndir = 9;
const size_t mem_size_0dir   = NX*NY;
const size_t mem_size_n0dir  = NX*NY*(ndir-1);
const size_t mem_size_scalar = NX*NY;

const double w0 = 4.0/9.0;  // zero weight
const double ws = 1.0/9.0;  // adjacent weight
const double wd = 1.0/36.0; // diagonal weight

const double nu = 1.0/6.0;
const double tau = 3.0*nu+0.5;

// Taylor-Green parameters
const double u_max = 0.04/scale;
const double rho0 = 1.0;

const unsigned int NSTEPS = 200*scale*scale;
const unsigned int NSAVE  =  50*scale*scale;
const unsigned int NMSG   =  50*scale*scale;

// compute L2 error and energy?
// disable for speed testing
const bool computeFlowProperties = false;

// suppress verbose output
const bool quiet = true;

void lid_driven_cavity(std::vector<double> &r, std::vector<double> &u, std::vector<double> &v);
void stream_collide_save(std::vector<double>&f0, std::vector<double>&f1, std::vector<double>&f2, std::vector<double>&r, std::vector<double> &u, std::vector<double>&v, bool save);
void init_equilibrium(std::vector<double>&f0, std::vector<double>&f1, std::vector<double> &r, std::vector<double> &u, std::vector<double>&v);
void compute_flow_properties(unsigned int t, std::vector<double>&r, std::vector<double> &u, std::vector<double> &v, std::vector<double> &prop);
void report_flow_properties(unsigned int t, std::vector<double>&rho, std::vector<double> &ux, std::vector<double> &uy);
void apply_lid_boundary(std::vector<double>&f1, std::vector<double>&rho, double u_lid);
void apply_bounce_back(std::vector<double>&f1);
void save_to_csv(const char* filename, unsigned int t, std::vector<double>&rho, std::vector<double> &ux, std::vector<double> &uy);

inline size_t field0_index(unsigned int x, unsigned int y)
{
    return NX*y+x;
}

inline size_t scalar_index(unsigned int x, unsigned int y)
{
    return NX*y+x;
}

inline size_t fieldn_index(unsigned int x, unsigned int y, unsigned int d)
{
    return (ndir-1)*(NX*y+x)+(d-1);
}

#endif /* __LBM_H */

