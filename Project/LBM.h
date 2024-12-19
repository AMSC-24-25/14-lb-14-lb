#ifndef __LBM_H
#define __LBM_H

const unsigned int scale = 3;
const unsigned int NX = 128*scale;
const unsigned int NY = NX;

const unsigned int ndir = 9;
const size_t mem_size_0dir   = sizeof(double)*NX*NY;
const size_t mem_size_n0dir  = sizeof(double)*NX*NY*(ndir-1);
const size_t mem_size_scalar = sizeof(double)*NX*NY;

const double w0 = 4.0/9.0;  // zero weight
const double ws = 1.0/9.0;  // adjacent weight
const double wd = 1.0/36.0; // diagonal weight

//Re = u*N/(cs^2 * (tau - 1/2))
//cs^2 = 1/3
const double Re = 100.0;

const double nu = 1.0/6.0;
const double tau = 3.0*nu+0.5;
const double cs = 1.0/1.732;

//Having large u_max destabilizes the simulation, higher NX is better
const double u_max = Re/(6*NX);
const double rho0 = 1.0;

const unsigned int NSTEPS = 2000*scale*scale;
const unsigned int NSAVE  =  50*scale*scale;
const unsigned int NMSG   =  50*scale*scale;

// compute L2 error and energy?
// disable for speed testing
const bool computeFlowProperties = false;

// suppress verbose output
const bool quiet = true;

void lid_driven_cavity(double *r, double *u, double *v);
void stream_collide_save(double*,double*,double*,double*,double*,double*,bool);
void init_equilibrium(double*,double*,double*,double*,double*);
void compute_flow_properties(unsigned int,double*,double*,double*,double*);
void report_flow_properties(unsigned int,double*,double*,double*);
void apply_lid_boundary(double *f1, double *rho, double u_lid);
void apply_bounce_back(double *f1);
void save_scalar(const char* name, double *scalar, unsigned int n);

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

