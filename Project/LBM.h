#ifndef __LBM_H
#define __LBM_H

const unsigned long scale = 2;
const unsigned long NX = 1024 * scale;
const unsigned long NY = NX;

const unsigned long ndir = 9;
const long mem_size_0dir = sizeof(double) * NX * NY;
const long mem_size_n0dir = sizeof(double) * NX * NY * (ndir - 1);
const long mem_size_scalar = sizeof(double) * NX * NY;

// Size for the populations and scalar values
const long size_0dir = NX * NY;
const long size_n0dir = NX * NY * (ndir - 1);
const long size_scalar = NX * NY;
const double w0 = 4.0 / 9.0;  // zero weight
const double ws = 1.0 / 9.0;  // adjacent weight
const double wd = 1.0 / 36.0; // diagonal weight

// Re = u*N/(cs^2 * (tau - 1/2))
// cs^2 = 1/3
const double Re = 100.0;

// kinematic viscosity nu and the corresponding relaxation parameter tau
const double nu = 1.0 / 6.0;
const double tau = 3.0 * nu + 0.5;
const double cs = 1.0 / 1.732;

// Having large u_max destabilizes the simulation, higher NX is better
const double u_max = Re / (6 * NX);
const double rho0 = 1.0;

const unsigned long NSTEPS = 3000 * scale * scale;
const unsigned long NSAVE = 100 * scale * scale;
const unsigned long NMSG = 50 * scale * scale;

// compute L2 error and energy?
// disable for speed testing
const bool computeFlowProperties = false;

// suppress verbose output
const bool quiet = true;

__global__ void lid_driven_cavity(double *r, double *u, double *v);
__global__ void stream_collide_save(double *, double *, double *, double *, double *, double *, bool);
__global__ void init_equilibrium(double *, double *, double *, double *, double *);
__global__ void compute_flow_properties(unsigned long, double *, double *, double *, double *);
void report_flow_properties(unsigned long, double *, double *, double *);
void apply_lid_boundary(double *f1, double *rho, double u_lid);
void apply_bounce_back(double *f1);
void save_scalar(const char *name, double *scalar, unsigned long n);

#endif /* __LBM_H */
