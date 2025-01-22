#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "LBM.h"
#include <omp.h>

inline void lid_driven_cavity(unsigned int x, unsigned int y, double &r, double &u, double &v) {

    
    r = rho0; //Set initial density constant
    
    // Inizialize velocity
    if (y == NY - 1) // Top lid
    {
        u = u_max; // Constant value on top lid
    }
    else // All other points
    {
        u = 0.0;
    }

    v = 0.0;
}


void lid_driven_cavity(double *r, double *u, double *v)
{
    #pragma omp parallel for default(none) shared(r,u,v) schedule(static)
    for (unsigned int y = 0; y < NY; ++y)
    {
        for (unsigned int x = 0; x < NX; ++x)
        {
            size_t sidx = scalar_index(x, y);
            lid_driven_cavity(x, y, r[sidx], u[sidx], v[sidx]);
        }
    }
}


void init_equilibrium(double *f0, double *f1, double *r, double *u, double *v)
{
    #pragma omp parallel for default(none) shared(f0,f1,r,u,v) schedule(static)
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            size_t si = scalar_index(x, y);
            double rho = r[si];
            double ux = u[si];
            double uy = v[si];

            // temporary variables
            double w0r = w0*rho;
            double wsr = ws*rho;
            double wdr = wd*rho;
            double omusq = 1.0 - 1.5*(ux*ux+uy*uy);
            
            double tux = 3.0*ux;
            double tuy = 3.0*uy;
            
            f0[field0_index(x,y)]    = w0r*(omusq);
            
            double cidot3u = tux;
            f1[fieldn_index(x,y,1)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy;
            f1[fieldn_index(x,y,2)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tux;
            f1[fieldn_index(x,y,3)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tuy;
            f1[fieldn_index(x,y,4)]  = wsr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            
            cidot3u = tux+tuy;
            f1[fieldn_index(x,y,5)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy-tux;
            f1[fieldn_index(x,y,6)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -(tux+tuy);
            f1[fieldn_index(x,y,7)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tux-tuy;
            f1[fieldn_index(x,y,8)]  = wdr*(omusq + cidot3u*(1.0+0.5*cidot3u));
        }
    }
}

void stream_collide_save(double *f0, double *f1, double *f2, double *r, double *u, double *v, bool save)
{
    #pragma omp parallel for default(none) \
       shared(f0,f1,f2,r,u,v,save) schedule(static)
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            unsigned int xp1 = (x+1)%NX;
            unsigned int yp1 = (y+1)%NY;
            unsigned int xm1 = (NX+x-1)%NX;
            unsigned int ym1 = (NY+y-1)%NY;
            
            // direction numbering scheme
            // 6 2 5
            // 3 0 1
            // 7 4 8

            double ft[9];
            
            ft[0] = f0[field0_index(x,y)];
            
            // load populations from adjacent nodes
            ft[1] = f1[fieldn_index(xm1,y,  1)];
            ft[2] = f1[fieldn_index(x,  ym1,2)];
            ft[3] = f1[fieldn_index(xp1,y,  3)];
            ft[4] = f1[fieldn_index(x,  yp1,4)];
            ft[5] = f1[fieldn_index(xm1,ym1,5)];
            ft[6] = f1[fieldn_index(xp1,ym1,6)];
            ft[7] = f1[fieldn_index(xp1,yp1,7)];
            ft[8] = f1[fieldn_index(xm1,yp1,8)];

             // Bounce back on Left wall
            if( x == 0)
            {
                ft[1] = f1[fieldn_index(x, y, 3)];
                ft[5] = f1[fieldn_index(x, y, 7)];
                ft[8] = f1[fieldn_index(x, y, 6)];
            }
            
            // Bounce back on Right wall
            if (x == NX - 1)
            { 
                ft[3] = f1[fieldn_index(x, y, 1)];
                ft[7] = f1[fieldn_index(x, y, 5)];
                ft[6] = f1[fieldn_index(x, y, 8)];
            }

            // Bounce back on Bottom wall

            if (y == 0)
            {
            ft[2] = f1[fieldn_index(x, y, 4)];
            ft[5] = f1[fieldn_index(x, y, 7)];
            ft[6] = f1[fieldn_index(x, y, 8)];
            }

            //  Top lid (moving wall)
            if (y == NX - 1)
            {
            double u_lid = u_max;
            double rho_top = r[scalar_index(x, y)];
        
            // Conditions for top lid

            double coeff = 2.0 * wd * (1.0/cs) * (1.0/cs) * rho_top;
            
            ft[4] = f1[fieldn_index(x, y, 2)];
            ft[7] = f1[fieldn_index(x, y, 5)] - coeff * u_lid;
            ft[8] = f1[fieldn_index(x, y, 6)] + coeff * u_lid;
            }            
            // compute moments
            double rho = ft[0]+ft[1]+ft[2]+ft[3]+ft[4]+ft[5]+ft[6]+ft[7]+ft[8];
            double rhoinv = 1.0/rho;
            
            double ux = rhoinv*(ft[1]+ft[5]+ft[8]-(ft[3]+ft[6]+ft[7]));
            double uy = rhoinv*(ft[2]+ft[5]+ft[6]-(ft[4]+ft[7]+ft[8]));
            
            // ALWAYS write to memory when needed
           
                r[scalar_index(x,y)] = rho;
                u[scalar_index(x,y)] = ux;
                v[scalar_index(x,y)] = uy;
            
            // temporary variables
            double omusq = 1.0 - 1.5*(ux*ux+uy*uy); // 1-(3/2)u.u
            
            double tux = 3.0*ux;
            double tuy = 3.0*uy;

            double feq[9];
            feq[0] = w0*rho*omusq;
            double cidot3u = tux;
            feq[1] = ws*rho*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy;
            feq[2] = ws*rho*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tux;
            feq[3] = ws*rho*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -tuy;
            feq[4] = ws*rho*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tux+tuy;
            feq[5] = wd*rho*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tuy-tux;
            feq[6] = wd*rho*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = -(tux+tuy);
            feq[7] = wd*rho*(omusq + cidot3u*(1.0+0.5*cidot3u));
            cidot3u = tux-tuy;
            feq[8] = wd*rho*(omusq + cidot3u*(1.0+0.5*cidot3u));

            double fSym = ft[0];
            double fAnti = 0;
            double feqSym = feq[0];
            double feqAnti = 0;
            double symmetric = omgSym*(fSym-feqSym);
            double antisymmetric = 0;
            f0[field0_index(x,y)]    = ft[0] - symmetric;
            
            int idx[4] = {1, 2, 5, 6};
            for (int j = 0; j < 4; j++) {
                int i = idx[j];
                int inv = invert_index(i);
                fSym = (ft[i] + ft[inv])/2;
                fAnti = (ft[i] - ft[inv])/2;
                feqSym = (feq[i] + feq[inv])/2;
                feqAnti = (feq[i] - feq[inv])/2;
                symmetric = omgSym*(fSym-feqSym);
                antisymmetric = omgAnti*(fAnti - feqAnti);
                f2[fieldn_index(x,y,i)] = ft[i] - symmetric - antisymmetric;

                f2[fieldn_index(x,y,inv)] = ft[inv] - symmetric + antisymmetric;
            }
        }
    }
}

void compute_flow_properties(unsigned int t, double *r, double *u, double *v, double *prop)
{
    // prop must point to space for 4 doubles:
    // 0: energy
    // 1: L2 error in rho
    // 2: L2 error in ux
    // 3: L2 error in uy
    
    double E = 0.0; // kinetic energy
    
    double sumrhoe2 = 0.0; // sum of error squared in rho
    double sumuxe2 = 0.0;  //                         ux
    double sumuye2 = 0.0;  //                         uy
    
    double sumrhoa2 = 0.0; // sum of analytical rho squared
    double sumuxa2 = 0.0;  //                   ux
    double sumuya2 = 0.0;  //                   uy
    
    #pragma omp parallel for default(none) shared(t,r,u,v) \
        reduction(+:E,sumrhoe2,sumuxe2,sumuye2, \
                      sumrhoa2,sumuxa2,sumuya2) \
        schedule(static)
    for(unsigned int y = 0; y < NY; ++y)
    {
        for(unsigned int x = 0; x < NX; ++x)
        {
            double rho = r[scalar_index(x,y)];
            double ux  = u[scalar_index(x,y)];
            double uy  = v[scalar_index(x,y)];
            E += rho*(ux*ux + uy*uy);
            
            double rhoa, uxa, uya;
            lid_driven_cavity(x, y, rhoa, uxa, uya);
            
            sumrhoe2 += (rho-rhoa)*(rho-rhoa);
            sumuxe2  += (ux-uxa)*(ux-uxa);
            sumuye2  += (uy-uya)*(uy-uya);

            sumrhoa2 += (rhoa-rho0)*(rhoa-rho0);
            sumuxa2  += uxa*uxa;
            sumuya2  += uya*uya;
        }
    }
    
    prop[0] = E;
    prop[1] = sqrt(sumrhoe2/sumrhoa2);
    prop[2] = sqrt(sumuxe2/sumuxa2);
    prop[3] = sqrt(sumuye2/sumuya2);
}

void report_flow_properties(unsigned int t, double *rho, double *ux, double *uy)
{
    double prop[4];
    compute_flow_properties(t,rho,ux,uy,prop);
    printf("%u,%g,%g,%g,%g\n",t,prop[0],prop[1],prop[2],prop[3]);
}

void save_scalar(const char* name, double *scalar, unsigned int n)
{
    // assume reasonably-sized file names
    char filename[128];
    char format[32];
    
    // compute maximum number of digits
    int ndigits = floor(log10((double)NSTEPS)+1.0);
    
    // generate format string
    // file name format is name0000nnn.bin
    sprintf(format,"../bin_results/%%s%%0%dd.bin",ndigits);
    sprintf(filename,format,name,n);
    
    // open file for writing
    FILE *fout = fopen(filename, "wb+");
    
    // write data
    fwrite(scalar,1,mem_size_scalar,fout);
    
    // close file
    fclose(fout);
    
    if(ferror(fout) && !quiet)
    {
        fprintf(stderr,"Error saving to %s\n",filename);
        perror("");
    }
    else
    {
        if(!quiet)
            printf("Saved to %s\n",filename);
    }
}
