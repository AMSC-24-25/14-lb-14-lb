#include <stdio.h>
#include <stdlib.h>

#include "seconds.h"
#include "LBM.h"
#include <utility>

int main(int argc, char* argv[])
{
    printf("Simulating Lid Driven cavity\n");
    printf("      domain size: %ux%u\n",NX,NY);
    printf("               nu: %g\n",nu);
    printf("              tau: %g\n",tau);
    printf("            u_max: %g\n",u_max);
    printf("             rho0: %g\n",rho0);
    printf("        timesteps: %u\n",NSTEPS);
    printf("       save every: %u\n",NSAVE);
    printf("    message every: %u\n",NMSG);
    printf("\n");
    
    double bytesPerMiB = 1024.0*1024.0;
    double bytesPerGiB = 1024.0*1024.0*1024.0;
    
    double *f0  = (double*) malloc(mem_size_0dir);
    double *f1  = (double*) malloc(mem_size_n0dir);
    double *f2  = (double*) malloc(mem_size_n0dir);
    double *rho = (double*) malloc(mem_size_scalar);
    double *ux  = (double*) malloc(mem_size_scalar);
    double *uy  = (double*) malloc(mem_size_scalar);
    
    size_t total_mem_bytes = mem_size_0dir + 2*mem_size_n0dir + 3*mem_size_scalar;
    
    if(f0 == NULL || f1 == NULL || f2 == NULL || rho == NULL || ux == NULL || uy == NULL)
    {
        fprintf(stderr,"Error: unable to allocate required memory (%.1f MiB).\n",total_mem_bytes/bytesPerMiB);
        exit(-1);
    }
    
    // compute lid-driven cavity flow at t=0 
    // to initialise rho, ux, uy fields.
    lid_driven_cavity(rho, ux, uy);
    
    // initialise f1 as equilibrium for rho, ux, uy
    init_equilibrium(f0,f1,rho,ux,uy);

    // Name of CSV
    const char *csv_filename = "simulation_parameters.csv";

    // Write head of CSV
    FILE *csv_file = fopen(csv_filename, "w");
    if (csv_file != NULL)
    {
        fprintf(csv_file, "NX,NY,NSTEPS,NSAVE,UMAX\n");
        fprintf(csv_file, "%u,%u,%u,%u,%lf\n", NX, NY, NSTEPS, NSAVE, u_max);
        fclose(csv_file);
    }
    else
    {
        fprintf(stderr, "Errore nell'apertura del file CSV %s\n", csv_filename);
    }

    save_scalar("rho", rho, 0);
    save_scalar("ux", ux, 0);
    save_scalar("uy", uy, 0);
    

    if(computeFlowProperties)
    {
        report_flow_properties(0,rho,ux,uy);
    }
    
    double start = seconds();
    #pragma omp parallel master 
    {
// main simulation loop; take NSTEPS time steps
    for(unsigned int n = 0; n < NSTEPS; ++n)
    {
        bool save = (n+1)%NSAVE == 0;
        bool msg  = (n+1)%NMSG == 0;
        bool need_scalars = save || (msg && computeFlowProperties);
        
        // stream and collide from f1 storing to f2
        // optionally compute and save moments
        stream_collide_save(f0,f1,f2,rho,ux,uy,need_scalars);

        
        if(save)
        {
            #pragma omp taskwait
            #pragma omp task 
            {
                    save_scalar("rho", rho, n + 1);
                    save_scalar("ux", ux, n + 1);
                    save_scalar("uy", uy, n + 1);
            }
        }
        
        // swap pointers
        std::swap(f1,f2);
        
        if(msg)
        {
            if(computeFlowProperties)
            {
                report_flow_properties(n+1,rho,ux,uy);
            }
            
            if(!quiet)
                printf("completed timestep %d\n",n+1);
        }
    }
    }
    
    double end = seconds();
    double runtime = end-start;

    size_t doubles_read = ndir; // per node every time step
    size_t doubles_written = ndir;
    size_t doubles_saved = 3; // per node every NSAVE time steps
    
    // note NX*NY overflows when NX=NY=65536
    size_t nodes_updated = NSTEPS*size_t(NX*NY);
    size_t nodes_saved   = (NSTEPS/NSAVE)*size_t(NX*NY);
    double speed = nodes_updated/(1e6*runtime);
    
    double bandwidth = (nodes_updated*(doubles_read + doubles_written)+nodes_saved*(doubles_saved))*sizeof(double)/(runtime*bytesPerGiB);
    
    printf(" ----- performance information -----\n");
    printf(" memory allocated: %.1f (MiB)\n",total_mem_bytes/bytesPerMiB);
    printf("        timesteps: %u\n",NSTEPS);
    printf("          runtime: %.3f (s)\n",runtime);
    printf("            speed: %.2f (Mlups)\n",speed);
    printf("        bandwidth: %.1f (GiB/s)\n",bandwidth);
    
    // deallocate memory
    free(f0);  free(f1); free(f2);
    free(rho); free(ux); free(uy);
    
    return 0;
}

