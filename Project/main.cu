#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "seconds.h"
#include "LBM.h"
#include <utility>

int main(int argc, char *argv[])
{
    printf("Simulating Lid Driven cavity\n");
    printf("      domain size: %ux%u\n", NX, NY);
    printf("               nu: %g\n", nu);
    printf("              tau: %g\n", tau);
    printf("            u_max: %g\n", u_max);
    printf("             rho0: %g\n", rho0);
    printf("        timesteps: %u\n", NSTEPS);
    printf("       save every: %u\n", NSAVE);
    printf("    message every: %u\n", NMSG);
    printf("\n");

    double bytesPerMiB = 1024.0 * 1024.0;
    double bytesPerGiB = 1024.0 * 1024.0 * 1024.0;

    double *f0;
    double *f1;
    double *f2;
    double *rho;
    double *ux;
    double *uy;

    cudaError_t f0_error = cudaMallocManaged(&f0, mem_size_0dir);
    cudaError_t f1_error = cudaMallocManaged(&f1, mem_size_n0dir);
    cudaError_t f2_error = cudaMallocManaged(&f2, mem_size_n0dir);
    cudaError_t rho_error = cudaMallocManaged(&rho, mem_size_scalar);
    cudaError_t ux_error = cudaMallocManaged(&ux, mem_size_scalar);
    cudaError_t uy_error = cudaMallocManaged(&uy, mem_size_scalar);

    long total_mem_bytes = mem_size_0dir + 2 * mem_size_n0dir + 3 * mem_size_scalar;

    if (
        f0_error != cudaSuccess ||
        f1_error != cudaSuccess ||
        f2_error != cudaSuccess ||
        rho_error != cudaSuccess ||
        ux_error != cudaSuccess ||
        uy_error != cudaSuccess)
    {
        fprintf(stderr, "Error: unable to allocate required memory (%.1f MiB).\n", total_mem_bytes / bytesPerMiB);
        exit(-1);
    }
    else
    {
        printf("Allocated required memory (%.1f MiB).\n", total_mem_bytes / bytesPerMiB);
    }

    // compute lid-driven cavity flow at t=0
    // to initialise rho, ux, uy fields.
    lid_driven_cavity<<<3072, 64>>>(rho, ux, uy);
    cudaDeviceSynchronize();
    
    // initialise f1 as equilibrium for rho, ux, uy
    init_equilibrium<<<3072, 64>>>(f0, f1, rho, ux, uy);
    cudaDeviceSynchronize();
    
    // Name of CSV
    const char *csv_filename = "simulation_parameters.csv";

    // Write head of CSV
    FILE *csv_file = fopen(csv_filename, "w");
    if (csv_file != nullptr)
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

    if (computeFlowProperties)
    {
        report_flow_properties(0, rho, ux, uy);
    }

    double start = seconds();

    // main simulation loop; take NSTEPS time steps
    for (unsigned long n = 0; n < NSTEPS; ++n)
    {
        bool save = (n + 1) % NSAVE == 0;
        bool msg = (n + 1) % NMSG == 0;
        bool need_scalars = save || (msg && computeFlowProperties);

        // stream and collide from f1 storing to f2
        // optionally compute and save moments
        stream_collide_save<<<3072, 64>>>(f0, f1, f2, rho, ux, uy, need_scalars);
        cudaDeviceSynchronize();
        if (save)
        {

            save_scalar("rho", rho, n + 1);
            save_scalar("ux", ux, n + 1);
            save_scalar("uy", uy, n + 1);
        }

        // swap pointers
        std::swap(f1, f2);

        if (msg)
        {
            if (computeFlowProperties)
            {
                report_flow_properties(n + 1, rho, ux, uy);
            }

            if (!quiet)
                printf("completed timestep %d\n", n + 1);
        }
    }

    double end = seconds();
    double runtime = end - start;

    long doubles_read = ndir; // per node every time step
    long doubles_written = ndir;
    long doubles_saved = 3; // per node every NSAVE time steps

    // note NX*NY overflows when NX=NY=65536
    long nodes_updated = NSTEPS * long(NX * NY);
    long nodes_saved = (NSTEPS / NSAVE) * long(NX * NY);
    double speed = nodes_updated / (1e6 * runtime);

    double bandwidth = (nodes_updated * (doubles_read + doubles_written) + nodes_saved * (doubles_saved)) * sizeof(double) / (runtime * bytesPerGiB);

    printf(" ----- performance information -----\n");
    printf(" memory allocated: %.1f (MiB)\n", total_mem_bytes / bytesPerMiB);
    printf("        timesteps: %u\n", NSTEPS);
    printf("          runtime: %.3f (s)\n", runtime);
    printf("            speed: %.2f (Mlups)\n", speed);
    printf("        bandwidth: %.1f (GiB/s)\n", bandwidth);

    // deallocate memory
    cudaFree(f0);
    cudaFree(f1);
    cudaFree(f2);
    cudaFree(rho);
    cudaFree(ux);
    cudaFree(uy);

    return 0;
}