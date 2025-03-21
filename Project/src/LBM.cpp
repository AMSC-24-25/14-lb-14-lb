#include <LatticeBoltzmannMethod/LBM.hpp>
#include <iostream>
#include <filesystem> 
#include <iomanip>
#include <fstream>
#include <csignal>
#include <mpi.h>

using Eigen::VectorXd;
using Eigen::ArrayXd;

namespace LatticeBoltzmannMethod{

    void LBM::init_equilibrium()
    {
        if (this->node_id == 0){
            std::cout << "Initializing simulation..." << std::endl;
        }

        #pragma omp parallel for default(none) schedule(static)
        for (unsigned int x = x_loc.start; x < x_loc.end; ++x)
        {
            for(unsigned int y = 0; y < this->N.y; ++y)
            {
                for(unsigned int z = 0; z < this->N.z; ++z)
                {
                    applyInitial(x,y,z);

                    double rho = get_rho(x,y,z);
                    VectorXd u = get_u(x,y,z);
                    double unorm = u.norm();
                    double omusq = 1.0 - 1.5*(unorm * unorm);
                    ArrayXd c3u = this->v->get_c() * (u * 3.0);
                    VectorXd f = (rho * this->v->get_w()).array() * (1.0 * omusq + c3u * (c3u*0.5 + 1.0));

                    savePopulationInit(x,y,z, f);
                }
            }
        }

        if (this->node_id == 0){
            std::cout << "Initialization completed." << std::endl;
        }        
    }

    void LBM::stream_collide_save()
    {
        if (this->node_id == 0){
            std::cout << "Simulating step " << this->step << "...  ";
        }
        // useful constants

        // sync with other nodes
        int population_offset = x_len*N.y*N.z * v->getQ() * (step & 1);
        int overlap_size = N.y * N.z * v->getQ();
        double *right_p = population.get() + population_offset + index_f(x_loc.end-1, 0, 0); // last slice of computed population
        double *left_p  = population.get() + population_offset + index_f(x_loc.start, 0, 0); // first slice of computed population
        MPI_Request request;
        
        if (x_loc.left_pad)
            MPI_Isend(left_p, overlap_size, MPI_DOUBLE, node_id - 1, 0, MPI_COMM_WORLD, &request);
        if (x_loc.right_pad)
            MPI_Isend(right_p, overlap_size, MPI_DOUBLE, node_id + 1, 0, MPI_COMM_WORLD, &request);  
        
        //@note this is a blocking operation while the previous one is not
        if (x_loc.left_pad)
            MPI_Recv(left_p - overlap_size, overlap_size, MPI_DOUBLE, node_id - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (x_loc.right_pad)
            MPI_Recv(right_p + overlap_size, overlap_size, MPI_DOUBLE, node_id + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
        
        #pragma omp parallel for default(none) schedule(static)
        for (unsigned int x = x_loc.start; x < x_loc.end; ++x)
        {
            for(unsigned int y = 0; y < this->N.y; ++y)
            {
                for(unsigned int z = 0; z < this->N.z; ++z)
                {
                    //classical stream
                    VectorXd stream = populationAdjacent(x,y,z);
                    //modify stream applying boundaries
                    applyBoundaryAndObstacle(x,y,z, stream);

                        //compute moments
                    double rho = stream.sum();
                    set_rho(x,y,z, rho );

                    double rhoinv = 1.0/rho;

                    VectorXd cf_rhoinv = (this->v->get_c().transpose() * stream) * rhoinv;
                    VectorXd& u = cf_rhoinv;
                    set_u(x,y,z, u);

                    //collision
                    double unorm = u.norm();
                    double omusq = 1.0 - 1.5*(unorm * unorm);
                    ArrayXd c3u = this->v->get_c() * (u * 3.0);
                    VectorXd twr = this->v->get_w() * rho;

                    int dim = v->getQ()/2+1;

                    VectorXd feq = (twr.array() * ( omusq + c3u * (1.0 + 0.5*c3u) ));
                    VectorXd ftsym(dim);
                    VectorXd ftanti(dim);
                    VectorXd feqsym(dim);
                    VectorXd feqanti(dim);

                    for (int i = 0; i < dim; i++) {
                        int idx = v->directionIndexBase(i);
                        int inv = v->directionIndexInvert(idx);
                        ftsym[i] = stream[idx];
                        ftanti[i] = stream[inv];
                        feqsym[i] = feq[idx];
                        feqanti[i] = feq[inv];
                    }

                    VectorXd sym = omgSym * (ftsym + ftanti - feqsym - feqanti) / 2.0;
                    VectorXd anti = omgAnti * (ftsym - ftanti - feqsym + feqanti) / 2.0;

                    for (int i = 0; i < dim; i++) {
                        int idx = v->directionIndexBase(i);
                        int inv = v->directionIndexInvert(idx);
                        stream[idx] -= sym[i] + anti[i];
                        stream[inv] -= sym[i] - anti[i];
                    }

                    //stream *= omtauinv;
                    //stream += (twr.array() * ( omusq + c3u * (1.0 + 0.5*c3u) )).matrix(); 

                    this->savePopulation(x,y,z, stream);
                }
            }
        }

        this->step++;
        if (this->node_id == 0){
            std::cout << " Completed." << std::endl;
        }
    }


    LBM::LBM(int world_rank, int world_size, LBM::VelocitySet::StandardSet vSet, LBM::dimensions d, double nu) : N(d), nu(nu)
    {
        v = std::make_unique<VelocitySet>(vSet);
        if(N.x == 0 || N.y == 0 || N.z == 0)
        {
            throw std::invalid_argument(std::string{"Setting problem size to zero!!\n"});
        }
        if((v->getD() == 1 && (N.y != 1 || N.z != 1)) || (v->getD() == 2 && N.z != 1))
        {
            throw std::invalid_argument(std::string{"Incompatible size and velocity set dimension!!\n"});
        }

        double partition_x_size = double(N.x) / world_size;
        x_loc = {std::lround(partition_x_size * world_rank), std::lround(partition_x_size * (world_rank + 1)), world_rank != 0, world_rank != world_size - 1};
        x_len = x_loc.end - x_loc.start + x_loc.left_pad + x_loc.right_pad;
        x_len_no_pad = x_loc.end - x_loc.start;
        node_id = world_rank;


        population = std::make_unique<double[]>(x_len * N.y * N.z * v->getQ() * 2);
        rho = std::make_unique<double[]>(x_len * N.y * N.z);
        u = std::make_unique<double[]>(x_len * N.y * N.z * v->getD());
        
        if(!quiet && world_rank == 0)
        {
            std::cout << "\n------------------------------------------------------" << std::endl;
            std::cout << "Lattice Boltzmann Simulation configured with a " << N.x << "x" << N.y << "x" << N.z << " lattice." << std::endl;
            std::cout << "The velocity has " << v->getQ() << " components." << std::endl;
            std::cout << "Kinematic viscosity was set to " << nu << "." << std::endl;
            std::cout << "------------------------------------------------------\n" << std::endl; 
        }
    }

    void LBM::applyInitial(unsigned int x, unsigned int y, unsigned int z)
    {
        InitialCondition(x,y,z, *this );
    }

    void LBM::applyBoundaryAndObstacle(unsigned int x, unsigned int y, unsigned int z, VectorXd& f)
    {
        for(std::function<void (unsigned int, unsigned int, unsigned int, Eigen::VectorXd &, LBM &)> bc : BoundaryConditions) bc(x,y,z, f, *this);
        for(auto o : obstacles) (*o)(x, y, z, f, *this);
    }

    void LBM::setInitialCondition(std::function<void(unsigned int,unsigned int,unsigned int, LBM&)> initial_condition)
    {
        this->InitialCondition = initial_condition;
    }

    void LBM::addBoundaryCondition(std::function<void(unsigned int,unsigned int,unsigned int, Eigen::VectorXd&, LBM&)> boundary)
    {
        this->BoundaryConditions.push_back(boundary);
    }

    void LBM::addObstacle(std::shared_ptr<Obstacle> o)
    {
        this->obstacles.push_back(o);
    }

    void LBM::saveToBin(unsigned int step)
    {
        const std::string result_base_dir = "./results/bin_results/partition_" + std::to_string(x_loc.start) + "-" + std::to_string(x_loc.end) + "/";
        std::filesystem::create_directories(result_base_dir);
        
        // Allocate arrays for the entire domain
        std::vector<double> u_x(x_len_no_pad * N.y * N.z, 0.0);
        std::vector<double> u_y(x_len_no_pad * N.y * N.z, 0.0);
        std::vector<double> u_z(x_len_no_pad * N.y * N.z, 0.0);
        unsigned int idx = 0;

        for (unsigned int z = 0; z < N.z; ++z)
        {
            for (unsigned int y = 0; y < N.y; ++y)
            {
                for (unsigned int x = x_loc.start; x < x_loc.end; ++x)
                {
                    unsigned int u_index = index_u(x, y, z);

                        u_x[idx] = this->u[u_index];
                        if (v->getD() > 1)
                        {
                            u_y[idx] = this->u[u_index + 1];
                        }
                        if (v->getD() == 3)
                        {
                            u_z[idx] = this->u[u_index + 2];
                        }
                        idx++;
                    }
                }
            }

        // Save u_x to a binary file
        std::ofstream file_u_x(result_base_dir + "/u_x_" + std::to_string(step) + ".bin", std::ios::binary);
        file_u_x.write(reinterpret_cast<char*>(u_x.data()), u_x.size() * sizeof(double));
        file_u_x.close();

        // Save u_y to a binary file if applicable
        if (v->getD() > 1)
        {
            std::ofstream file_u_y(result_base_dir + "/u_y_" + std::to_string(step) + ".bin", std::ios::binary);
            file_u_y.write(reinterpret_cast<char*>(u_y.data()), u_y.size() * sizeof(double));
            file_u_y.close();
        }

        // Save u_z to a binary file if applicable
        if (v->getD() == 3)
        {
            std::ofstream file_u_z(result_base_dir + "/u_z_" + std::to_string(step) + ".bin", std::ios::binary);
            file_u_z.write(reinterpret_cast<char*>(u_z.data()), u_z.size() * sizeof(double));
            file_u_z.close();
        }
    }
}