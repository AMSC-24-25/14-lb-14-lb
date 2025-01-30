#include "LBM.hpp"
#include <iostream>
#include <filesystem> 
#include <iomanip>
#include <fstream>
#include <csignal>

using Eigen::VectorXd;
using Eigen::ArrayXd;

void LBM::init_equilibrium()
{
    std::cout << "Initializing simulation..." << std::endl;

    #pragma omp parallel for default(none) schedule(static)
    for(unsigned int x = 0; x < this->N.x; ++x)
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

    std::cout << "Initialization completed." << std::endl;
}

void LBM::stream_collide_save()
{
    std::cout << "Simulating step " << this->step << "...  ";
        // useful constants
    const double tauinv = 2.0/(6.0*this->nu+1.0); // 1/tau
    const double omtauinv = 1.0-tauinv;     // 1 - 1/tau

    #pragma omp parallel for default(none) shared(tauinv, omtauinv) schedule(static)
    for(unsigned int x = 0; x < this->N.x; ++x)
    {
        for(unsigned int y = 0; y < this->N.y; ++y)
        {
            for(unsigned int z = 0; z < this->N.z; ++z)
            {
                //classical stream
                VectorXd stream = populationAdjacent(x,y,z);
                //modify stream applying boundaries
                applyBoundary(x,y,z, stream);

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
                VectorXd twr = this->v->get_w() * tauinv * rho;

                stream *= omtauinv;
                stream += (twr.array() * ( omusq + c3u * (1.0 + 0.5*c3u) )).matrix(); 

                this->savePopulation(x,y,z, stream);
            }
        }
    }

    this->step++;
    std::cout << " Completed." << std::endl;
}


LBM::LBM(LBM::VelocitySet::StandardSet vSet, LBM::dimensions d,  double nu) : N(d), nu(nu)
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
    population = std::make_unique<double[]>(N.x * N.y * N.z * v->getQ() * 2);
    rho = std::make_unique<double[]>(N.x * N.y * N.z);
    u = std::make_unique<double[]>(N.x * N.y * N.z * v->getD());
    
    if(!quiet)
    {
        std::cout << "Lattice Boltzmann Simulation configured with a " << N.x << "x" << N.y << "x" << N.z << " lattice." << std::endl;
        std::cout << "The velocity has " << v->getQ() << " components." << std::endl;
        std::cout << "Kinematic viscosity was set to " << nu << "." << std::endl;
    }
}

void LBM::applyInitial(unsigned int x, unsigned int y, unsigned int z)
{
    InitialCondition(x,y,z, *this );
}

void LBM::applyBoundary(unsigned int x, unsigned int y, unsigned int z, VectorXd& f)
{
    for(std::function<void (unsigned int, unsigned int, unsigned int, Eigen::VectorXd &, LBM &)> bc : BoundaryConditions) bc(x,y,z, f, *this);
}

void LBM::setInitialCondition(std::function<void(unsigned int,unsigned int,unsigned int, LBM&)> initial_condition)
{
    this->InitialCondition = initial_condition;
}

void LBM::addBoundaryCondition(std::function<void(unsigned int,unsigned int,unsigned int, Eigen::VectorXd&, LBM&)> boundary)
{
    this->BoundaryConditions.push_back(boundary);
}


void LBM::saveToBin(unsigned int step)
{
    std::filesystem::create_directories("./bin_results");

    // Allocate arrays for the entire domain
    std::vector<double> u_x(N.x * N.y * N.z, 0.0);
    std::vector<double> u_y(N.x * N.y * N.z, 0.0);
    std::vector<double> u_z(N.x * N.y * N.z, 0.0);
    unsigned int idx = 0;

    for (unsigned int z = 0; z < N.z; ++z)
    {
        for (unsigned int y = 0; y < N.y; ++y)
        {
            for (unsigned int x = 0; x < N.x; ++x)
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
    std::ofstream file_u_x("./bin_results/u_x_" + std::to_string(step) + ".bin", std::ios::binary);
    file_u_x.write(reinterpret_cast<char*>(u_x.data()), u_x.size() * sizeof(double));
    file_u_x.close();

    // Save u_y to a binary file if applicable
    if (v->getD() > 1)
    {
        std::ofstream file_u_y("./bin_results/u_y_" + std::to_string(step) + ".bin", std::ios::binary);
        file_u_y.write(reinterpret_cast<char*>(u_y.data()), u_y.size() * sizeof(double));
        file_u_y.close();
    }

    // Save u_z to a binary file if applicable
    if (v->getD() == 3)
    {
        std::ofstream file_u_z("./bin_results/u_z_" + std::to_string(step) + ".bin", std::ios::binary);
        file_u_z.write(reinterpret_cast<char*>(u_z.data()), u_z.size() * sizeof(double));
        file_u_z.close();
    }
}