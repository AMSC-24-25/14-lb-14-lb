#include <memory>
#include <functional>
#include <vector>
#include <exception>
#include <Eigen/Dense>
#include <iostream>
#include <LatticeBoltzmannMethod/Obstacle.hpp>

#ifndef __LBM_H
#define __LBM_H

namespace LatticeBoltzmannMethod{

    //a class for customized simulations using Lattice Boltzmann Method
    class LBM 
    {
        public:
            class VelocitySet
            {
                public: 
                    enum StandardSet { D1Q3, D2Q9, D3Q15, D3Q19, D3Q27 };
                
                protected: 
                    //dimension of the domain's vector space
                    const unsigned int D;
                    //directions considered for the mesoscopic simulation of the fluid
                    const unsigned int Q;
                    Eigen::MatrixXd c;
                    Eigen::VectorXd w;
                    std::array<int, 27> inverse;
                    std::vector<int> base;
                public:
                    VelocitySet(StandardSet set);
                    const unsigned int getD();
                    const unsigned int getQ();
                    const Eigen::MatrixXd& get_c();
                    const Eigen::VectorXd& get_w();
                    // Given the i-th direction of the velocity set, 
                    //the method returns the index of the direction opposed to it
                    const unsigned int directionIndexInvert(unsigned int i);
                    const unsigned int directionIndexBase(unsigned int i);
                private:
                    static const unsigned int fromStdD(StandardSet std);
                    static const unsigned int fromStdQ(StandardSet std);
            };

            std::function<void(unsigned int,unsigned int,unsigned int, LBM&)> InitialCondition;
            std::vector<std::function<void(unsigned int,unsigned int,unsigned int, Eigen::VectorXd&, LBM&)>> BoundaryConditions;
            std::vector<std::shared_ptr<Obstacle>> obstacles;
            
            struct dimensions {
                unsigned int x;
                unsigned int y;
                unsigned int z;
            } const N;
            struct x_location

            {
                long start;
                long end;
                bool left_pad;
                bool right_pad;
            } x_loc; 


            const double nu;
            const double tau = 3.0*nu+0.5;
            const double cs = 1.0 / 1.73205080756887729352;
            const double omgSym = 1/(nu/(cs*cs)+0.5);
            const double lambda = 0.25; //Provides the most stable simulations according to the book
            const double omgAnti = 1/(lambda*cs*cs/nu + 0.5);
            std::unique_ptr<VelocitySet> v;
        
        private:
            bool computeFlowProperties;
            bool quiet = false;
            unsigned int step = 0;
            int node_id;
            int x_len_no_pad;
            int x_len;

            std::unique_ptr<double[]> population;
            std::unique_ptr<double[]> rho;
            std::unique_ptr<double[]> u;

        public:
            void init_equilibrium();
            void stream_collide_save();

            LBM(int world_rank, int world_size, LBM::VelocitySet::StandardSet vSet, LBM::dimensions d,  double nu);

            void setInitialCondition(std::function<void(unsigned int,unsigned int,unsigned int, LBM&)> InitialCondition);
            void addBoundaryCondition(std::function<void(unsigned int,unsigned int,unsigned int, Eigen::VectorXd&, LBM&)>);
            void addObstacle(std::shared_ptr<Obstacle> o);

            inline const double get_rho(unsigned int x, unsigned int y, unsigned int z);
            inline void set_rho(unsigned int x, unsigned int y, unsigned int z, double rho);
            inline const Eigen::VectorXd get_u(unsigned int x, unsigned int y, unsigned int z);
            inline void set_u(unsigned int x, unsigned int y, unsigned int z, const Eigen::VectorXd& u);
            inline Eigen::VectorXd getPopulation(unsigned int x, unsigned int y, unsigned int z);

            inline void setVerbose();
            inline void setQuiet();
            
            void saveToBin(unsigned int step);


        
        private:
        
            inline Eigen::VectorXd populationAdjacent(unsigned int x, unsigned int y, unsigned int z);
            inline void savePopulation(unsigned int x, unsigned int y, unsigned int z, const Eigen::VectorXd& population);
            inline void savePopulationInit(unsigned int x, unsigned int y, unsigned int z, const Eigen::VectorXd& population);

            void applyInitial(unsigned int x, unsigned int y, unsigned int z);
            void applyBoundaryAndObstacle(unsigned int x, unsigned int y, unsigned int z, Eigen::VectorXd& f);

            inline void check_coordinates(unsigned int x, unsigned int y, unsigned int z);
            inline unsigned int index_r(unsigned int x, unsigned int y, unsigned int z);
            inline unsigned int index_u(unsigned int x, unsigned int y, unsigned int z);
            inline unsigned int index_f(unsigned int x, unsigned int y, unsigned int z);
    };



    const double LBM::get_rho(unsigned int x, unsigned int y, unsigned int z) { return rho[index_r(x,y,z)]; }
    void LBM::set_rho(unsigned int x, unsigned int y, unsigned int z, double rho) { this->rho[index_r(x,y,z)] = rho; }

    const Eigen::VectorXd LBM::get_u(unsigned int x, unsigned int y, unsigned int z)
    {
        Eigen::VectorXd u(v->getD());
        for(unsigned int i = 0; i < v->getD(); ++i)
        {
            u(i) = this->u[index_u(x,y,z) + i];
        } 

        return u;
    }

    void LBM::set_u(unsigned int x, unsigned int y, unsigned int z, const Eigen::VectorXd& u)
    {
        for(unsigned int i = 0; i < v->getD(); ++i)
        {
            this->u[index_u(x,y,z) + i] = u(i);
        } 
    }

    Eigen::VectorXd LBM::getPopulation(unsigned int x, unsigned int y, unsigned int z)
    {
        Eigen::VectorXd p( v->getQ() );
        //const Eigen::MatrixXd& c = v->get_c();
        for(unsigned int i = 0; i < v->getQ(); ++i)
        {
            p(i) = this->population[index_f(x, y, z) + x_len*N.y*N.z * (step & 1) * v->getQ() + i];
        } 
        return p;
    }


    Eigen::VectorXd LBM::populationAdjacent(unsigned int x, unsigned int y, unsigned int z) //PLEASE CHECK THIS FUNCTION
    {
        Eigen::VectorXd adj( v->getQ() );
        const Eigen::MatrixXd& c = v->get_c();
        for(unsigned int i = 0; i < v->getQ(); ++i)
        {
            int x_adj = x - c(i, 0);
            int y_adj = y - ((v->getD() > 1) ? c(i, 1) : 0);
            int z_adj = z - ((v->getD() == 3) ? c(i, 2) : 0);
            if (x_adj >= 0 && x_adj < (int)N.x &&
                y_adj >= 0 && y_adj < (int)N.y &&
                z_adj >= 0 && z_adj < (int)N.z)
            {
                adj(i) = this->population[index_f(x_adj, y_adj, z_adj) + x_len*N.y*N.z * (step & 1) * v->getQ() + i];
            }
            else
            {
                adj(i) = 0.0;
            }
        }  

            return adj;
        }

    void LBM::savePopulation(unsigned int x, unsigned int y, unsigned int z, const Eigen::VectorXd& population)
    {
        for(unsigned int i = 0; i < v->getQ(); ++i)
        {
            this->population[index_f(x,y,z) + x_len*N.y*N.z * (~step & 1) * v->getQ() + i] = population(i);
        }     
    }

    //ON THE ZERO STEP, THE POPULATION IS SAVED TO THE SECOND HALF OF THE ARRAY 
    void LBM::savePopulationInit(unsigned int x, unsigned int y, unsigned int z, const Eigen::VectorXd& population)
    {
        for(unsigned int i = 0; i < v->getQ(); ++i)
        {
            this->population[index_f(x,y,z) + x_len*N.y*N.z * (step & 1) * v->getQ() + i] = population(i);
        }     
    }

    unsigned int LBM::index_r(unsigned int x, unsigned int y, unsigned int z) { 
        check_coordinates(x,y,z);
        x = x - x_loc.start + x_loc.left_pad;
        return (N.y * x + y) * N.z + z;
    }

    unsigned int LBM::index_u(unsigned int x, unsigned int y, unsigned int z) { 
        check_coordinates(x,y,z);
        return index_r(x,y,z) * v->getD(); 
    }

    unsigned int LBM::index_f(unsigned int x, unsigned int y, unsigned int z) { 
        check_coordinates(x,y,z);
        return index_r(x,y,z) * v->getQ(); 
    }

    void LBM::check_coordinates(unsigned int x, unsigned int y, unsigned int z)
    {
        if(x >= N.x || x < 0)
        {
            throw std::out_of_range(std::string{"Invalid x coordinate\n"});
        }
        if(y < 0 || y >= N.y)
        {
            throw std::out_of_range(std::string{"Invalid y coordinate\n"});
        }
        if(z < 0 || z >= N.z)
        {
            throw std::out_of_range(std::string{"Invalid z coordinate\n"});
        }
    }


    inline void LBM::setVerbose() { quiet = false; }
    inline void LBM::setQuiet() { quiet = true; }

}
#endif /* __LBM_H */