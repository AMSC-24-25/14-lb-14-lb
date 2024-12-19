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
#include <functional>
#include <vector>
#include <Eigen/Dense>

#ifndef __LBM_H
#define __LBM_H

class LBM 
{
    public:
        class VelocitySet
        {
            public: 
                enum StandardSet { D1Q3, D2Q9, D3Q15, D3Q19, D3Q27 };
            
            protected: 
                const int D;
                const int Q;
                Eigen::MatrixXd c;
                Eigen::VectorXd w;

            public:
                VelocitySet(StandardSet set);
                const int getD();
                const int getQ();
                const Eigen::MatrixXd& get_c();
                const Eigen::VectorXd& get_w();

            private:
                static const int fromStdD(StandardSet std);
                static const int fromStdQ(StandardSet std);
        };

    std::function<void(int,int,int, LBM&)> InitialCondition;
    std::vector<std::function<void(int,int,int, Eigen::VectorXd&, LBM&)>> BoundaryConditions;
        struct dimensions {
            int x;
            int y;
            int z;
        } const N;
        
        const double nu;
        const double tau = 3.0*nu+0.5;
        const double cs = 1.0/1.732;
    
    private:
        std::unique_ptr<VelocitySet> v;
        
        

        bool computeFlowProperties;
        bool quiet;
        int step = 0;

        std::unique_ptr<double[]> population;
        std::unique_ptr<double[]> rho;
        std::unique_ptr<double[]> u;

    public:
        void init_equilibrium();
        void stream_collide_save();

        LBM::LBM(LBM::VelocitySet::StandardSet vSet, LBM::dimensions d,  double nu);

        inline const double get_rho(int x, int y, int z);
        inline void set_rho(int x, int y, int z, double rho);
        inline const Eigen::VectorXd get_u(int x, int y, int z);
        inline void set_u(int x, int y, int z, const Eigen::VectorXd& u);
        inline Eigen::VectorXd getPopulation(int x, int y, int z);
    
    private:
    
        inline Eigen::VectorXd populationAdjacent(int x, int y, int z);
        inline void savePopulation(int x, int y, int z, const Eigen::VectorXd& population);

        void applyInitial(int x, int y, int z);
        void applyBoundary(int x, int y, int z, Eigen::VectorXd& f);

        inline int index_r(int x, int y, int z);
        inline int index_u(int x, int y, int z);
        inline int index_f(int x, int y, int z);



};

#endif /* __LBM_H */

const double LBM::get_rho(int x, int y, int z) { return rho[index_r(x,y,z)]; }
void LBM::set_rho(int x, int y, int z, double rho) { this->rho[index_r(x,y,z)] = rho; }

const Eigen::VectorXd LBM::get_u(int x, int y, int z)
{
    Eigen::VectorXd u(v->getD());
    for(int i = 0; i < v->getD(); ++i)
    {
        u(i) = this->u[index_u(x,y,z) + i];
    } 

    return u;
}

void LBM::set_u(int x, int y, int z, const Eigen::VectorXd& u)
{
    for(int i = 0; i < v->getD(); ++i)
    {
        this->u[index_u(x,y,z) + i] = u(i);
    } 
}

Eigen::VectorXd LBM::getPopulation(int x, int y, int z)
{
    Eigen::VectorXd p( v->getQ() );
    const Eigen::MatrixXd& c = v->get_c();
    for(int i = 0; i < v->getQ(); ++i)
    {
        p(i) = this->population[index_f(x, y, z) + N.x*N.y*N.z * (~step & 1) + i];
    } 
    return p;
}


Eigen::VectorXd LBM::populationAdjacent(int x, int y, int z)
{
    Eigen::VectorXd adj( v->getQ() );
    const Eigen::MatrixXd& c = v->get_c();
    for(int i = 0; i < v->getQ(); ++i)
    {
        int x_adj = x + c(i, 0);
        int y_adj = y + (v->getD() > 1)? c(i, 1) : 0;
        int z_adj = z + (v->getD() == 3)? c(i, 2) : 0;
        adj(i) = this->population[index_f(x_adj, y_adj, z_adj) + N.x*N.y*N.z * (~step & 1) + i];
    } 
    return adj;
}

void LBM::savePopulation(int x, int y, int z, const Eigen::VectorXd& population)
{
    for(int i = 0; i < v->getQ(); ++i)
    {
        this->population[index_f(x,y,z) + N.x*N.y*N.z * (step & 1) + i] = population(i);
    }     
}

int LBM::index_r(int x, int y, int z) { return (N.y * z + y) * N.x + x; }
int LBM::index_u(int x, int y, int z) { return index_r(x,y,z) * v->getD(); }
int LBM::index_f(int x, int y, int z) { return index_r(x,y,z) * v->getQ(); }
