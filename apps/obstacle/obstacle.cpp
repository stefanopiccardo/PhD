/*
 *       /\        Matteo Cicuttin (C) 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    This is a prototype implementation of the Obstacle
 *  /_\/_\/_\/_\   problem using HHO.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Implementation of Discontinuous Skeletal methods on arbitrary-dimensional,
 * polytopal meshes using generic programming.
 * M. Cicuttin, D. A. Di Pietro, A. Ern.
 * Journal of Computational and Applied Mathematics.
 * DOI: 10.1016/j.cam.2017.09.017
 */

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/SparseExtra>

using namespace Eigen;

#include "core/core"
#include "core/solvers"
#include "dataio/silo_io.hpp"

#include "methods/hho"



template<typename Mesh>
void run_hho_obstacle(const Mesh& msh, size_t degree)
{

	hho_degree_info hdi(0, degree);

	using T = typename Mesh::coordinate_type;

	auto rhs_fun = [](const typename Mesh::point_type& pt) -> T {
        return 2.0 * M_PI * M_PI * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    auto sol_fun = [](const typename Mesh::point_type& pt) -> T {
        return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };


    auto bcs_fun = [&](const typename Mesh::point_type& pt) -> T {
        return sol_fun(pt);
    };

    auto obstacle_fun = [&](const typename Mesh::point_type& pt) -> T {
        return 0.5;
    };


    Matrix<T, Dynamic, 1> alpha = Matrix<T, Dynamic, 1>::Zero( msh.cells.size() );
	Matrix<T, Dynamic, 1> beta  = Matrix<T, Dynamic, 1>::Ones( msh.cells.size() );
	Matrix<T, Dynamic, 1> gamma = Matrix<T, Dynamic, 1>::Zero( msh.cells.size() );
	T c = 1.0;

	size_t iter = 0;
    bool has_to_iterate = true;
    while ( has_to_iterate && iter < 50 )
    {
    	std::cout << bold << "Iteration " << iter << reset << std::endl;
    	std::stringstream ss;
    	ss << "obstacle_cycle_" << iter << ".silo";

    	silo_database silo;
    	silo.create( ss.str() );
    	silo.add_mesh(msh, "mesh");

    	for (size_t i = 0; i < msh.cells.size(); i++)
    		gamma(i) = 0.5;

    	Matrix<T, Dynamic, 1> diff = beta + c * ( alpha - gamma );
    	std::vector<bool> in_A;
    	in_A.resize(diff.size());
    	Matrix<T, Dynamic, 1> active = Matrix<T, Dynamic, 1>::Zero(diff.size());

    	for (size_t i = 0; i < diff.size(); i++)
    	{
    		in_A.at(i) = (diff(i) >= 0);
    		active(i) = (diff(i) >= 0);
    	}

    	timecounter tc;
    	tc.tic();

    	auto assembler = make_obstacle_assembler(msh, in_A, hdi);

    	//assembler.dump_tables();

    	for (auto& cl : msh.cells)
		{
			auto gr = make_hho_laplacian(msh, cl, hdi);
        	Matrix<T, Dynamic, Dynamic> stab = make_hho_fancy_stabilization(msh, cl, gr.first, hdi);
        	Matrix<T, Dynamic, Dynamic> lc = gr.second + stab;
        	Matrix<T, Dynamic, 1> f = Matrix<T, Dynamic, 1>::Zero(lc.rows());
        	f = make_rhs(msh, cl, hdi.cell_degree(), rhs_fun);
        	assembler.assemble_A(msh, cl, lc, f, gamma, sol_fun);
		}

		assembler.finalize();

		tc.toc();
		std::cout << bold << yellow << "Assembly: " << tc << " seconds" << reset << std::endl;


		ss.str("");
		ss << "matrix_" << iter << ".txt";

		dump_sparse_matrix(assembler.LHS, ss.str());


		silo.add_variable("mesh", "difference", diff.data(), diff.size(), zonal_variable_t);
		silo.add_variable("mesh", "active", active.data(), active.size(), zonal_variable_t);


		tc.tic();

		SparseLU<SparseMatrix<T>>  solver;
    	solver.analyzePattern(assembler.LHS);
    	solver.factorize(assembler.LHS);

    	Matrix<T, Dynamic, 1> sol;
    	sol = solver.solve(assembler.RHS);
    	tc.toc();
		std::cout << bold << yellow << "Solver: " << tc << " seconds" << reset << std::endl;

		Matrix<T, Dynamic, Dynamic> alpha_prev = alpha;
		assembler.expand_solution(msh, sol, sol_fun, gamma, alpha, beta);


		silo.add_variable("mesh", "alpha", alpha.data(), alpha.size(), zonal_variable_t);
		silo.add_variable("mesh", "beta", beta.data(), beta.size(), zonal_variable_t);

		silo.close();

		if ( (alpha_prev - alpha).norm() < 1e-7 )
			break;

		iter++;
    }

    
   




	
	

}

int main(int argc, char **argv)
{
    using T = double;
    size_t degree = 0;

    mesh_init_params<T> mip;
    mip.Nx = 5;
    mip.Ny = 5;

    bool displace_nodes = true;

    int ch;
    while ( (ch = getopt(argc, argv, "k:N:")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                if (degree != 0 || degree != 1)
                {
                	std::cout << "Degree can be 0 or 1. Falling back to 1" << std::endl;
                	degree = 1;
                }
                break;

            case 'N':
            	mip.Nx = atoi(optarg);
                mip.Ny = atoi(optarg);
                break;

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;


    timecounter tc;

    /************** BUILD MESH **************/
    tc.tic();
    quad_mesh<T> msh(mip);
    tc.toc();
    std::cout << bold << yellow << "Mesh generation: " << tc << " seconds" << reset << std::endl;


    run_hho_obstacle(msh, degree);

    return 0;
}




