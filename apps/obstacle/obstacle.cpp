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

	for (auto& cl : msh.cells)
	{
		auto gr = make_hho_laplacian(msh, cl, hdi);
        Matrix<T, Dynamic, Dynamic> stab = make_hho_fancy_stabilization(msh, cl, gr.first, hdi);
        Matrix<T, Dynamic, Dynamic> lc = gr.second + stab;
        Matrix<T, Dynamic, 1> f = Matrix<T, Dynamic, 1>::Zero(lc.rows());
        f = make_rhs(msh, cl, hdi.cell_degree(), rhs_fun);
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




