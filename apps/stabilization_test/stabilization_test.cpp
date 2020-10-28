/*
 *       /\        Matteo Cicuttin (C) 2017,2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    This is ProtoN, a library for fast Prototyping of
 *  /_\/_\/_\/_\   Numerical methods.
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
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/SparseExtra>

using namespace Eigen;

#include "core/core"
#include "core/solvers"
#include "methods/hho"


double test_stabilization(size_t N, size_t k)
{
    using RealType = double;

    hho_degree_info hdi(k, k);


    mesh_init_params<RealType> mip;
    mip.Nx = N;
    mip.Ny = N;

    quad_mesh<RealType> msh(mip);


    auto rhs_fun = [](const typename quad_mesh<RealType>::point_type& pt) -> RealType {
        return 2.0 * M_PI * M_PI * std::sin(2*M_PI*pt.x()) * std::sin(2*M_PI*pt.y());
        //return pt.x();
    };


    RealType error = 0.0;

    for (auto& cl : msh.cells)
    {
        auto gr = make_hho_laplacian(msh, cl, hdi);

        Matrix<RealType, Dynamic, Dynamic> stab;   
        stab = make_hho_fancy_stabilization(msh, cl, gr.first, hdi);

        Matrix<RealType, Dynamic, 1> proj = project_function(msh, cl, hdi, rhs_fun);

        error += proj.dot(stab*proj);

        break;
    }

    return std::sqrt(error);
}

int main(void)
{

    for (size_t k = 0; k < 6; k++)
    {
        std::vector<double> errors;

        for (size_t N = 2; N < 64; N *= 2)
        {
            auto err = test_stabilization(N,k);
            errors.push_back(err);
        }

        for (size_t i = 1; i < errors.size(); i++)
            std::cout << std::setprecision(2) << std::log(errors.at(i-1)/errors.at(i))/std::log(2.0) << "  ";
        std::cout << std::endl;

    }

}





