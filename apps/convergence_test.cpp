/*
 *       /\        Matteo Cicuttin (C) 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    This is a prototype implementation of the CutHHO method,
 *  /_\/_\/_\/_\   an unfitted Hybrid High-order method.
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
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <sstream>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/SparseExtra>



using namespace Eigen;

#include "core/core"
#include "core/solvers"
#include "methods/hho"

struct convergence_test_params
{
    size_t  deg_min;
    size_t  deg_max;
    size_t  min_N;
    size_t  steps;
    bool    preconditioner;
    bool    direct;

    convergence_test_params() :
        deg_min(0),
        deg_max(6),
        min_N(4),
        steps(5),
        preconditioner(true),
        direct(false)
    {}
};

int test_method_convergence(const convergence_test_params& ctp)
{
    using RealType = double;

    size_t deg_min = ctp.deg_min;
    size_t deg_max = ctp.deg_max;

    size_t min_elems = ctp.min_N;
    size_t steps = ctp.steps;

    bool preconditioner = ctp.preconditioner;
    bool direct = ctp.direct;

    auto rhs_fun = [](const typename mesh<RealType>::point_type& pt) -> RealType {
        return 2.0 * M_PI * M_PI * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    auto sol_fun = [](const typename mesh<RealType>::point_type& pt) -> RealType {
        return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    for (size_t k = deg_min; k < deg_max+1; k++)
    {
        std::cout << "Testing degree " << k << ". Expected rate " << k+1 << std::endl;

        std::vector<RealType> errors_int, errors_mm;
        errors_int.resize(steps);
        errors_mm.resize(steps);

        for (size_t i = 0; i < steps; i++)
        {
            errors_int.at(i) = 0.0;
            errors_mm.at(i) = 0.0;
        }

        std::stringstream ss;
        if (preconditioner)
            ss << "hho_history_precond_" << k << ".txt";
        else
            ss << "hho_history_" << k << ".txt";

        std::ofstream hho_conv_ofs(ss.str());

        for (size_t i = 0, N = min_elems; i < steps; i++, N *= 2)
        {
            mesh_init_params<RealType> mip;
            mip.Nx = N;
            mip.Ny = N;

            mesh<RealType> msh(mip);

            for (auto& fc : msh.faces)
            {
                if (fc.is_boundary)
                    fc.bndtype = boundary::DIRICHLET;
            }

            auto assembler = make_assembler(msh, k);
            for (auto& cl : msh.cells)
            {
                auto gr = make_hho_laplacian(msh, cl, k);
                Matrix<RealType, Dynamic, Dynamic> stab = make_hho_naive_stabilization(msh, cl, k);
                Matrix<RealType, Dynamic, Dynamic> lc = gr.second + stab;
                Matrix<RealType, Dynamic, 1> f = make_rhs(msh, cl, k, rhs_fun);
                assembler.assemble(msh, cl, lc, f, sol_fun);
            }

            assembler.finalize();

            Matrix<RealType, Dynamic, 1> sol;

            if (direct)
            {
                SparseLU<SparseMatrix<RealType>>  solver;

                solver.analyzePattern(assembler.LHS);
                solver.factorize(assembler.LHS);
                sol = solver.solve(assembler.RHS);
            }
            else
            {
                std::stringstream ss;
                if (preconditioner)
                    ss << "cg_history_precond_" << N << "_" << k << ".txt";
                else
                    ss << "cg_history_" << N << "_" << k << ".txt";

                cg_params<RealType> cgp;
                cgp.apply_preconditioner = preconditioner;
                cgp.convergence_threshold = 1e-12;
                cgp.max_iter = assembler.LHS.rows();
                cgp.histfile = ss.str();
                auto exit_reason = conjugated_gradient(assembler.LHS, assembler.RHS, sol, cgp);

                if ( exit_reason != cg_exit_reason::CONVERGED )
                    std::cout << "Warning! Solver didn't converge..." << std::endl;
            }

            size_t cell_i = 0;
            for (auto& cl : msh.cells)
            {
                cell_basis<mesh<RealType>, RealType> cb(msh, cl, k);
                auto cbs = cb.size();

                Matrix<RealType, Dynamic, 1> cdofs = sol.block(cell_i*cbs, 0, cbs, 1);

                auto qps = integrate(msh, cl, 2*k);
                for (auto& qp : qps)
                {
                    Matrix<RealType, Dynamic, 1> phi = cb.eval_basis(qp.first);
                    auto val = cdofs.dot( phi );
                    auto real_val = sol_fun( qp.first );
                    errors_int.at(i) += qp.second*(real_val - val)*(real_val - val);

                    Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, k);
                    Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cl, k, sol_fun);
                    Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                    Matrix<RealType, Dynamic, 1> diff = real_dofs - cdofs;
                    errors_mm.at(i) += diff.dot(mass*diff);
                }

                cell_i++;
            }

            hho_conv_ofs << errors_int.at(i) << " " << errors_mm.at(i) << std::endl;

            if (i > 0)
            {
                auto error_prev_int = sqrt(errors_int.at(i-1));
                auto error_cur_int = sqrt(errors_int.at(i));
                std::cout << log10(error_prev_int/error_cur_int) / log10(2) << "\t\t";

                auto error_prev_mm = sqrt(errors_mm.at(i-1));
                auto error_cur_mm = sqrt(errors_mm.at(i));
                std::cout << log10(error_prev_mm/error_cur_mm) / log10(2) << std::endl;
            }
        }

        hho_conv_ofs.close();
    }

    return 0;
}

int main(int argc, char **argv)
{
    convergence_test_params ctp;
    test_method_convergence(ctp);
    return 0;
}
