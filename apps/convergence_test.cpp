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
#include <random>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/SparseExtra>



using namespace Eigen;

#include "core/core"
#include "core/solvers"
#include "methods/hho"

#include "sol2/sol.hpp"

#include "dataio/silo_io.hpp"

enum class method_type
{
    MIXED_ORDER_INF,
    EQUAL_ORDER,
    MIXED_ORDER_SUP
};

struct convergence_test_params
{
    size_t      deg_min;
    size_t      deg_max;
    size_t      min_N;
    size_t      steps;
    bool        preconditioner;
    bool        direct;
    bool        stab_hho;
    method_type mt;

    convergence_test_params() :
        deg_min(0),
        deg_max(6),
        min_N(4),
        steps(5),
        preconditioner(true),
        direct(false),
        stab_hho(true),
        mt( method_type::EQUAL_ORDER )
    {}
};

int test_method_convergence(const convergence_test_params& ctp)
{
    using RealType = double;


    std::random_device r;
    std::default_random_engine e1(r());


    size_t deg_min = ctp.deg_min;
    size_t deg_max = ctp.deg_max;

    size_t min_elems = ctp.min_N;
    size_t steps = ctp.steps;

    bool preconditioner = ctp.preconditioner;
    bool direct = ctp.direct;

    auto rhs_fun = [](const typename quad_mesh<RealType>::point_type& pt) -> RealType {
        return 2.0 * M_PI * M_PI * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    auto sol_fun = [](const typename quad_mesh<RealType>::point_type& pt) -> RealType {
        return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    auto sol_grad = [](const typename quad_mesh<RealType>::point_type& pt) -> Matrix<RealType, 1, 2> {
        Matrix<RealType, 1, 2> ret;
        ret(0) = M_PI*std::cos(M_PI*pt.x())*std::sin(M_PI*pt.y());
        ret(1) = M_PI*std::sin(M_PI*pt.x())*std::cos(M_PI*pt.y());
        return ret;
    };

    for (size_t k = deg_min; k < deg_max+1; k++)
    {
        std::cout << "Testing degree " << k << std::endl;

        std::vector<RealType> errors_int, errors_mm, errors_energy;
        errors_int.resize(steps);
        errors_mm.resize(steps);
        errors_energy.resize(steps);

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

        hho_degree_info hdi(k+1, k);

        for (size_t i = 0, N = min_elems; i < steps; i++, N *= 2)
        {
            mesh_init_params<RealType> mip;
            mip.Nx = N;
            mip.Ny = N;

            quad_mesh<RealType> msh(mip);

            std::stringstream ss;
            ss << "convergence_test_N_" << N << "_k_" << k << ".silo";

            RealType delta = 0.1/N;
            std::uniform_real_distribution<RealType> uniform_dist(-delta, delta);
            
            for (size_t i = 0; i < msh.points.size(); i++)
            {
                RealType dx = uniform_dist(e1);
                RealType dy = uniform_dist(e1);
                auto pt = msh.points[i];
                auto rx = pt.x() + dx;
                auto ry = pt.y() + dy;
                msh.points[i] = typename quad_mesh<RealType>::point_type({rx, ry});
            }
            

            silo_database silo;
            silo.create(ss.str());
            silo.add_mesh(msh, "mesh");

            for (auto& fc : msh.faces)
            {
                if (fc.is_boundary)
                    fc.bndtype = boundary::DIRICHLET;
            }

            auto assembler = make_assembler(msh, hdi);
            for (auto& cl : msh.cells)
            {
                auto gr = make_hho_laplacian(msh, cl, hdi);

                Matrix<RealType, Dynamic, Dynamic> stab;
                if ( ctp.stab_hho )
                    stab = make_hho_fancy_stabilization(msh, cl, gr.first, hdi);
                else
                    stab = make_hho_naive_stabilization(msh, cl, hdi);

                Matrix<RealType, Dynamic, Dynamic> lc = gr.second + stab;
                Matrix<RealType, Dynamic, 1> f = make_rhs(msh, cl, hdi.cell_degree(), rhs_fun);
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
                cgp.max_iter = 3*assembler.LHS.rows();
                cgp.histfile = ss.str();
                auto exit_reason = conjugated_gradient(assembler.LHS, assembler.RHS, sol, cgp);

                if ( exit_reason != cg_exit_reason::CONVERGED )
                    std::cout << "Warning! Solver didn't converge..." << std::endl;
            }

            size_t cell_i = 0;
            for (auto& cl : msh.cells)
            {
                cell_basis<quad_mesh<RealType>, RealType> cb(msh, cl, hdi.cell_degree());
                auto cbs = cb.size();

                Matrix<RealType, Dynamic, 1> cdofs = sol.block(cell_i*cbs, 0, cbs, 1);

                auto qps = integrate(msh, cl, 2*hdi.cell_degree());
                for (auto& qp : qps)
                {
                    Matrix<RealType, Dynamic, 1> phi = cb.eval_basis(qp.first);
                    auto val = cdofs.dot( phi );
                    auto real_val = sol_fun( qp.first );
                    errors_int.at(i) += qp.second*(real_val - val)*(real_val - val);

                    Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, hdi.cell_degree());
                    Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cl, hdi.cell_degree(), sol_fun);
                    Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
                    Matrix<RealType, Dynamic, 1> diff = real_dofs - cdofs;
                    errors_mm.at(i) += diff.dot(mass*diff);
                }

                cell_basis<quad_mesh<RealType>, RealType> rb(msh, cl, hdi.reconstruction_degree());
                auto rbs = rb.size();

                auto qps2 = integrate(msh, cl, 2*hdi.reconstruction_degree());
                for (auto& qp : qps)
                {
                    auto alldofs = assembler.take_local_data(msh, cl, sol, sol_fun);
                    auto gr = make_hho_laplacian(msh, cl, hdi);

                    Matrix<RealType, Dynamic, 1> recdofs = gr.first * alldofs;
                    Matrix<RealType, Dynamic, 2> dphi = rb.eval_gradients(qp.first);
                    Matrix<RealType, 1, 2> gval = Matrix<RealType, 1, 2>::Zero();

                    for (size_t i = 1; i < dphi.rows(); i++)
                        gval = gval + recdofs(i-1)*dphi.block(i, 0, 1, 2);

                    auto real_gval = sol_grad( qp.first );
                    errors_energy.at(i) += qp.second*(real_gval - gval).dot(real_gval - gval);
                }

                cell_i++;
            }

            auto mesh_h = diameter(msh, msh.cells.front());
            hho_conv_ofs << mesh_h << " " << errors_int.at(i) << " " << errors_mm.at(i) << std::endl;

            if (i > 0)
            {
                auto error_prev_int = sqrt(errors_int.at(i-1));
                auto error_cur_int = sqrt(errors_int.at(i));
                std::cout << log10(error_prev_int/error_cur_int) / log10(2) << "\t\t";

                auto error_prev_mm = sqrt(errors_mm.at(i-1));
                auto error_cur_mm = sqrt(errors_mm.at(i));
                std::cout << log10(error_prev_mm/error_cur_mm) / log10(2) << "\t\t";

                auto error_prev_energy = sqrt(errors_energy.at(i-1));
                auto error_cur_energy = sqrt(errors_energy.at(i));
                std::cout << log10(error_prev_energy/error_cur_energy) / log10(2) << std::endl;
            }
        }

        hho_conv_ofs.close();
    }

    return 0;
}

int main(int argc, char **argv)
{
    sol::state lua;
    lua.open_libraries(sol::lib::math, sol::lib::base);

    if (argc < 2)
    {
        convergence_test_params ctp;
        test_method_convergence(ctp);
        return 0;
    }

    auto r = lua.do_file(argv[1]);
    if ( !r.valid() )
    {
        std::cout << "Problems opening configuration file" << std::endl;
        return 1;
    }

    convergence_test_params ctp;

    auto deg_min    = lua["deg_min"];   if (deg_min.valid())    ctp.deg_min = deg_min;
    auto deg_max    = lua["deg_max"];   if (deg_max.valid())    ctp.deg_max = deg_max;
    auto min_N      = lua["min_N"];     if (min_N.valid())      ctp.min_N = min_N;
    auto steps      = lua["steps"];     if (steps.valid())      ctp.steps = steps;
    auto precond    = lua["precond"];   if (precond.valid())    ctp.preconditioner = precond;
    auto direct     = lua["direct"];    if (direct.valid())     ctp.direct = direct;
    auto stab_hho   = lua["stab_hho"];  if (stab_hho.valid())   ctp.stab_hho = stab_hho;

    test_method_convergence(ctp);
    return 0;
}


/*
size_t      deg_min;
size_t      deg_max;
size_t      min_N;
size_t      steps;
bool        preconditioner;
bool        direct;
bool        stab_hho;
method_type mt;
 */
