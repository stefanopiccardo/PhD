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
#include <numeric>
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
#include "dataio/silo_io.hpp"

#include "methods/hho"
#include "methods/cuthho"






template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_operator(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
{
    return make_hho_laplacian(msh, cl, degree);
}



template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
{
    return make_hho_naive_stabilization(msh, cl, degree);
}

/*****************************************************************************
 *   Test stuff
 *****************************************************************************/
template<typename Mesh>
void
plot_basis_functions(const Mesh& msh)
{
    using T = typename Mesh::coordinate_type;

    std::ofstream c_ofs("cell_basis_check.dat");

    for (auto cl : msh.cells)
    {
        cell_basis<Mesh, T> cb(msh, cl, 3);

        auto tps = make_test_points(msh, cl);

        for (auto& tp : tps)
        {
            c_ofs << tp.x() << " " << tp.y() << " ";

            auto vals = cb.eval_basis(tp);
            for(size_t i = 0; i < cb.size(); i++)
                c_ofs << vals(i) << " ";

            c_ofs << std::endl;
        }
    }

    c_ofs.close();

    std::ofstream f_ofs("face_basis_check.dat");

    for (auto fc : msh.faces)
    {
        face_basis<Mesh, T> fb(msh, fc, 2);

        auto tps = make_test_points(msh, fc);

        for (auto& tp : tps)
        {
            f_ofs << tp.x() << " " << tp.y() << " ";

            auto vals = fb.eval_basis(tp);
            for(size_t i = 0; i < fb.size(); i++)
                f_ofs << vals(i) << " ";

            f_ofs << std::endl;
        }
    }

    f_ofs.close();
}

template<typename Mesh>
void
plot_quadrature_points(const Mesh& msh, size_t degree)
{
    std::ofstream c_ofs("cell_quadrature_check.dat");

    for (auto& cl : msh.cells)
    {
        auto qps = integrate(msh, cl, degree);

        for (auto& qp : qps)
        {
            c_ofs << qp.first.x() << " " << qp.first.y();
            c_ofs << " " << qp.second << std::endl;
        }
    }

    c_ofs.close();

    std::ofstream f_ofs("face_quadrature_check.dat");

    for (auto& fc : msh.faces)
    {
        auto qps = integrate(msh, fc, degree);

        for (auto& qp : qps)
        {
            f_ofs << qp.first.x() << " " << qp.first.y();
            f_ofs << " " << qp.second << std::endl;
        }
    }

    f_ofs.close();
}


template<typename Mesh>
void
test_mass_matrices(const Mesh& msh, size_t degree)
{
    using T = typename Mesh::coordinate_type;

    auto rhs_fun = [](const typename Mesh::point_type& pt) -> T {
        return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    std::ofstream c_ofs("cell_mass_check.dat");

    cell_basis<Mesh, T>::print_structure(degree);

    for (auto& cl : msh.cells)
    {
        Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, degree);
        Matrix<T, Dynamic, 1> rhs = make_rhs(msh, cl, degree, rhs_fun);
        Matrix<T, Dynamic, 1> sol = mass.llt().solve(rhs);

        cell_basis<T,T> cb(msh, cl, degree);

        auto tps = make_test_points(msh, cl);
        for (auto& tp : tps)
        {
            auto phi = cb.eval_basis(tp);
            auto val = sol.dot(phi);
            c_ofs << tp.x() << " " << tp.y() << " " << val << std::endl;
        }

    }

    c_ofs.close();


    std::ofstream f_ofs("face_mass_check.dat");

    for (auto& fc : msh.faces)
    {
        Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, degree);
        Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, degree, rhs_fun);
        Matrix<T, Dynamic, 1> sol = mass.llt().solve(rhs);

        face_basis<T,T> fb(msh, fc, degree);

        auto tps = make_test_points(msh, fc);
        for (auto& tp : tps)
        {
            auto phi = fb.eval_basis(tp);
            auto val = sol.dot(phi);
            f_ofs << tp.x() << " " << tp.y() << " " << val << std::endl;
        }

    }

    f_ofs.close();
}


















template<typename T>
struct temp_tri
{
    std::array< point<T,2>, 3 > pts;
};

template<typename T>
std::ostream&
operator<<(std::ostream& os, const temp_tri<T>& t)
{
    os << "line([" << t.pts[0].x() << "," << t.pts[1].x() << "],[" << t.pts[0].y() << "," << t.pts[1].y() << "]);" << std::endl;
    os << "line([" << t.pts[1].x() << "," << t.pts[2].x() << "],[" << t.pts[1].y() << "," << t.pts[2].y() << "]);" << std::endl;
    os << "line([" << t.pts[2].x() << "," << t.pts[0].x() << "],[" << t.pts[2].y() << "," << t.pts[0].y() << "]);" << std::endl;
    return os;
}

template<typename T>
std::vector<temp_tri<T>>
triangulate(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::cell_type& cl,
            element_location where)
{
    assert( is_cut(msh, cl) );

    auto tp = collect_triangulation_points(msh, cl, where);
    auto bar = barycenter(tp);

    std::vector<temp_tri<T>> tris;

    for (size_t i = 0; i < tp.size(); i++)
    {
        temp_tri<T> t;
        t.pts[0] = bar;
        t.pts[1] = tp[i];
        t.pts[2] = tp[(i+1)%tp.size()];

        tris.push_back(t);
    }

    return tris;
}

template<typename T>
void test_triangulation(const cuthho_mesh<T>& msh)
{
    std::ofstream ofs("triangulation_dump.m");
    for (auto& cl : msh.cells)
    {
        if ( !is_cut(msh, cl) )
            continue;

        auto tris = triangulate(msh, cl, element_location::IN_NEGATIVE_SIDE);

        for (auto& tri : tris)
            ofs << tri << std::endl;
    }

    ofs.close();
}

template<typename T>
std::vector< std::pair<point<T,2>, T> >
integrate(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::cell_type& cl,
          size_t degree, const element_location& where)
{
    if ( !is_cut(msh, cl) ) /* Element is not cut, use std. integration */
        return integrate(msh, cl, degree);

    std::vector< std::pair<point<T,2>, T> > ret;
    auto tris = triangulate(msh, cl, where);
    for (auto& tri : tris)
    {
        auto qpts = triangle_quadrature(tri.pts[0], tri.pts[1], tri.pts[2], degree);
        ret.insert(ret.end(), qpts.begin(), qpts.end());
    }

    return ret;
}

template<typename T>
std::vector< std::pair<point<T,2>, T> >
integrate_interface(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::cell_type& cl,
                    size_t degree)
{
    assert( is_cut(msh, cl) );

    std::vector< std::pair<point<T,2>, T> > ret;

    auto qps = edge_quadrature<T>(degree);  // <-- This has to be changed! Slows down everything!

    for (size_t i = 1; i < cl.user_data.interface.size(); i++)
    {
        auto p0 = cl.user_data.interface.at(i-1);
        auto p1 = cl.user_data.interface.at(i);

        auto scale = p1 - p0;
        auto meas = scale.to_vector().norm();

        for (auto itor = qps.begin(); itor != qps.end(); itor++)
        {
            auto qp = *itor;
            auto p = qp.first.x() * scale + p0;
            auto w = qp.second * meas;

            ret.push_back( std::make_pair(p, w) );
        }
    }

    return ret;
}

template<typename T, typename Function>
std::pair<T, T>
test_integration(const cuthho_mesh<T>& msh, const Function& f)
{
    T int_val = 0.0;
    T circ_int_val = 0.0;

    for (auto& cl : msh.cells)
    {
        bool in_negative_side = (location(msh, cl) == element_location::IN_NEGATIVE_SIDE);
        bool on_interface = (location(msh, cl) == element_location::ON_INTERFACE);
        if ( !in_negative_side && !on_interface )
            continue;

        auto qpts = integrate(msh, cl, 1, element_location::IN_NEGATIVE_SIDE);

        for (auto& qp : qpts)
            int_val += qp.second * f(qp.first);

        if (on_interface)
        {
            auto iqpts = integrate_interface(msh, cl, 1);
            for (auto& qp : iqpts)
                circ_int_val += qp.second * f(qp.first);
        }

    }

    return std::make_pair(int_val, circ_int_val);
}





int main(int argc, char **argv)
{
    using RealType = double;
    size_t degree = 0;

    mesh_init_params<RealType> mip;
    mip.Nx = 5;
    mip.Ny = 5;

    int ch;
    while ( (ch = getopt(argc, argv, "k:M:N:t")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                break;

            case 'M':
                mip.Nx = atoi(optarg);
                break;

            case 'N':
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

    cuthho_mesh<RealType> msh(mip);

    silo_database silo;
    silo.create("cuthho_square.silo");
    silo.add_mesh(msh, "mesh");

    RealType radius = 0.35;

    auto level_set_function = [&](const typename cuthho_mesh<RealType>::point_type pt) -> RealType {
        auto x = pt.x();
        auto y = pt.y();

        auto alpha = 0.5;
        auto beta = 0.5;

        return (x-alpha)*(x-alpha) + (y-beta)*(y-beta) - radius*radius;
        //return std::abs(x+y) + std::abs(x-y) - (1./3.);
    };

    detect_node_position(msh, level_set_function);
    detect_cut_faces(msh, level_set_function);
    detect_cut_cells(msh, level_set_function);
    refine_interface(msh, level_set_function, 3);
    dump_mesh(msh);
    test_triangulation(msh);

    auto intfunc = [](const point<RealType,2>& pt) -> RealType {
        return 1;
    };
    auto ints = test_integration(msh, intfunc);
    auto expval = radius*radius*M_PI;
    std::cout << "Integral relative error: " << 100*std::abs(ints.first-expval)/expval << "%" <<std::endl;
    expval = 2*M_PI*radius;
    std::cout << "Integral relative error: " << 100*std::abs(ints.second-expval)/expval << "%" <<std::endl;


    std::vector<RealType> cut_cell_markers;
    for (auto& cl : msh.cells)
    {
        if ( location(msh, cl) == element_location::IN_POSITIVE_SIDE )
            cut_cell_markers.push_back(1.0);
        else if ( location(msh, cl) == element_location::IN_NEGATIVE_SIDE )
            cut_cell_markers.push_back(-1.0);
        else if ( location(msh, cl) == element_location::ON_INTERFACE )
            cut_cell_markers.push_back(0.0);
        else
            throw std::logic_error("shouldn't have arrived here...");
    }
    silo.add_variable("mesh", "cut_cells", cut_cell_markers.data(), cut_cell_markers.size(), zonal_variable_t);


    std::vector<RealType> level_set_vals;
    for (auto& pt : msh.points)
        level_set_vals.push_back( level_set_function(pt) );
    silo.add_variable("mesh", "level_set", level_set_vals.data(), level_set_vals.size(), nodal_variable_t);

    std::vector<RealType> node_pos;
    for (auto& n : msh.nodes)
        node_pos.push_back( location(msh, n) == element_location::IN_POSITIVE_SIDE ? +1.0 : -1.0 );
    silo.add_variable("mesh", "node_pos", node_pos.data(), node_pos.size(), nodal_variable_t);


    auto rhs_fun = [](const typename cuthho_mesh<RealType>::point_type& pt) -> RealType {
        return 2.0 * M_PI * M_PI * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    auto sol_fun = [](const typename cuthho_mesh<RealType>::point_type& pt) -> RealType {
        return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    for (auto& fc : msh.faces)
    {
        if (fc.is_boundary)
            fc.bndtype = boundary::DIRICHLET;
    }

    auto assembler = make_assembler(msh, degree);
    for (auto& cl : msh.cells)
    {
        auto gr = make_operator(msh, cl, degree);
        Matrix<RealType, Dynamic, Dynamic> stab = make_stabilization(msh, cl, degree);
        Matrix<RealType, Dynamic, Dynamic> lc = gr.second + stab;
        Matrix<RealType, Dynamic, 1> f = make_rhs(msh, cl, degree, rhs_fun);
        assembler.assemble(msh, cl, lc, f, sol_fun);
    }

    assembler.finalize();

    std::cout << "System unknowns: " << assembler.LHS.rows() << std::endl;
    std::cout << "Cells: " << msh.cells.size() << std::endl;
    std::cout << "Faces: " << msh.faces.size() << std::endl;

//#if 0
    SparseLU<SparseMatrix<RealType>>  solver;

    solver.analyzePattern(assembler.LHS);
    solver.factorize(assembler.LHS);
    Matrix<RealType, Dynamic, 1> sol = solver.solve(assembler.RHS);
//#endif
#if 0
    Matrix<RealType, Dynamic, 1> sol = Matrix<RealType, Dynamic, 1>::Zero(assembler.RHS.rows());
    conjugated_gradient(assembler.LHS, assembler.RHS, sol);
#endif

    std::ofstream ofs("solution.dat");
    RealType error_int = 0.0;
    RealType error_mm = 0.0;

    std::vector<RealType>   error_int_percell;
    std::vector<RealType>   error_mm_percell;
    std::vector<RealType>   solution_val;
    std::vector<RealType>   stab_val;

    size_t cell_i = 0;
    for (auto& cl : msh.cells)
    {
        cell_basis<cuthho_mesh<RealType>, RealType> cb(msh, cl, degree);
        auto cbs = cb.size();

        Matrix<RealType, Dynamic, 1> cdofs = sol.block(cell_i*cbs, 0, cbs, 1);

        auto tps = make_test_points(msh, cl);
        for (auto& tp : tps)
        {
            Matrix<RealType, Dynamic, 1> phi = cb.eval_basis(tp);
            auto val = cdofs.dot( phi );
            ofs << tp.x() << " " << tp.y() << " " << val << std::endl;
        }

        RealType e_int_c = 0.0;
        auto qps = integrate(msh, cl, 2*degree);
        for (auto& qp : qps)
        {
            Matrix<RealType, Dynamic, 1> phi = cb.eval_basis(qp.first);
            auto val = cdofs.dot( phi );
            auto real_val = sol_fun( qp.first );
            e_int_c += qp.second*(real_val - val)*(real_val - val);
        }
        error_int_percell.push_back(e_int_c);
        error_int += e_int_c;

        Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, degree);
        Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cl, degree, sol_fun);
        Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
        Matrix<RealType, Dynamic, 1> diff = real_dofs - cdofs;
        RealType e_mm_c = diff.dot(mass*diff);
        error_mm_percell.push_back(e_mm_c);
        error_mm += e_mm_c;

        {
            auto bar = barycenter(msh, cl);
            Matrix<RealType, Dynamic, 1> phi = cb.eval_basis(bar);
            auto val = cdofs.dot( phi );
            solution_val.push_back(val);
        }

        Matrix<RealType, Dynamic, Dynamic> stab = make_stabilization(msh, cl, degree);
        Matrix<RealType, Dynamic, 1> c_alldofs = assembler.take_local_data(msh, cl, sol, sol_fun);
        auto sv = c_alldofs.dot(stab*c_alldofs);
        stab_val.push_back(sv);

        cell_i++;
    }
    ofs.close();

    silo.add_variable("mesh", "err_int_c", error_int_percell.data(), error_int_percell.size(), zonal_variable_t);
    silo.add_variable("mesh", "err_mm_c", error_mm_percell.data(), error_mm_percell.size(), zonal_variable_t);
    silo.add_variable("mesh", "u", solution_val.data(), solution_val.size(), zonal_variable_t);
    silo.add_variable("mesh", "stab_val", stab_val.data(), stab_val.size(), zonal_variable_t);

    std::cout << "L2-error (int): " << sqrt(error_int) << std::endl;
    std::cout << "L2-error (mm) : " << sqrt(error_mm) << std::endl;

    return 0;
}
