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
#include "dataio/silo_io.hpp"

#include "methods/hho"
#include "methods/cuthho"




template<typename T, typename Function>
point<T, 2>
find_zero_crossing(const point<T,2>& p0, const point<T,2>& p1, const Function& level_set_function,
                   const T& threshold)
{
    /* !!! We assume that the level set function *has* a zero crossing
     * between p0 and p1 !!! */
    auto pa = p0;
    auto pb = p1;
    auto pm = (pa+pb)/2.0;
    auto pm_prev = pm;

    T x_diff_sq, y_diff_sq;

    /* A threshold of 1/1000 the diameter of the element is considered
     * acceptable. Since with 20 iterations we reduce the error by 1024
     * and the worst case is that the two points are at the opposite sides
     * of the element, we put this limit. */
    size_t max_iter = 20;

    do {
        auto la = level_set_function(pa);
        auto lb = level_set_function(pb);
        auto lm = level_set_function(pm);

        if ( (lb >= 0 && lm >= 0) || (lb < 0 && lm < 0) )
        {   /* intersection is between pa and pm */
            pm_prev = pm;
            pb = pm;
            pm = (pa+pb)/2.0;
        }
        else
        {   /* intersection is between pm and pb */
            pm_prev = pm;
            pa = pm;
            pm = (pa+pb)/2.0;
        }

        x_diff_sq = (pm_prev.x() - pm.x()) * (pm_prev.x() - pm.x());
        y_diff_sq = (pm_prev.y() - pm.y()) * (pm_prev.y() - pm.y());

    } while ( (sqrt(x_diff_sq + y_diff_sq) > threshold) && max_iter-- );

    return pm;

    /* Affine zero crossing was like that: */
    //auto t = l0/(l0-l1);
    //auto ip = (pts[1] - pts[0]) * t + pts[0];
}

template<typename T, typename Function>
void
detect_cut_faces(cuthho_mesh<T>& msh, const Function& level_set_function)
{
    for (auto& fc : msh.faces)
    {
        auto pts = points(msh, fc);
        auto l0 = level_set_function(pts[0]);
        auto l1 = level_set_function(pts[1]);
        if (l0 >= 0 && l1 >= 0)
        {
            fc.user_data.location = element_location::IN_POSITIVE_SIDE;
            continue;
        }

        if (l0 < 0 && l1 < 0)
        {
            fc.user_data.location = element_location::IN_NEGATIVE_SIDE;
            continue;
        }

        auto threshold = diameter(msh, fc) / 1000.0;
        auto pm = find_zero_crossing(pts[0], pts[1], level_set_function, threshold);

        /* If node 0 is in the negative region, mark it as node inside, otherwise mark node 1 */
        fc.user_data.node_inside = ( l0 < 0 ) ? 0 : 1;
        fc.user_data.location = element_location::ON_INTERFACE;
        fc.user_data.intersection_point = pm;
    }
}

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
bool
is_cut(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::cell_type& cl)
{
    assert(cl.user_data.location != element_location::UNDEF);
    return cl.user_data.location == element_location::ON_INTERFACE;
}

template<typename T>
bool
is_cut(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::face_type& fc)
{
    assert(fc.user_data.location != element_location::UNDEF);
    return fc.user_data.location == element_location::ON_INTERFACE;
}

template<typename T>
element_location
location(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::cell_type& cl)
{
    assert(cl.user_data.location != element_location::UNDEF);
    return cl.user_data.location;
}

template<typename T>
element_location
location(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::face_type& fc)
{
    assert(fc.user_data.location != element_location::UNDEF);
    return fc.user_data.location;
}



template<typename T, typename Function>
void
detect_cut_cells(cuthho_mesh<T>& msh, const Function& level_set_function)
{
    typedef typename cuthho_mesh<T>::face_type  face_type;
    typedef typename cuthho_mesh<T>::point_type point_type;

    size_t cell_i = 0;
    for (auto& cl : msh.cells)
    {
        auto fcs = faces(msh, cl);
        std::array< std::pair<size_t, point_type>, 2 >  cut_faces;

        size_t k = 0;
        for (size_t i = 0; i < fcs.size(); i++)
        {
            if ( is_cut(msh, fcs[i]) )
                cut_faces.at(k++) = std::make_pair(i, fcs[i].user_data.intersection_point);
        }

        /* If a face is cut, the cells that own the face are cut. Is this
         * unconditionally true? It should...fortunately this isn't avionics
         * software */

        if (k == 0)
        {
            auto is_positive = [&](const point_type& pt) -> bool {
                return level_set_function(pt) > 0;
            };

            auto pts = points(msh, cl);
            if ( std::all_of(pts.begin(), pts.end(), is_positive) )
                cl.user_data.location = element_location::IN_POSITIVE_SIDE;
            else
                cl.user_data.location = element_location::IN_NEGATIVE_SIDE;
        }

        if (k == 2)
        {
            cl.user_data.location = element_location::ON_INTERFACE;
            auto p0 = cut_faces[0].second;
            auto p1 = cut_faces[1].second;
            auto pt = p1 - p0;
            auto pn = p0 + point<T,2>(-pt.y(), pt.x());

            if ( level_set_function(pn) >= 0 )
            {
                cl.user_data.p0 = p1;
                cl.user_data.p1 = p0;
            }
            else
            {
                cl.user_data.p0 = p0;
                cl.user_data.p1 = p1;
            }

            cl.user_data.interface.push_back(cl.user_data.p0);
            cl.user_data.interface.push_back(cl.user_data.p1);
        }

        if ( k != 0 && k != 2 )
            throw std::logic_error("invalid number of cuts in cell");

        cell_i++;
    }
}

template<typename T, typename Function>
void
refine_interface(cuthho_mesh<T>& msh, typename cuthho_mesh<T>::cell_type& cl,
                 const Function& level_set_function, size_t min, size_t max)
{
    if ( (max-min) < 2 )
        return;

    typedef typename cuthho_mesh<T>::point_type     point_type;

    size_t mid = (max+min)/2;
    auto p0 = cl.user_data.interface.at(min);
    auto p1 = cl.user_data.interface.at(max);
    auto pm = (p0+p1)/2.0;
    auto pt = p1 - p0;
    auto pn = point_type(-pt.y(), pt.x());
    auto ps1 = pm + pn;
    auto ps2 = pm - pn;

    auto lm = level_set_function(pm);
    auto ls1 = level_set_function(ps1);
    auto ls2 = level_set_function(ps2);

    point_type ip;

    if ( !((lm >= 0 && ls1 >= 0) || (lm < 0 && ls1 < 0)) )
    {
        auto threshold = diameter(msh, cl) / 1000.0;
        ip = find_zero_crossing(pm, ps1, level_set_function, threshold);
    }
    else if ( !((lm >= 0 && ls2 >= 0) || (lm < 0 && ls2 < 0)) )
    {
        auto threshold = diameter(msh, cl) / 1000.0;
        ip = find_zero_crossing(pm, ps2, level_set_function, threshold);
    }
    else
        throw std::logic_error("interface not found in search range");

    cl.user_data.interface.at(mid) = ip;

    refine_interface(msh, cl, level_set_function, min, mid);
    refine_interface(msh, cl, level_set_function, mid, max);
}

template<typename T, typename Function>
void
refine_interface(cuthho_mesh<T>& msh, const Function& level_set_function, size_t levels)
{
    if (levels == 0)
        return;

    size_t interface_points = iexp_pow(2, levels);

    for (auto& cl : msh.cells)
    {
        if ( !is_cut(msh, cl) )
            continue;

        cl.user_data.interface.resize(interface_points+1);
        cl.user_data.interface.at(0)                = cl.user_data.p0;
        cl.user_data.interface.at(interface_points) = cl.user_data.p1;

        refine_interface(msh, cl, level_set_function, 0, interface_points);
    }
}


template<typename T>
void dump_mesh(const cuthho_mesh<T>& msh)
{
    std::ofstream ofs("mesh_dump.m");
    size_t i = 0;
    ofs << "hold on;" << std::endl;
    for (auto& fc : msh.faces)
    {
        auto pts = points(msh, fc);
        if (fc.is_boundary)
            ofs << "line([" << pts[0].x() << ", " << pts[1].x() << "], [" << pts[0].y() << ", " << pts[1].y() << "], 'Color', 'r');" << std::endl;
        else if ( is_cut(msh, fc) )
            ofs << "line([" << pts[0].x() << ", " << pts[1].x() << "], [" << pts[0].y() << ", " << pts[1].y() << "], 'Color', 'g');" << std::endl;
        else
            ofs << "line([" << pts[0].x() << ", " << pts[1].x() << "], [" << pts[0].y() << ", " << pts[1].y() << "], 'Color', 'k');" << std::endl;

        auto bar = barycenter(msh, fc);
        ofs << "text(" << bar.x() << ", " << bar.y() << ", '" << i << "');" << std::endl;

        i++;
    }

    for (auto& cl : msh.cells)
    {
        if ( is_cut(msh, cl) )
        {
            auto p0 = cl.user_data.p0;
            auto p1 = cl.user_data.p1;
            auto q = p1 - p0;
            ofs << "quiver(" << p0.x() << ", " << p0.y() << ", " << q.x() << ", " << q.y() << ", 0)" << std::endl;

            for (auto& ip : cl.user_data.interface)
                ofs << "plot(" << ip.x() << ", " << ip.y() << ", '*k');" << std::endl;
        }
    }

    ofs << "t = linspace(0,2*pi,1000);plot(sqrt(0.15)*cos(t)+0.5, sqrt(0.15)*sin(t)+0.5)" << std::endl;

    ofs.close();
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

    auto level_set_function = [](const typename cuthho_mesh<RealType>::point_type pt) -> RealType {
        auto x = pt.x();
        auto y = pt.y();
        auto alpha = 0.5;
        auto beta = 0.5;

        return (x-alpha)*(x-alpha) + (y-beta)*(y-beta) - 0.15;
    };

    detect_cut_faces(msh, level_set_function);
    detect_cut_cells(msh, level_set_function);
    refine_interface(msh, level_set_function, 3);
    dump_mesh(msh);


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
