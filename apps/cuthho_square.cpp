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




template<typename T>
struct circle_level_set
{
    T radius, alpha, beta;

    circle_level_set(T r, T a, T b)
        : radius(r), alpha(a), beta(b)
    {}

    T operator()(const point<T,2>& pt) const
    {
        auto x = pt.x();
        auto y = pt.y();

        return (x-alpha)*(x-alpha) + (y-beta)*(y-beta) - radius*radius;
    }

    Eigen::Matrix<T,2,1> gradient(const point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;
        ret(0) = 2*pt.x() - 2*alpha;
        ret(1) = 2*pt.y() - 2*beta;
        return ret;
    }

    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;

        ret = gradient(pt);
        return ret/ret.norm();
    }

};


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




















template<typename T, typename ET>
void test_triangulation(const cuthho_mesh<T, ET>& msh)
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
struct material_parameters
{
    T kappa_1, kappa_2;

    material_parameters() : kappa_1(1.0), kappa_2(1.0) {}
};

template<typename T, typename ET>
T
cell_eta(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl)
{
    return 40*measure(msh, cl);
}


template<typename T, typename ET, typename Function>
std::pair<   Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>  >
make_hho_cut_laplacian(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
                       const Function& level_set_function, hho_degree_info di)
{
    assert( is_cut(msh, cl) );

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    cell_basis<cuthho_mesh<T, ET>,T>     cb(msh, cl, recdeg);

    auto rbs = cell_basis<cuthho_mesh<T, ET>,T>::size(recdeg);
    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);

    Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    Matrix<T, Dynamic, Dynamic> gr_lhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
    Matrix<T, Dynamic, Dynamic> gr_rhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, cbs + 4*fbs);

    /* Cell term (cut) */
    auto qps = integrate(msh, cl, 2*recdeg, element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : qps)
    {
        auto dphi = cb.eval_gradients(qp.first);
        stiff += qp.second * dphi * dphi.transpose();
    }

    auto hT = measure(msh, cl);

    /* Interface term */
    auto iqps = integrate_interface(msh, cl, 2*recdeg);
    for (auto& qp : iqps)
    {
        auto phi    = cb.eval_basis(qp.first);
        auto dphi   = cb.eval_gradients(qp.first);
        auto n      = level_set_function.normal(qp.first);

        stiff -= qp.second * phi * (dphi * n).transpose();
        stiff -= qp.second * (dphi * n) * phi.transpose();
        stiff += qp.second * phi * phi.transpose();// * cell_eta(msh, cl) / hT;
    }

    gr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
    gr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

    auto ns = normals(msh, cl);
    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];
        auto n = ns[i];

        face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg);
        /* Terms on faces */
        auto qps = integrate(msh, fc, 2*recdeg, element_location::IN_NEGATIVE_SIDE);
        for (auto& qp : qps)
        {
            auto c_phi = cb.eval_basis(qp.first);
            auto f_phi = fb.eval_basis(qp.first);
            auto r_dphi_tmp = cb.eval_gradients(qp.first);
            auto r_dphi = r_dphi_tmp.block(1, 0, rbs-1, 2);
            gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.second * (r_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs-1, cbs) -= qp.second * (r_dphi * n) * c_phi.transpose();
        }
    }

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}


template<typename T, typename ET>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>
make_hho_cut_stabilization(const cuthho_mesh<T, ET>& msh,
                           const typename cuthho_mesh<T, ET>::cell_type& cl,
                           const hho_degree_info& di)
{
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);

    auto fcs = faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+4*fbs, cbs+4*fbs);
    Matrix<T, Dynamic, Dynamic> If = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);

    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, celdeg);

    auto hT = measure(msh, cl);

    for (size_t i = 0; i < 4; i++)
    {
        auto fc = fcs[i];
        face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg);

        Matrix<T, Dynamic, Dynamic> oper = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs+4*fbs);
        Matrix<T, Dynamic, Dynamic> mass = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> trace = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs);

        oper.block(0, cbs+i*fbs, fbs, fbs) = -If;

        auto qps = integrate(msh, fc, 2*facdeg, element_location::IN_NEGATIVE_SIDE);
        for (auto& qp : qps)
        {
            auto c_phi = cb.eval_basis(qp.first);
            auto f_phi = fb.eval_basis(qp.first);

            mass += qp.second * f_phi * f_phi.transpose();
            trace += qp.second * f_phi * c_phi.transpose();
        }

        if (qps.size() == 0) /* Avoid to invert a zero matrix */
            continue;

        oper.block(0, 0, fbs, cbs) = mass.llt().solve(trace);

        data += oper.transpose() * mass * oper * (1./hT);
    }

    return data;
}


template<typename T, typename ET, typename Function>
std::pair<   Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>  >
make_operator(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
              const Function& level_set_function, hho_degree_info hdi)
{
    if ( is_cut(msh, cl) )
        return make_hho_cut_laplacian(msh, cl, level_set_function, hdi);
    else
        return make_hho_laplacian(msh, cl, hdi);
}



template<typename T, typename ET>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>
make_stabilization(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl, hho_degree_info hdi)
{
    if ( is_cut(msh, cl) )
        return make_hho_cut_stabilization(msh, cl, hdi);
    else
        return make_hho_naive_stabilization(msh, cl, hdi);
}


template<typename T, typename ET, typename F1, typename F2, typename F3>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, 1>
make_rhs(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
         size_t degree, const F1& f, const element_location where, const F2& level_set_function, const F3& bcs)
{
    if ( location(msh, cl) == where )
        return make_rhs(msh, cl, degree, f);
    else if ( location(msh, cl) == element_location::ON_INTERFACE )
    {
        cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, degree);
        auto cbs = cb.size();

        auto hT = measure(msh, cl);

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

        auto qps = integrate(msh, cl, 2*degree, where);
        for (auto& qp : qps)
        {
            auto phi = cb.eval_basis(qp.first);
            ret += qp.second * phi * f(qp.first);
        }

        /*
        auto qpsi = integrate_interface(msh, cl, degree);
        for (auto& qp : qpsi)
        {
            auto phi = cb.eval_basis(qp.first);
            auto dphi = cb.eval_gradients(qp.first);
            auto n = level_set_function.normal(qp.first);

            ret += qp.second * bcs(qp.first) * ( phi * cell_eta(msh, cl)/hT - dphi*n);
        }
        */

        return ret;  
    }
    else
    {
        auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(degree);
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);
        return ret;
    }
}



template<typename T>
std::string quiver(const point<T,2>& p, const Eigen::Matrix<T,2,1>& v)
{
    std::stringstream ss;

    ss << "quiver(" << p.x() << ", " << p.y() << ", ";
    ss << v(0) << ", " << v(1) << ", 0);";

    return ss.str();
}

template<typename T, typename ET, typename Function1, typename Function2>
std::pair<T, T>
test_integration(const cuthho_mesh<T, ET>& msh, const Function1& f, const Function2& level_set_function)
{
    T surf_int_val = 0.0;
    T line_int_val = 0.0;

    std::ofstream ofs("normals.m");

    for (auto& cl : msh.cells)
    {
        bool in_negative_side = (location(msh, cl) == element_location::IN_NEGATIVE_SIDE);
        bool on_interface = (location(msh, cl) == element_location::ON_INTERFACE);
        if ( !in_negative_side && !on_interface )
            continue;

        auto qpts = integrate(msh, cl, 1, element_location::IN_NEGATIVE_SIDE);

        for (auto& qp : qpts)
            surf_int_val += qp.second * f(qp.first);

        if (on_interface)
        {
            auto iqpts = integrate_interface(msh, cl, 1);
            for (auto& qp : iqpts)
            {
                line_int_val += qp.second * f(qp.first);
                auto n = level_set_function.normal(qp.first);
                ofs << quiver(qp.first, n) << std::endl;
            }
        }

    }

    ofs.close();

    std::ofstream ofs_int("face_ints.m");

    for (auto& fc : msh.faces)
    {
        auto qpts = integrate(msh, fc, 2, element_location::IN_NEGATIVE_SIDE);
        for (auto& qp : qpts)
        {
            ofs_int << "hold on;" << std::endl;
            ofs_int << "plot(" << qp.first.x() << ", " << qp.first.y() << ", 'ko');" << std::endl;
        }
    }

    ofs_int.close();

    return std::make_pair(surf_int_val, line_int_val);
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

    cuthho_quad_mesh<RealType> msh(mip);

    RealType radius = 0.35;

    auto level_set_function = circle_level_set<RealType>(radius, 0.5, 0.5);


    detect_node_position(msh, level_set_function);
    detect_cut_faces(msh, level_set_function);
    
    move_nodes(msh, level_set_function);
    detect_cut_faces(msh, level_set_function);
    move_nodes(msh, level_set_function);
    detect_cut_faces(msh, level_set_function);
    
    detect_cut_cells(msh, level_set_function);
    refine_interface(msh, level_set_function, 2);
    dump_mesh(msh);
    test_triangulation(msh);

    silo_database silo;
    silo.create("cuthho_square.silo");
    silo.add_mesh(msh, "mesh");

    auto intfunc = [](const point<RealType,2>& pt) -> RealType {
        return 1;
    };
    auto ints = test_integration(msh, intfunc, level_set_function);
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


    auto rhs_fun = [](const typename cuthho_quad_mesh<RealType>::point_type& pt) -> RealType {
        return 2.0 * M_PI * M_PI * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
        //return 1;
    };

    auto sol_fun = [](const typename cuthho_quad_mesh<RealType>::point_type& pt) -> RealType {
        return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
        //return 0;
    };

    for (auto& fc : msh.faces)
    {
        if (fc.is_boundary)
            fc.bndtype = boundary::DIRICHLET;
    }

    hho_degree_info hdi(degree+1, degree);

    auto assembler = make_assembler(msh, hdi);
    for (auto& cl : msh.cells)
    {
        auto gr = make_operator(msh, cl, level_set_function, hdi);

        Matrix<RealType, Dynamic, Dynamic> stab = make_stabilization(msh, cl, hdi);
        Matrix<RealType, Dynamic, Dynamic> lc = gr.second + stab;

        Matrix<RealType, Dynamic, 1> f = Matrix<RealType, Dynamic, 1>::Zero(lc.rows());

        f = make_rhs(msh, cl, hdi.cell_degree(), rhs_fun, element_location::IN_NEGATIVE_SIDE, level_set_function, sol_fun);

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

    auto optest_fun = [](const typename cuthho_quad_mesh<RealType>::point_type& pt) -> RealType {
        return pt.x() * pt.y();
        //return 0;
    };

    std::vector<RealType>   solution_uT, solution_Ru, optest, optest_r;

    std::ofstream ofs("basis_check.dat");

    size_t cell_i = 0;
    for (auto& cl : msh.cells)
    {
        cell_basis<cuthho_quad_mesh<RealType>, RealType> cb(msh, cl, hdi.cell_degree());
        auto cbs = cb.size();

        cell_basis<cuthho_quad_mesh<RealType>, RealType> rb(msh, cl, hdi.reconstruction_degree());
        auto rbs = rb.size();

        auto gr = make_operator(msh, cl, level_set_function, hdi);
        Matrix<RealType, Dynamic, 1> locdata = assembler.take_local_data(msh, cl, sol, sol_fun);
        Matrix<RealType, Dynamic, 1> cell_dofs = locdata.head(cbs);
        Matrix<RealType, Dynamic, 1> rec_dofs = gr.first * locdata;

        auto bar = barycenter(msh, cl/*, element_location::IN_NEGATIVE_SIDE*/);
        
        Matrix<RealType, Dynamic, 1> c_phi = cb.eval_basis(bar);
        auto c_val = cell_dofs.dot( c_phi );
        solution_uT.push_back(c_val);

        Matrix<RealType, Dynamic, 1> r_phi = rb.eval_basis(bar);
        auto r_val = rec_dofs.dot( r_phi.tail(rbs-1) ) + locdata(0);
        solution_Ru.push_back(r_val);

        Matrix<RealType, Dynamic, 1> proj = project_function(msh, cl, hdi, element_location::IN_NEGATIVE_SIDE, optest_fun);
        Matrix<RealType, Dynamic, 1> prec = gr.first * proj;
        
        auto p_val = proj.head(cbs).dot( c_phi );
        optest.push_back(p_val);

        auto Rp_val = prec.dot( r_phi.tail(rbs-1) ) + proj(0);
        optest_r.push_back(Rp_val);

        std::cout << ( is_cut(msh, cl) ? "\x1b[31mCUT:\x1b[0m" : "\x1b[32mNOT CUT:\x1b[0m" );
        std::cout << " Cell " << cell_i << std::endl;
        std::cout << gr.first << std::endl << std::endl;
        std::cout << proj.transpose() << std::endl << std::endl;
        std::cout << prec.transpose() << std::endl;

        auto tps = make_test_points(msh, cl);
        for (auto& tp : tps)
        {
            auto t_phi = rb.eval_basis(tp);

            ofs << tp.x() << " " << tp.y() << " ";
            for (size_t i = 0; i < rbs; i++)
                ofs << t_phi(i) << " ";

            ofs << prec.dot( t_phi.tail(rbs-1) ) + proj(0) << std::endl;
        }

        cell_i++;
    }

    ofs.close();

    silo.add_variable("mesh", "uT", solution_uT.data(), solution_uT.size(), zonal_variable_t);
    silo.add_variable("mesh", "Ru", solution_Ru.data(), solution_Ru.size(), zonal_variable_t);
    silo.add_variable("mesh", "optest", optest.data(), optest.size(), zonal_variable_t);
    silo.add_variable("mesh", "optest_r", optest_r.data(), optest_r.size(), zonal_variable_t);


    return 0;
}
