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





template<typename Mesh, typename Function>
std::pair<bool, typename Mesh::point_type>
find_face_interface_intersection(const Mesh& msh, const typename Mesh::face_type& fc,
                                 const Function& level_set_function)
{
    auto pts = points(msh, fc);
    auto l0 = level_set_function(pts[0]);
    auto l1 = level_set_function(pts[1]);
    if ( (l0 >= 0 && l1 >= 0) || (l0 < 0 && l1 < 0) )
        return std::make_pair( false, typename Mesh::point_type() );

    auto t = l0/(l0-l1);
    auto ip = (pts[1] - pts[0]) * t + pts[0];
    return std::make_pair(true, ip);
}

#if 0
template<typename T, typename Function>
std::shared_ptr<cell_unfitted_info<T>>
make_cell_unfitted_info(mesh<T>& msh, const typename mesh<T>::cell_type& cl,
                        const Function& level_set_function)
{
    auto pts = points(msh, cl);
    auto fcs = faces(msh, cl);

    std::array< std::pair<bool, point<T,2>>, 4> intersections;

    auto tf = [&](const typename mesh<T>::face_type& fc) -> auto {
        return find_face_interface_intersection(msh, fc, level_set_function);
    };
    std::transform(fcs.begin(), fcs.end(), intersections.begin(), tf);

    auto is_intersected = [](const std::pair<bool, typename mesh<T>::point_type>& int_pt) -> auto {
        return int_pt.first;
    };
    auto num_cuts = std::count_if(intersections.begin(), intersections.end(), is_intersected);

    if (num_cuts == 0)  /* The cell is not cut, don't return anything */
        return nullptr;

    if (num_cuts != 2)
        throw std::logic_error("invalid number of cuts in cell");

    auto ui = std::make_shared<cell_unfitted_info<T>>();
    ui->intersections = intersections;

    return ui;
}


template<typename T, typename Function>
void
fill_unfitted_info(mesh<T>& msh, const Function& level_set_function)
{
    /*
    for (auto& cl : msh.cells)
        cl.unfitted_info = make_cell_unfitted_info(msh, cl, level_set_function);
    */
}
#endif

template<typename T>
struct cell_info
{
    //bool        cut_by_interface;

    std::array< std::pair<bool, point<T,2>>, 4 >        cuts;
    /*std::vector<*/ std::array<point<T,2>, 2> /*>*/            additional_faces;

    cell_info() //:
        //cut_by_interface(false)
    {
        for (auto& c : cuts)
            c.first = false;
    }

    bool is_cut_by_interface(void) const
    {
        auto is_cut = [&](const std::pair<bool, point<T,2>>& c) -> bool {
            return c.first;
        };
        return std::any_of(cuts.begin(), cuts.end(), is_cut);
    }
};

template<typename T>
struct face_info
{
    bool                    intersected_by_interface;
    point<T, 2>             intersection_point;

    face_info() :
        intersected_by_interface(false)
    {}
};

template<typename T>
class element_info
{
    typedef typename mesh<T>::cell_type     cell_type;
    typedef typename mesh<T>::face_type     face_type;
    typedef typename mesh<T>::point_type    point_type;

    std::vector<cell_info<T>>               cell_information;
    std::vector<face_info<T>>               face_information;

    template<typename Function>
    std::pair<bool, point_type>
    find_face_interface_intersection(const mesh<T>& msh, const face_type& fc, const Function& level_set_function)
    {
        auto pts = points(msh, fc);
        auto l0 = level_set_function(pts[0]);
        auto l1 = level_set_function(pts[1]);
        if ( (l0 >= 0 && l1 >= 0) || (l0 < 0 && l1 < 0) )
            return std::make_pair( false, point_type() );

        auto t = l0/(l0-l1);
        auto ip = (pts[1] - pts[0]) * t + pts[0];
        return std::make_pair(true, ip);
    }

    template<typename Function>
    void
    detect_cut_cells(const mesh<T>& msh, const Function& level_set_function)
    {
        size_t cell_i = 0;
        for (auto& cl : msh.cells)
        {
            auto pts = points(msh, cl);
            std::vector<int> signs;
            signs.resize(pts.size());

            auto fcs = faces(msh, cl);
            auto tf = [&](const face_type& fc) -> auto {
                return find_face_interface_intersection(msh, fc, level_set_function);
            };
            std::transform(fcs.begin(), fcs.end(), cell_information.at(cell_i).cuts.begin(), tf);

            auto count_cut = [](const std::pair<bool, point_type>& c) -> auto {
                return c.first;
            };
            auto num_cuts = std::count_if(cell_information.at(cell_i).cuts.begin(),
                                          cell_information.at(cell_i).cuts.end(), count_cut);

            if (num_cuts != 0 && num_cuts != 2)
            {
                std::cout << num_cuts << std::endl;
                throw std::logic_error("invalid number of cuts in cell");
            }

            if (num_cuts == 2)
            {
                for (size_t i = 0, k = 0; i < 4; i++)
                {
                    if ( cell_information.at(cell_i).cuts[i].first )
                        cell_information.at(cell_i).additional_faces[k++] = cell_information.at(cell_i).cuts[i].second;
                }

                auto p0 = cell_information.at(cell_i).additional_faces[0];
                auto p1 = cell_information.at(cell_i).additional_faces[1];

                auto pt = p1 - p0;
                auto pn = p0 + point<T,2>(-pt.y(), pt.x());

                if ( level_set_function(pn) >= 0 )
                    std::swap(cell_information.at(cell_i).additional_faces[0], cell_information.at(cell_i).additional_faces[1]);
            }

            cell_i++;
        }
    }

    template<typename Function>
    void
    detect_face_intersection_points(const mesh<T>& msh, const Function& level_set_function)
    {
        size_t face_i = 0;
        for (auto& fc : msh.faces)
        {
            auto fi = find_face_interface_intersection(msh, fc, level_set_function);

            face_information.at(face_i).intersected_by_interface = fi.first;
            face_information.at(face_i).intersection_point = fi.second;

            face_i++;
        }
    }

public:
    template<typename Function>
    element_info(const mesh<T>& msh, const Function& level_set_function)
    {
        cell_information.resize(msh.cells.size());
        face_information.resize(msh.faces.size());
        detect_cut_cells(msh, level_set_function);
        detect_face_intersection_points(msh, level_set_function);
    }

    cell_info<T>
    info(const mesh<T>& msh, const cell_type& cl) const
    {
        return cell_information.at( offset(msh, cl) );
    }

    face_info<T>
    info(const mesh<T>& msh, const face_type& fc) const
    {
        return face_information.at( offset(msh, fc) );
    }
};

template<typename T, typename Function>
auto make_element_info(const mesh<T>& msh, const Function& level_set_function)
{
    return element_info<T>(msh, level_set_function);
}

/*****************************************************************************
 *   SILO stuff
 *****************************************************************************/



/*****************************************************************************
 *   Geometry
 *****************************************************************************/



/*****************************************************************************
 *   Bases
 *****************************************************************************/


/*****************************************************************************
 *   Quadrature rules
 *****************************************************************************/




template<typename T>
std::vector<std::pair<point<T,2>, T>>
integrate_cut(const mesh<T>& msh, const typename mesh<T>::cell_type& cl, const cell_info<T>& ci, size_t degree)
{
    typedef typename mesh<T>::point_type    point_type;

    auto qps = edge_quadrature<T>(degree);
    auto pts = points(msh, cl);

    auto scale_x = (pts[1] - pts[0]).x();
    auto scale_y = (pts[3] - pts[0]).y();
    auto meas = scale_x * scale_y;

    std::vector<std::pair<point<T,2>, T>> ret;

    for (auto jtor = qps.begin(); jtor != qps.end(); jtor++)
    {
        for (auto itor = qps.begin(); itor != qps.end(); itor++)
        {
            auto qp_x = *itor;
            auto qp_y = *jtor;

            auto px = qp_x.first.x() * scale_x + pts[0].x();
            auto py = qp_y.first.x() * scale_y + pts[0].y();

            auto w = qp_x.second * qp_y.second * meas;

            ret.push_back( std::make_pair( point_type(px, py), w ) );
        }
    }

    return ret;
}

template<typename T>
std::vector<std::pair<point<T,2>, T>>
integrate(const mesh<T>& msh, const typename mesh<T>::cell_type& cl, const cell_info<T>& ci, size_t degree)
{
    if ( ci.is_cut_by_interface() )
        return integrate_cut(msh, cl, ci, degree);

    return integrate(msh, cl, degree);
}



/*****************************************************************************
 *   Utilities
 *****************************************************************************/


/*****************************************************************************
 *   HHO stuff
 *****************************************************************************/


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
                auto gr = make_operator(msh, cl, k);
                Matrix<RealType, Dynamic, Dynamic> stab = make_stabilization(msh, cl, k);
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






template<typename T>
bool
is_cut(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::cell_type& cl)
{
    return cl.user_data.num_cuts != 0;
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
        auto pts = points(msh, cl);
        std::vector<int> signs;
        signs.resize(pts.size());

        auto fcs = faces(msh, cl);
        auto tf = [&](const face_type& fc) -> auto {
            return find_face_interface_intersection(msh, fc, level_set_function);
        };

        std::transform(fcs.begin(), fcs.end(), cl.user_data.cuts.begin(), tf);


        auto count_cut = [](const std::pair<bool, point_type>& c) -> auto {
            return c.first;
        };
        cl.user_data.num_cuts = std::count_if(cl.user_data.cuts.begin(),
                                              cl.user_data.cuts.end(),
                                              count_cut);


        if (cl.user_data.num_cuts != 0 && cl.user_data.num_cuts != 2)
        {
            std::cout << cl.user_data.num_cuts << std::endl;
            throw std::logic_error("invalid number of cuts in cell");
        }

        if (cl.user_data.num_cuts == 2)
        {
            for (size_t i = 0, k = 0; i < 4; i++)
            {
                if ( cl.user_data.cuts[i].first )
                    cl.user_data.additional_faces[k++] = cl.user_data.cuts[i].second;
            }

            auto p0 = cl.user_data.additional_faces[0];
            auto p1 = cl.user_data.additional_faces[1];

            auto pt = p1 - p0;
            auto pn = p0 + point<T,2>(-pt.y(), pt.x());

            if ( level_set_function(pn) >= 0 )
                std::swap(cl.user_data.additional_faces[0], cl.user_data.additional_faces[1]);
        }

        cell_i++;
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
        else
            ofs << "line([" << pts[0].x() << ", " << pts[1].x() << "], [" << pts[0].y() << ", " << pts[1].y() << "]);" << std::endl;

        auto bar = barycenter(msh, fc);
        ofs << "text(" << bar.x() << ", " << bar.y() << ", '" << i << "');" << std::endl;

        i++;
    }

    for (auto& cl : msh.cells)
    {
        if ( is_cut(msh, cl) )
        {
            auto p0 = cl.user_data.additional_faces[0];
            auto p1 = cl.user_data.additional_faces[1];
            auto q = p1 - p0;
            ofs << "quiver(" << p0.x() << ", " << p0.y() << ", " << q.x() << ", " << q.y() << ", 0)" << std::endl;;
        }
    }

    ofs.close();
}



int main(int argc, char **argv)
{
    using RealType = double;
    size_t degree = 0;

    mesh_init_params<RealType> mip;
    mip.Nx = 5;
    mip.Ny = 5;

    bool testing_only = false;

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

            case 't':
                testing_only = true;
                break;

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    if (testing_only)
    {
        convergence_test_params ctp;
        ctp.deg_min = 5;
        test_method_convergence(ctp);
        return 0;
    }

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

    detect_cut_cells(msh, level_set_function);
    dump_mesh(msh);


    std::vector<RealType> cut_cell_markers;
    for (auto& cl : msh.cells)
    {
        if ( is_cut(msh, cl) )
            cut_cell_markers.push_back(1.0);
        else
            cut_cell_markers.push_back(0.0);
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
