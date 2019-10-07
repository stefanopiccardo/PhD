/*
 *       /\        Guillaume Delay 2018,2019
 *      /__\       guillaume.delay@enpc.fr
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
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <memory>
#include <sstream>
#include <list>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/SparseExtra>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

using namespace Eigen;

#include "core/core"
#include "core/solvers"
#include "dataio/silo_io.hpp"

#include "methods/hho"
#include "methods/cuthho"



///////////////////////   DIVERGENCE RECONSTRUCTION  ///////////////////////////
template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_divergence_reconstruction(const Mesh& msh, const typename Mesh::cell_type& cl,
                                   const hho_degree_info& di)
{
    using T = typename Mesh::coordinate_type;
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();
    const auto pdeg = di.face_degree();

    cell_basis<Mesh,T>                   pb(msh, cl, pdeg);
    vector_cell_basis<Mesh,T>            cb(msh, cl, celdeg);

    auto pbs = cell_basis<Mesh,T>::size(pdeg);
    auto fbs = vector_face_basis<Mesh,T>::size(facdeg);
    auto cbs = vector_cell_basis<Mesh,T>::size(celdeg);

    const auto fcs = faces(msh, cl);
    const auto num_faces = fcs.size();
    const auto ns = normals(msh, cl);

    matrix_type dr_lhs = matrix_type::Zero(pbs, pbs);
    matrix_type dr_rhs = matrix_type::Zero(pbs, cbs + num_faces*fbs);


    const auto qps = integrate(msh, cl, celdeg + pdeg - 1);
    for (auto& qp : qps)
    {
        const auto s_phi  = pb.eval_basis(qp.first);
        const auto s_dphi = pb.eval_gradients(qp.first);
        const auto v_phi  = cb.eval_basis(qp.first);

        dr_lhs += qp.second * s_phi * s_phi.transpose();
        dr_rhs.block(0, 0, pbs, cbs) -= qp.second * s_dphi * v_phi.transpose();
    }


    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc     = fcs[i];
        const auto n      = ns[i];
        vector_face_basis<Mesh,T>            fb(msh, fc, facdeg);

        const auto qps_f = integrate(msh, fc, facdeg + pdeg);
        for (auto& qp : qps_f)
        {
            const auto p_phi = pb.eval_basis(qp.first);
            const auto f_phi = fb.eval_basis(qp.first);

            const Matrix<T, Dynamic, 2> p_phi_n = (p_phi * n.transpose());
            dr_rhs.block(0, cbs + i * fbs, pbs, fbs) += qp.second * p_phi_n * f_phi.transpose();
        }
    }


    assert(dr_lhs.rows() == pbs && dr_lhs.cols() == pbs);
    assert(dr_rhs.rows() == pbs && dr_rhs.cols() == cbs + num_faces * fbs);

    matrix_type oper = dr_lhs.ldlt().solve(dr_rhs);
    matrix_type data = dr_rhs;

    return std::make_pair(oper, data);
}



template<typename T, size_t ET, typename Function>
std::pair<   Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>  >
make_hho_divergence_reconstruction(const cuthho_mesh<T, ET>& msh,
                                   const typename cuthho_mesh<T, ET>::cell_type& cl,
                                   const Function& level_set_function, const hho_degree_info& di,
                                   element_location where, const T coeff)
{
    typedef Matrix<T, Dynamic, Dynamic> matrix_type;

    if ( !is_cut(msh, cl) )
        return make_hho_divergence_reconstruction(msh, cl, di);

    const auto celdeg = di.cell_degree();
    const auto facdeg = di.face_degree();
    const auto pdeg = di.face_degree();

    cell_basis<cuthho_mesh<T, ET>,T>                   pb(msh, cl, pdeg);
    vector_cell_basis<cuthho_mesh<T, ET>,T>            cb(msh, cl, celdeg);

    auto pbs = cell_basis<cuthho_mesh<T, ET>,T>::size(pdeg);
    auto fbs = vector_face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);
    auto cbs = vector_cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);

    const auto fcs = faces(msh, cl);
    const auto num_faces = fcs.size();
    const auto ns = normals(msh, cl);

    matrix_type dr_lhs = matrix_type::Zero(pbs, pbs);
    matrix_type dr_rhs = matrix_type::Zero(pbs, cbs + num_faces*fbs);


    const auto qps = integrate(msh, cl, celdeg + pdeg - 1, where);
    for (auto& qp : qps)
    {
        const auto s_phi  = pb.eval_basis(qp.first);
        const auto s_dphi = pb.eval_gradients(qp.first);
        const auto v_phi  = cb.eval_basis(qp.first);

        dr_lhs += qp.second * s_phi * s_phi.transpose();
        dr_rhs.block(0, 0, pbs, cbs) -= qp.second * s_dphi * v_phi.transpose();
    }


    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc     = fcs[i];
        const auto n      = ns[i];
        cut_vector_face_basis<cuthho_mesh<T, ET>,T>            fb(msh, fc, facdeg, where);

        const auto qps_f = integrate(msh, fc, facdeg + pdeg, where);
        for (auto& qp : qps_f)
        {
            const auto p_phi = pb.eval_basis(qp.first);
            const auto f_phi = fb.eval_basis(qp.first);

            const Matrix<T, Dynamic, 2> p_phi_n = (p_phi * n.transpose());
            dr_rhs.block(0, cbs + i * fbs, pbs, fbs) += qp.second * p_phi_n * f_phi.transpose();
        }
    }


    // interface term
    auto iqp = integrate_interface(msh, cl, celdeg + pdeg, element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : iqp)
    {
        const auto v_phi = cb.eval_basis(qp.first);
        const auto s_phi = pb.eval_basis(qp.first);
        const auto n = level_set_function.normal(qp.first);

        const Matrix<T, Dynamic, 2> p_phi_n = (s_phi * n.transpose());
        dr_rhs.block(0, 0, pbs, cbs) += (1.0 - coeff) * qp.second * p_phi_n * v_phi.transpose();
    }


    assert(dr_lhs.rows() == pbs && dr_lhs.cols() == pbs);
    assert(dr_rhs.rows() == pbs && dr_rhs.cols() == cbs + num_faces * fbs);

    matrix_type oper = dr_lhs.ldlt().solve(dr_rhs);
    matrix_type data = dr_rhs;

    return std::make_pair(oper, data);
}



//// make_hho_divergence_reconstruction_interface
// return the divergence reconstruction for the interface pb
// coeff -> scales the interface term
template<typename T, size_t ET, typename Function>
std::pair<   Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>  >
make_hho_divergence_reconstruction_interface
(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
 const Function& level_set_function, const hho_degree_info& di, element_location where, T coeff)
{

    if ( !is_cut(msh, cl) )
        throw std::invalid_argument("The cell is not cut");

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto celdeg  = di.cell_degree();
    const auto facdeg  = di.face_degree();
    const auto pdeg    = facdeg;

    vector_cell_basis<cuthho_mesh<T, ET>,T>   cb(msh, cl, celdeg);
    cell_basis<cuthho_mesh<T, ET>,T>          pb(msh, cl, pdeg);


    auto cbs = vector_cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = vector_face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);
    auto pbs = cell_basis<cuthho_mesh<T, ET>,T>::size(pdeg);

    const auto num_faces = faces(msh, cl).size();

    matrix_type       rhs_tmp = matrix_type::Zero(pbs, cbs + num_faces * fbs);
    matrix_type        dr_lhs = matrix_type::Zero(pbs, pbs);
    matrix_type        dr_rhs = matrix_type::Zero(pbs, 2*cbs + 2*num_faces * fbs);

    const auto qps = integrate(msh, cl, celdeg - 1 + pdeg, where);
    for (auto& qp : qps)
    {
        const auto v_phi     = cb.eval_basis(qp.first);
        const auto s_phi   = pb.eval_basis(qp.first);
        const auto s_dphi  = pb.eval_gradients(qp.first);

        dr_lhs += qp.second * s_phi * s_phi.transpose();
        rhs_tmp.block(0, 0, pbs, cbs) -= qp.second * s_dphi * v_phi.transpose();
    }

    const auto fcs = faces(msh, cl);
    const auto ns = normals(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = ns[i];
        cut_vector_face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg, where);

        const auto qps_f = integrate(msh, fc, pdeg + std::max(facdeg, celdeg), where);
        for (auto& qp : qps_f)
        {
            const matrix_type f_phi      = fb.eval_basis(qp.first);
            const auto        s_phi      = pb.eval_basis(qp.first);
            const matrix_type qp_s_phi_n = qp.second * s_phi * n.transpose();

            rhs_tmp.block(0, cbs + i * fbs, pbs, fbs) += qp_s_phi_n * f_phi.transpose();
        }
    }

    // term on the interface
    matrix_type        interface_term = matrix_type::Zero(pbs, cbs);
    const auto iqps = integrate_interface(msh, cl, celdeg + pdeg, element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : iqps)
    {
        const auto v_phi        = cb.eval_basis(qp.first);
        const auto s_phi        = pb.eval_basis(qp.first);

        Matrix<T,2,1> n = level_set_function.normal(qp.first);
        const matrix_type qp_s_phi_n = qp.second * s_phi * n.transpose();

        interface_term += qp_s_phi_n * v_phi.transpose();
    }

    // finalize rhs
    if(where == element_location::IN_NEGATIVE_SIDE)
    {
        dr_rhs.block(0, 0, pbs, cbs) += rhs_tmp.block(0, 0, pbs, cbs);
        dr_rhs.block(0, 2*cbs, pbs, num_faces*fbs)
            += rhs_tmp.block(0, cbs, pbs, num_faces*fbs);
        dr_rhs.block(0, 0, pbs, cbs) += (1.0 - coeff) * interface_term;
        dr_rhs.block(0, cbs, pbs, cbs) += coeff * interface_term;
    }
    else if( where == element_location::IN_POSITIVE_SIDE)
    {
        dr_rhs.block(0, cbs, pbs, cbs) += rhs_tmp.block(0, 0, pbs, cbs);
        dr_rhs.block(0, 2*cbs + num_faces*fbs, pbs, num_faces*fbs)
                     += rhs_tmp.block(0, cbs, pbs, num_faces*fbs);
        dr_rhs.block(0, 0, pbs, cbs) -= coeff * interface_term;
        dr_rhs.block(0, cbs, pbs, cbs) += (coeff-1.0) * interface_term;
    }

    matrix_type oper = dr_lhs.ldlt().solve(dr_rhs);

    return std::make_pair(oper, dr_rhs);
}

///////////////////////   Stokes Stabilization

template<typename T, size_t ET, typename F1>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>
make_stokes_interface_stabilization
(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
 const hho_degree_info& di, const F1& level_set_function)
{
    auto hT = diameter(msh, cl);

    auto celdeg = di.cell_degree();
    auto pdeg = di.face_degree();

    vector_cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, celdeg);
    cell_basis<cuthho_mesh<T, ET>,T>        pb(msh, cl, pdeg);

    auto cbs = cb.size();
    auto pbs = pb.size();

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero( 2*(cbs+pbs), 2*(cbs+pbs));
    Matrix<T, Dynamic, Dynamic> ret_temp = Matrix<T, Dynamic, Dynamic>::Zero( cbs+pbs, cbs+pbs  );

    auto qpsi = integrate_interface(msh, cl, celdeg - 1 + pdeg ,
                                    element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : qpsi)
    {
        const auto v_dphi = cb.eval_gradients(qp.first);
        const auto s_phi  = pb.eval_basis(qp.first);
        const auto n = level_set_function.normal(qp.first);
        const Matrix<T, Dynamic, Dynamic> v_dphi_n = outer_product(v_dphi, n);
        const Matrix<T, Dynamic, Dynamic> s_phi_n = s_phi * n.transpose();

        ret_temp.block(0, 0, cbs, cbs) += qp.second * v_dphi_n * v_dphi_n.transpose();
        ret_temp.block(cbs, 0, pbs, cbs) += qp.second * s_phi_n * v_dphi_n.transpose();
        ret_temp.block(cbs, cbs, pbs, pbs) += qp.second * s_phi_n * s_phi_n.transpose();
    }

    // vel -- vel
    ret.block(0, 0, cbs, cbs)     += hT * ret_temp.block(0, 0, cbs, cbs);
    ret.block(0, cbs, cbs, cbs)   -= hT * ret_temp.block(0, 0, cbs, cbs);
    ret.block(cbs, 0, cbs, cbs)   -= hT * ret_temp.block(0, 0, cbs, cbs);
    ret.block(cbs, cbs, cbs, cbs) += hT * ret_temp.block(0, 0, cbs, cbs);

    // vel -- p
    ret.block(0, 2*cbs, cbs, pbs)       -= hT * ret_temp.block(cbs, 0, pbs, cbs).transpose();
    ret.block(0, 2*cbs + pbs, cbs, pbs) += hT * ret_temp.block(cbs, 0, pbs, cbs).transpose();
    ret.block(cbs, 2*cbs, cbs, pbs)     += hT * ret_temp.block(cbs, 0, pbs, cbs).transpose();
    ret.block(cbs, 2*cbs+pbs, cbs, pbs) -= hT * ret_temp.block(cbs, 0, pbs, cbs).transpose();

    // p -- vel
    ret.block(2*cbs, 0, pbs, cbs)       -= hT * ret_temp.block(cbs, 0, pbs, cbs);
    ret.block(2*cbs, cbs, pbs, cbs)     += hT * ret_temp.block(cbs, 0, pbs, cbs);
    ret.block(2*cbs+pbs, 0, pbs, cbs)   += hT * ret_temp.block(cbs, 0, pbs, cbs);
    ret.block(2*cbs+pbs, cbs, pbs, cbs) -= hT * ret_temp.block(cbs, 0, pbs, cbs);

    // p -- p
    ret.block(2*cbs, 2*cbs, pbs, pbs)         += hT * ret_temp.block(cbs, cbs, pbs, pbs);
    ret.block(2*cbs, 2*cbs+pbs, pbs, pbs)     -= hT * ret_temp.block(cbs, cbs, pbs, pbs);
    ret.block(2*cbs+pbs, 2*cbs, pbs, pbs)     -= hT * ret_temp.block(cbs, cbs, pbs, pbs);
    ret.block(2*cbs+pbs, 2*cbs+pbs, pbs, pbs) += hT * ret_temp.block(cbs, cbs, pbs, pbs);

    return ret;
}

///////////////////////   RHS  pressure
template<typename T, size_t ET, typename F1, typename F2>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, 1>
make_pressure_rhs(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
                  size_t degree, const element_location where,
                  const F1& level_set_function, const F2& bcs_fun)
{
    if( location(msh, cl) != element_location::ON_INTERFACE )
    {
        auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(degree);
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);
        return ret;
    }

    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, degree);
    auto cbs = cb.size();
    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

    auto qpsi = integrate_interface(msh, cl, 2*degree, element_location::IN_NEGATIVE_SIDE );
    for (auto& qp : qpsi)
    {
        auto phi = cb.eval_basis(qp.first);
        auto n = level_set_function.normal(qp.first);

        ret -= qp.second * bcs_fun(qp.first).dot(n) * phi;
    }

    return ret;
}


template<typename T, size_t ET, typename F1, typename F2>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, 1>
make_stokes_interface_stabilization_RHS
(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
 const hho_degree_info& di, const F1& level_set_function, const F2& neumann_jump)
{
    auto hT = diameter(msh, cl);

    auto celdeg = di.cell_degree();
    auto pdeg = di.face_degree();

    vector_cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, celdeg);
    cell_basis<cuthho_mesh<T, ET>,T>        pb(msh, cl, pdeg);

    auto cbs = cb.size();
    auto pbs = pb.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero( 2*(cbs + 2*pbs) );

    auto qpsi = integrate_interface(msh, cl, 2*celdeg, element_location::IN_NEGATIVE_SIDE );
    for (auto& qp : qpsi)
    {
        auto v_dphi = cb.eval_gradients(qp.first);
        auto s_phi  = pb.eval_basis(qp.first);
        auto n = level_set_function.normal(qp.first);
        const Matrix<T, Dynamic, Dynamic> v_dphi_n = outer_product(v_dphi, n);
        const Matrix<T, Dynamic, Dynamic> s_phi_n = s_phi * n.transpose();

        ret.head(cbs) += hT * qp.second * v_dphi_n * neumann_jump(qp.first);
        ret.tail(pbs) -= hT * qp.second * s_phi_n  * neumann_jump(qp.first);
    }

    return ret;
}

///////////////////////   FICTITIOUS DOMAIN METHODS  ///////////////////////////

template<typename T, size_t ET, typename testType>
class stokes_fictdom_method
{
    using Mat  = Matrix<T, Dynamic, Dynamic>;
    using Vect = Matrix<T, Dynamic, 1>;
    using Mesh = cuthho_mesh<T, ET>;

protected:
    stokes_fictdom_method(){}

    virtual std::pair< std::pair<Mat,Mat>, std::pair<Vect,Vect> >
    make_contrib_cut(const Mesh& msh, const typename Mesh::cell_type& cl,
                     const testType test_case, const hho_degree_info hdi,
                     const element_location where = element_location::IN_NEGATIVE_SIDE,
                     const params<T>& parms = params<T>())
    {
    }

public:
    std::pair< std::pair<Mat,Mat>, std::pair<Vect,Vect> >
    make_contrib_uncut(const Mesh& msh, const typename Mesh::cell_type& cl,
                       const hho_degree_info hdi, const testType test_case)
    {
        auto gr = make_hho_gradrec_matrix(msh, cl, hdi);
        Mat stab = make_hho_vector_naive_stabilization(msh, cl, hdi);
        Mat lc = gr.second + stab;
        auto dr = make_hho_divergence_reconstruction(msh, cl, hdi);
        Vect f = make_vector_rhs(msh, cl, hdi.cell_degree(), test_case.rhs_fun);
        Vect p_rhs = Vect::Zero( dr.first.rows() );
        return std::make_pair( std::make_pair(lc, dr.second) , std::make_pair(f,p_rhs) );
    }


    std::pair< std::pair<Mat,Mat>, std::pair<Vect,Vect> >
    make_contrib(const Mesh& msh, const typename Mesh::cell_type& cl,
                 const testType test_case, const hho_degree_info hdi,
                 const element_location where = element_location::IN_NEGATIVE_SIDE,
                 const params<T>& parms = params<T>())
    {
        if( location(msh, cl) == where )
            return make_contrib_uncut(msh, cl, hdi, test_case);
        else if( location(msh, cl) != element_location::ON_INTERFACE )
        {
            Mat lc;
            Vect f;
            return std::make_pair(std::make_pair(lc, lc) , std::make_pair(f,f) );
        }
        else // on interface
            return make_contrib_cut(msh, cl, test_case, hdi, where, parms);
    }
};

/////////////////////////  GRADREC_FICTITIOUS_METHOD

template<typename T, size_t ET, typename testType>
class gradrec_stokes_fictdom_method : public stokes_fictdom_method<T, ET, testType>
{
    using Mat = Matrix<T, Dynamic, Dynamic>;
    using Vect = Matrix<T, Dynamic, 1>;
    using Mesh = cuthho_mesh<T, ET>;

public:
    T eta;

    gradrec_stokes_fictdom_method(T eta_)
        : stokes_fictdom_method<T,ET,testType>(), eta(eta_) {}

    std::pair< std::pair<Mat,Mat>, std::pair<Vect,Vect> >
    make_contrib_cut(const Mesh& msh, const typename Mesh::cell_type& cl,
                     const testType test_case, const hho_degree_info hdi,
                     const element_location where = element_location::IN_NEGATIVE_SIDE,
                     const params<T>& parms = params<T>())
    {
        // LHS
        auto gr = make_hho_gradrec_matrix(msh, cl, test_case.level_set_, hdi, where, 1.0);
        Mat stab = make_hho_vector_cut_stabilization(msh, cl, hdi, where)
            + make_hho_cut_interface_vector_penalty(msh, cl, hdi, eta);
        Mat lc = gr.second + stab;
        auto dr = make_hho_divergence_reconstruction(msh, cl, test_case.level_set_,
                                                     hdi, where, 1.0);

        // RHS
        auto celdeg = hdi.cell_degree();
        auto cbs = vector_cell_basis<Mesh,T>::size(celdeg);

        Vect f = Vect::Zero(lc.rows());
        f.block(0, 0, cbs, 1) += make_vector_rhs(msh, cl, celdeg, test_case.rhs_fun, where);
        f.block(0, 0, cbs, 1) += make_vector_rhs_penalty(msh, cl, celdeg, test_case.bcs_vel, eta);
        f += make_vector_GR_rhs(msh, cl, celdeg, test_case.bcs_vel, test_case.level_set_, gr.first);
        auto p_rhs = make_pressure_rhs(msh, cl, hdi.face_degree(), where,
                                       test_case.level_set_, test_case.bcs_vel);

        return std::make_pair(std::make_pair(lc, dr.second), std::make_pair(f,p_rhs) );
    }
};



template<typename T, size_t ET, typename testType>
auto make_gradrec_stokes_fictdom_method(const cuthho_mesh<T, ET>& msh, const T eta_,
                                 const testType test_case)
{
    return gradrec_stokes_fictdom_method<T, ET, testType>(eta_);
}

///////////////////////////

template<typename Mesh, typename testType>
stokes_test_info<typename Mesh::coordinate_type>
run_cuthho_fictdom(const Mesh& msh, size_t degree, testType test_case)
{
    using RealType = typename Mesh::coordinate_type;

    auto level_set_function = test_case.level_set_;

    auto sol_vel = test_case.sol_vel;
    auto vel_grad = test_case.vel_grad;
    auto bcs_fun = test_case.bcs_vel;


    /************** OPEN SILO DATABASE **************/
    silo_database silo;
    silo.create("cuthho_fictdom.silo");
    silo.add_mesh(msh, "mesh");

    /************** MAKE A SILO VARIABLE FOR CELL POSITIONING **************/
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

    /************** MAKE A SILO VARIABLE FOR LEVEL SET FUNCTION **************/
    std::vector<RealType> level_set_vals;
    for (auto& pt : msh.points)
        level_set_vals.push_back( level_set_function(pt) );
    silo.add_variable("mesh", "level_set", level_set_vals.data(), level_set_vals.size(), nodal_variable_t);

    /************** MAKE A SILO VARIABLE FOR NODE POSITIONING **************/
    std::vector<RealType> node_pos;
    for (auto& n : msh.nodes)
        node_pos.push_back( location(msh, n) == element_location::IN_POSITIVE_SIDE ? +1.0 : -1.0 );
    silo.add_variable("mesh", "node_pos", node_pos.data(), node_pos.size(), nodal_variable_t);


    timecounter tc;

    bool sc = true; // static condensation

    /************** ASSEMBLE PROBLEM **************/
    hho_degree_info hdi(degree+1, degree);

    element_location where = element_location::IN_NEGATIVE_SIDE;

    tc.tic();
    auto assembler = make_stokes_fict_assembler(msh, hdi, where);
    auto assembler_sc = make_stokes_fict_condensed_assembler2(msh, hdi, where);
    // auto assembler_sc = make_stokes_fict_condensed_assembler(msh, hdi, where);

    // method with gradient reconstruction (penalty-free)
    auto class_meth = make_gradrec_stokes_fictdom_method(msh, 1.0, test_case);

    for (auto& cl : msh.cells)
    {
        if( !(location(msh, cl) == element_location::ON_INTERFACE || location(msh, cl) == where) )
            continue;
        auto contrib = class_meth.make_contrib(msh, cl, test_case, hdi,
                                               element_location::IN_NEGATIVE_SIDE);
        auto lc_A = contrib.first.first;
        auto lc_B = -contrib.first.second;
        auto rhs_A = contrib.second.first;
        auto rhs_B = -contrib.second.second;


        if( sc )
            assembler_sc.assemble(msh, cl, lc_A, lc_B, rhs_A, rhs_B, bcs_fun);
        else
            assembler.assemble(msh, cl, lc_A, lc_B, rhs_A, rhs_B, bcs_fun);
    }

    if( sc )
        assembler_sc.finalize();
    else
        assembler.finalize();

    tc.toc();
    std::cout << bold << yellow << "Matrix assembly: " << tc << " seconds" << reset << std::endl;

    if( sc )
        std::cout << "System unknowns: " << assembler_sc.LHS.rows() << std::endl;
    else
        std::cout << "System unknowns: " << assembler.LHS.rows() << std::endl;

    std::cout << "Cells: " << msh.cells.size() << std::endl;
    std::cout << "Faces: " << msh.faces.size() << std::endl;


    /************** SOLVE **************/
    tc.tic();
#if 1
    SparseLU<SparseMatrix<RealType>>  solver;
    Matrix<RealType, Dynamic, 1> sol;

    if( sc )
    {
        solver.analyzePattern(assembler_sc.LHS);
        solver.factorize(assembler_sc.LHS);
        sol = solver.solve(assembler_sc.RHS);
    }
    else
    {
        solver.analyzePattern(assembler.LHS);
        solver.factorize(assembler.LHS);
        sol = solver.solve(assembler.RHS);
    }
#endif
#if 0
    Matrix<RealType, Dynamic, 1> sol;
    cg_params<RealType> cgp;
    cgp.histfile = "cuthho_cg_hist.dat";
    cgp.verbose = true;
    cgp.apply_preconditioner = true;
    if( sc )
    {
        sol = Matrix<RealType, Dynamic, 1>::Zero(assembler_sc.RHS.rows());
        cgp.max_iter = assembler_sc.LHS.rows();
        conjugated_gradient(assembler_sc.LHS, assembler_sc.RHS, sol, cgp);
    }
    else
    {
        sol = Matrix<RealType, Dynamic, 1>::Zero(assembler.RHS.rows());
        cgp.max_iter = assembler.LHS.rows();
        conjugated_gradient(assembler.LHS, assembler.RHS, sol, cgp);
    }
#endif
    tc.toc();
    std::cout << bold << yellow << "Linear solver: " << tc << " seconds" << reset << std::endl;

    /************** POSTPROCESS **************/



    postprocess_output<RealType>  postoutput;

    auto uT1_gp  = std::make_shared< gnuplot_output_object<RealType> >("fictdom_uT1.dat");
    auto uT2_gp  = std::make_shared< gnuplot_output_object<RealType> >("fictdom_uT2.dat");
    auto p_gp    = std::make_shared< gnuplot_output_object<RealType> >("fictdom_p.dat");

    tc.tic();
    RealType    H1_error = 0.0;
    RealType    L2_error = 0.0;
    RealType    L2_pressure_error = 0.0;
    for (auto& cl : msh.cells)
    {
        bool hide_fict_dom = true; // hide the fictitious domain in the gnuplot outputs
        if (hide_fict_dom && location(msh,cl) == element_location::IN_POSITIVE_SIDE)
            continue;

        vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> cb(msh, cl, hdi.cell_degree());
        cell_basis<cuthho_poly_mesh<RealType>, RealType> s_cb(msh, cl, hdi.face_degree());

        auto cbs = cb.size();

        Matrix<RealType, Dynamic, 1> locdata_vel, locdata_p;
        if( sc )
        {
            locdata_vel = assembler_sc.take_velocity(msh, cl, sol, sol_vel);
            locdata_p   = assembler_sc.take_pressure(msh, cl, sol, sol_vel);
        }
        else
        {
            locdata_vel = assembler.take_velocity(msh, cl, sol, sol_vel);
            locdata_p   = assembler.take_pressure(msh, cl, sol);
        }

        Matrix<RealType, Dynamic, 1> cell_v_dofs = locdata_vel.head(cbs);

        auto bar = barycenter(msh, cl, element_location::IN_NEGATIVE_SIDE);

        if ( location(msh, cl) == element_location::IN_NEGATIVE_SIDE ||
             location(msh, cl) == element_location::ON_INTERFACE )
        {
            Matrix<RealType, 1, 2> real_grad_int = Matrix<RealType, 1, 2>::Zero();
            Matrix<RealType, 1, 2> comp_grad_int = Matrix<RealType, 1, 2>::Zero();
            auto qps = integrate(msh, cl, 2*hdi.cell_degree(), element_location::IN_NEGATIVE_SIDE);
            for (auto& qp : qps)
            {
                /* Compute H1-error */
                auto t_dphi = cb.eval_gradients( qp.first );
                Matrix<RealType, 2, 2> grad = Matrix<RealType, 2, 2>::Zero();

                for (size_t i = 0; i < cbs; i++ )
                    grad += cell_v_dofs(i) * t_dphi[i].block(0, 0, 2, 2);

                Matrix<RealType, 2, 2> grad_diff = vel_grad(qp.first) - grad;

                H1_error += qp.second * inner_product(grad_diff , grad_diff);


                /* L2 - error */
                auto t_phi = cb.eval_basis( qp.first );
                auto v = t_phi.transpose() * cell_v_dofs;
                Matrix<RealType, 2, 1> sol_diff = sol_vel(qp.first) - v;
                L2_error += qp.second * sol_diff.dot(sol_diff);

                uT1_gp->add_data( qp.first, v(0) );
                uT2_gp->add_data( qp.first, v(1) );

                /* L2 - pressure - error */
                auto s_cphi = s_cb.eval_basis( qp.first );
                RealType p_num = s_cphi.dot(locdata_p);
                RealType p_diff = test_case.sol_p( qp.first ) - p_num;
                L2_pressure_error += qp.second * p_diff * p_diff;

                p_gp->add_data( qp.first, p_num );
            }
        }
    }

    std::cout << bold << green << "Energy-norm absolute error:           " << std::sqrt(H1_error) << std::endl;
    std::cout << bold << green << "L2 - pressure - error:                " << std::sqrt(L2_pressure_error) << std::endl;

    postoutput.add_object(uT1_gp);
    postoutput.add_object(uT2_gp);
    postoutput.add_object(p_gp);
    postoutput.write();

    stokes_test_info<RealType> TI;
    TI.H1_vel = std::sqrt(H1_error);
    TI.L2_vel = std::sqrt(L2_error);
    TI.L2_p   = std::sqrt(L2_pressure_error);

    tc.toc();
    std::cout << bold << yellow << "Postprocessing: " << tc << " seconds" << reset << std::endl;

    return TI;
}



//////////////////////////////  INTERFACE METHODS  ///////////////////////////

template<typename T, size_t ET, typename testType>
class stokes_interface_method
{
    using Mat  = Matrix<T, Dynamic, Dynamic>;
    using Vect = Matrix<T, Dynamic, 1>;
    using Mesh = cuthho_mesh<T, ET>;

protected:
    stokes_interface_method(){}

    virtual std::pair<Mat, Vect>
    make_contrib_cut(const Mesh& msh, const typename Mesh::cell_type& cl,
                     const testType test_case, const hho_degree_info hdi)
    {
    }

public:
    std::pair<Mat, Vect>
    make_contrib_uncut(const Mesh& msh, const typename Mesh::cell_type& cl,
                       const hho_degree_info hdi, const testType test_case)
    {
        T kappa;
        if ( location(msh, cl) == element_location::IN_NEGATIVE_SIDE )
            kappa = test_case.parms.kappa_1;
        else
            kappa = test_case.parms.kappa_2;

        auto gr = make_hho_gradrec_matrix(msh, cl, hdi);
        Mat stab = make_hho_vector_naive_stabilization(msh, cl, hdi);
        Mat lc = kappa * (gr.second + stab);
        auto dr = make_hho_divergence_reconstruction(msh, cl, hdi);
        Mat f = make_vector_rhs(msh, cl, hdi.cell_degree(), test_case.rhs_fun);

        size_t v_size = gr.second.rows();
        size_t p_size = dr.first.rows();
        size_t loc_size = v_size + p_size;
        Mat lhs = Mat::Zero( loc_size, loc_size );
        Vect rhs = Vect::Zero( loc_size );

        lhs.block(0, 0, v_size, v_size) = lc;
        lhs.block(0, v_size, v_size, p_size) = -dr.second.transpose();
        lhs.block(v_size, 0, p_size, v_size) = -dr.second;

        rhs.head(f.rows()) = f;
        return std::make_pair(lhs, rhs);
    }


    std::pair<Mat, Vect>
    make_contrib(const Mesh& msh, const typename Mesh::cell_type& cl,
                 const testType test_case, const hho_degree_info hdi)
    {
        if( location(msh, cl) != element_location::ON_INTERFACE )
            return make_contrib_uncut(msh, cl, hdi, test_case);
        else // on interface
            return make_contrib_cut(msh, cl, test_case, hdi);
    }
};


////////////////////////  SYMMETRIC GRADREC INTERFACE METHOD


template<typename T, size_t ET, typename testType>
class Sym_gradrec_stokes_interface_method : public stokes_interface_method<T, ET, testType>
{
    using Mat = Matrix<T, Dynamic, Dynamic>;
    using Vect = Matrix<T, Dynamic, 1>;
    using Mesh = cuthho_mesh<T, ET>;

public:
    T eta, gamma_0;

    Sym_gradrec_stokes_interface_method(T eta_, T gamma_)
        : stokes_interface_method<T,ET,testType>(), eta(eta_), gamma_0(gamma_) {}

    std::pair<Mat, Vect>
    make_contrib_cut(const Mesh& msh, const typename Mesh::cell_type& cl,
                     const testType test_case, const hho_degree_info hdi)
    {
        auto parms = test_case.parms;
        auto level_set_function = test_case.level_set_;

        ///////////////   LHS
        auto celdeg = hdi.cell_degree();
        auto pdeg = hdi.face_degree();
        auto cbs = vector_cell_basis<Mesh,T>::size(celdeg);
        auto pbs = cell_basis<Mesh,T>::size(pdeg);

        // GR
        auto gr_n = make_hho_gradrec_matrix_interface(msh, cl, level_set_function, hdi,
                                                      element_location::IN_NEGATIVE_SIDE, 0.5);
        auto gr_p = make_hho_gradrec_matrix_interface(msh, cl, level_set_function, hdi,
                                                      element_location::IN_POSITIVE_SIDE, 0.5);

        // stab
        Mat stab = make_hho_vector_stabilization_interface(msh, cl, level_set_function, hdi,parms);

        Mat penalty = make_hho_cut_interface_vector_penalty(msh, cl, hdi, eta).block(0,0,cbs,cbs);
        stab.block(0, 0, cbs, cbs) += parms.kappa_2 * penalty;
        stab.block(0, cbs, cbs, cbs) -= parms.kappa_2 * penalty;
        stab.block(cbs, 0, cbs, cbs) -= parms.kappa_2 * penalty;
        stab.block(cbs, cbs, cbs, cbs) += parms.kappa_2 * penalty;

        Mat lc = stab + parms.kappa_1 * gr_n.second + parms.kappa_2 * gr_p.second;

        // DR
        auto dr_n = make_hho_divergence_reconstruction_interface
            (msh, cl, level_set_function, hdi, element_location::IN_NEGATIVE_SIDE, 0.5);
        auto dr_p = make_hho_divergence_reconstruction_interface
            (msh, cl, level_set_function, hdi, element_location::IN_POSITIVE_SIDE, 0.5);


        Mat lhs = Mat::Zero(lc.rows() + 2*pbs, lc.rows() + 2*pbs);
        lhs.block(0, 0, lc.rows(), lc.rows()) = lc;
        lhs.block(0, lc.rows(), lc.rows(), pbs) -= dr_n.second.transpose();
        lhs.block(0, lc.rows() + pbs, lc.rows(), pbs) -= dr_p.second.transpose();
        lhs.block(lc.rows(), 0, pbs, lc.rows()) -= dr_n.second;
        lhs.block(lc.rows() + pbs, 0, pbs, lc.rows()) -= dr_p.second;


        // stokes stabilization terms
        auto stokes_stab = make_stokes_interface_stabilization(msh, cl, hdi, level_set_function);
        lhs.block(0, 0, 2*cbs, 2*cbs) -= gamma_0 * stokes_stab.block(0, 0, 2*cbs, 2*cbs);
        lhs.block(0, lc.rows(), 2*cbs, 2*pbs) -= gamma_0 * stokes_stab.block(0,2*cbs,2*cbs,2*pbs);
        lhs.block(lc.rows(), 0, 2*pbs, 2*cbs) -= gamma_0 * stokes_stab.block(2*cbs,0,2*pbs,2*cbs);
        lhs.block(lc.rows(), lc.rows(), 2*pbs, 2*pbs)
            -= gamma_0 * stokes_stab.block(2*cbs, 2*cbs, 2*pbs, 2*pbs);



        ////////////////    RHS

        Vect f = Vect::Zero(lc.rows());
        // neg part
        f.block(0, 0, cbs, 1) += make_vector_rhs(msh, cl, celdeg, test_case.rhs_fun,
                                                 element_location::IN_NEGATIVE_SIDE);
        f.head(cbs) += 0.5*make_vector_flux_jump(msh,cl,celdeg, element_location::IN_NEGATIVE_SIDE,
                                                 test_case.neumann_jump);

        // pos part
        f.block(cbs, 0, cbs, 1) += make_vector_rhs(msh, cl, celdeg, test_case.rhs_fun,
                                                   element_location::IN_POSITIVE_SIDE);
        f.block(cbs, 0, cbs, 1)
            += 0.5 * make_vector_flux_jump(msh, cl, celdeg, element_location::IN_POSITIVE_SIDE,
                                           test_case.neumann_jump);


        Vect rhs = Vect::Zero(lc.rows() + 2*pbs);
        rhs.head(lc.rows()) = f;

        // stokes stabilization rhs
        auto stab_rhs = make_stokes_interface_stabilization_RHS
            (msh, cl, hdi, level_set_function, test_case.neumann_jump);

        rhs.head(2*cbs) -= gamma_0 * stab_rhs.head(2*cbs);
        rhs.tail(2*pbs) -= gamma_0 * stab_rhs.tail(2*pbs);

        return std::make_pair(lhs, rhs);
    }
};

template<typename T, size_t ET, typename testType>
auto make_sym_gradrec_stokes_interface_method(const cuthho_mesh<T, ET>& msh, const T eta_,
                                              const T gamma_, testType test_case)
{
    return Sym_gradrec_stokes_interface_method<T, ET, testType>(eta_, gamma_);
}

///////////////////////////////////////

template<typename Mesh, typename testType, typename meth>
stokes_test_info<typename Mesh::coordinate_type>
run_cuthho_interface(const Mesh& msh, size_t degree, meth method, testType test_case)
{
    using RealType = typename Mesh::coordinate_type;

    auto level_set_function = test_case.level_set_;

    auto rhs_fun = test_case.rhs_fun;
    auto sol_vel = test_case.sol_vel;
    auto sol_p = test_case.sol_p;
    auto vel_grad = test_case.vel_grad;
    auto bcs_vel = test_case.bcs_vel;
    auto neumann_jump = test_case.neumann_jump;
    struct params<RealType> parms = test_case.parms;

    timecounter tc;

    bool sc = true; // static condensation


    /************** ASSEMBLE PROBLEM **************/
    hho_degree_info hdi(degree+1, degree);

    tc.tic();
    auto assembler = make_stokes_interface_assembler(msh, hdi);
    // auto assembler_sc = make_interface_vector_condensed_assembler(msh, hdi);
    for (auto& cl : msh.cells)
    {
        auto contrib = method.make_contrib(msh, cl, test_case, hdi);
        auto lc = contrib.first;
        auto f = contrib.second;

        if (location(msh, cl) != element_location::ON_INTERFACE)
        {
            // if( sc )
            //     assembler_sc.assemble(msh, cl, lc, f, bcs_vel);
            // else
            assembler.assemble(msh, cl, lc, f, bcs_vel);
            // assembler.assemble(msh, cl, lhs_A, lhs_B, rhs_A, bcs_vel, location(msh,cl));
        }
        else
        {
            // if( sc )
            //     assembler_sc.assemble_cut(msh, cl, lc, f);
            // else
            assembler.assemble_cut(msh, cl, lc, f);
            // assembler.assemble_cut(msh, cl, lhs_A, lhs_B, lhs_C, rhs_A, rhs_B);
        }
    }

    // if( sc )
    //     assembler_sc.finalize();
    // else
    assembler.finalize();

    tc.toc();
    std::cout << bold << yellow << "Matrix assembly: " << tc << " seconds" << reset << std::endl;

    // if( sc )
    //     std::cout << "System unknowns: " << assembler_sc.LHS.rows() << std::endl;
    // else
    std::cout << "System unknowns: " << assembler.LHS.rows() << std::endl;

    std::cout << "Cells: " << msh.cells.size() << std::endl;
    std::cout << "Faces: " << msh.faces.size() << std::endl;

    /************** SOLVE **************/
    tc.tic();
#if 1
    SparseLU<SparseMatrix<RealType>>  solver;
    Matrix<RealType, Dynamic, 1> sol;

    // if( sc )
    // {
    //     solver.analyzePattern(assembler_sc.LHS);
    //     solver.factorize(assembler_sc.LHS);
    //     sol = solver.solve(assembler_sc.RHS);
    // }
    // else
    // {
    solver.analyzePattern(assembler.LHS);
    solver.factorize(assembler.LHS);
    sol = solver.solve(assembler.RHS);
    // }
#endif
#if 0
    Matrix<RealType, Dynamic, 1> sol;
    cg_params<RealType> cgp;
    cgp.histfile = "cuthho_cg_hist.dat";
    cgp.verbose = true;
    cgp.apply_preconditioner = true;
    // if( sc )
    // {
    //     sol = Matrix<RealType, Dynamic, 1>::Zero(assembler_sc.RHS.rows());
    //     cgp.max_iter = assembler_sc.LHS.rows();
    //     conjugated_gradient(assembler_sc.LHS, assembler_sc.RHS, sol, cgp);
    // }
    // else
    // {
    sol = Matrix<RealType, Dynamic, 1>::Zero(assembler.RHS.rows());
    cgp.max_iter = assembler.LHS.rows();
    conjugated_gradient(assembler.LHS, assembler.RHS, sol, cgp);
    // }
#endif
    tc.toc();
    std::cout << bold << yellow << "Linear solver: " << tc << " seconds" << reset << std::endl;

    /************** POSTPROCESS **************/


    postprocess_output<RealType>  postoutput;

    auto uT1_gp  = std::make_shared< gnuplot_output_object<RealType> >("interface_uT1.dat");
    auto uT2_gp  = std::make_shared< gnuplot_output_object<RealType> >("interface_uT2.dat");
    auto p_gp    = std::make_shared< gnuplot_output_object<RealType> >("interface_p.dat");

    tc.tic();
    RealType    H1_error = 0.0;
    RealType    L2_error = 0.0;
    RealType    L2_pressure_error = 0.0;
    for (auto& cl : msh.cells)
    {
        vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> cb(msh, cl, hdi.cell_degree());
        cell_basis<cuthho_poly_mesh<RealType>, RealType> pb(msh, cl, hdi.face_degree());
        auto cbs = cb.size();
        auto pbs = pb.size();

        Matrix<RealType, Dynamic, 1> vel_locdata_n, vel_locdata_p, vel_locdata;
        Matrix<RealType, Dynamic, 1> P_locdata_n, P_locdata_p, P_locdata;
        Matrix<RealType, Dynamic, 1> vel_cell_dofs_n, vel_cell_dofs_p, vel_cell_dofs;

        if (location(msh, cl) == element_location::ON_INTERFACE)
        {
            // if( sc )
            // {
            //     vel_locdata_n = assembler_sc.take_local_data(msh, cl, sol, bcs_vel, element_location::IN_NEGATIVE_SIDE);
            //     vel_locdata_p = assembler_sc.take_local_data(msh, cl, sol, bcs_vel, element_location::IN_POSITIVE_SIDE);
            // }
            // else
            // {
            vel_locdata_n = assembler.take_velocity(msh, cl, sol, bcs_vel, element_location::IN_NEGATIVE_SIDE);
            vel_locdata_p = assembler.take_velocity(msh, cl, sol, bcs_vel, element_location::IN_POSITIVE_SIDE);
            P_locdata_n = assembler.take_pressure(msh,cl, sol, element_location::IN_NEGATIVE_SIDE);
            P_locdata_p = assembler.take_pressure(msh,cl, sol, element_location::IN_POSITIVE_SIDE);
            // }

            vel_cell_dofs_n = vel_locdata_n.head(cbs);
            vel_cell_dofs_p = vel_locdata_p.head(cbs);


            auto qps_n = integrate(msh, cl, 2*hdi.cell_degree(), element_location::IN_NEGATIVE_SIDE);
            for (auto& qp : qps_n)
            {
                /* Compute H1-error */
                auto t_dphi = cb.eval_gradients( qp.first );
                Matrix<RealType, 2, 2> grad = Matrix<RealType, 2, 2>::Zero();

                for (size_t i = 1; i < cbs; i++ )
                    grad += vel_cell_dofs_n(i) * t_dphi[i].block(0, 0, 2, 2);

                Matrix<RealType, 2, 2> grad_diff = vel_grad(qp.first) - grad;
                H1_error += qp.second * inner_product(grad_diff , grad_diff);


                /* Compute L2-error */
                auto t_phi = cb.eval_basis( qp.first );
                auto v = t_phi.transpose() * vel_cell_dofs_n;
                Matrix<RealType, 2, 1> sol_diff = sol_vel(qp.first) - v;
                L2_error += qp.second * sol_diff.dot(sol_diff);

                uT1_gp->add_data( qp.first, v(0) );
                uT2_gp->add_data( qp.first, v(1) );

                /* L2 - pressure - error */
                auto p_phi = pb.eval_basis( qp.first );
                RealType p_num = p_phi.dot(P_locdata_n);
                RealType p_diff = test_case.sol_p( qp.first ) - p_num;
                L2_pressure_error += qp.second * p_diff * p_diff;

                p_gp->add_data( qp.first, p_num );
            }

            auto qps_p = integrate(msh, cl, 2*hdi.cell_degree(), element_location::IN_POSITIVE_SIDE);
            for (auto& qp : qps_p)
            {
                /* Compute H1-error */
                auto t_dphi = cb.eval_gradients( qp.first );
                Matrix<RealType, 2, 2> grad = Matrix<RealType, 2, 2>::Zero();

                for (size_t i = 1; i < cbs; i++ )
                    grad += vel_cell_dofs_p(i) * t_dphi[i].block(0, 0, 2, 2);

                Matrix<RealType, 2, 2> grad_diff = vel_grad(qp.first) - grad;
                H1_error += qp.second * inner_product(grad_diff , grad_diff);

                /* Compute L2-error */
                auto t_phi = cb.eval_basis( qp.first );
                auto v = t_phi.transpose() * vel_cell_dofs_p;
                Matrix<RealType, 2, 1> sol_diff = sol_vel(qp.first) - v;
                L2_error += qp.second * sol_diff.dot(sol_diff);

                uT1_gp->add_data( qp.first, v(0) );
                uT2_gp->add_data( qp.first, v(1) );

                /* L2 - pressure - error */
                auto p_phi = pb.eval_basis( qp.first );
                RealType p_num = p_phi.dot(P_locdata_p);
                RealType p_diff = test_case.sol_p( qp.first ) - p_num;
                L2_pressure_error += qp.second * p_diff * p_diff;

                p_gp->add_data( qp.first, p_num );
            }
        }
        else
        {
            // if( sc )
            // {
            //     vel_locdata = assembler_sc.take_local_data(msh, cl, sol, bcs_vel, element_location::IN_POSITIVE_SIDE);
            // }
            // else
            vel_locdata = assembler.take_velocity(msh, cl, sol, bcs_vel, element_location::IN_POSITIVE_SIDE);
            P_locdata = assembler.take_pressure(msh,cl, sol, element_location::IN_POSITIVE_SIDE);
            vel_cell_dofs = vel_locdata.head(cbs);

            auto qps = integrate(msh, cl, 2*hdi.cell_degree());
            for (auto& qp : qps)
            {
                /* Compute H1-error */
                auto t_dphi = cb.eval_gradients( qp.first );
                Matrix<RealType, 2, 2> grad = Matrix<RealType, 2, 2>::Zero();

                for (size_t i = 1; i < cbs; i++ )
                    grad += vel_cell_dofs(i) * t_dphi[i].block(0, 0, 2, 2);

                Matrix<RealType, 2, 2> grad_diff = vel_grad(qp.first) - grad;
                H1_error += qp.second * inner_product(grad_diff , grad_diff);

                /* Compute L2-error */
                auto t_phi = cb.eval_basis( qp.first );
                auto v = t_phi.transpose() * vel_cell_dofs;
                Matrix<RealType, 2, 1> sol_diff = sol_vel(qp.first) - v;
                L2_error += qp.second * sol_diff.dot(sol_diff);

                uT1_gp->add_data( qp.first, v(0) );
                uT2_gp->add_data( qp.first, v(1) );

                /* L2 - pressure - error */
                auto p_phi = pb.eval_basis( qp.first );
                RealType p_num = p_phi.dot(P_locdata);
                RealType p_diff = test_case.sol_p( qp.first ) - p_num;
                L2_pressure_error += qp.second * p_diff * p_diff;

                p_gp->add_data( qp.first, p_num );
            }
        }

    }

    std::cout << bold << green << "Energy-norm absolute error:           " << std::sqrt(H1_error) << std::endl;
    std::cout << bold << green << "L2-norm absolute error:               " << std::sqrt(L2_error) << std::endl;
    std::cout << bold << green << "Pressure L2-norm absolute error:      " << std::sqrt(L2_pressure_error) << std::endl;

    postoutput.add_object(uT1_gp);
    postoutput.add_object(uT2_gp);
    postoutput.add_object(p_gp);
    postoutput.write();



    stokes_test_info<RealType> TI;
    TI.H1_vel = std::sqrt(H1_error);
    TI.L2_vel = std::sqrt(L2_error);
    TI.L2_p = std::sqrt(L2_pressure_error);

    if (false)
    {
        /////////////// compute condition number
        SparseMatrix<RealType> Mat;
        // // Matrix<RealType, Dynamic, Dynamic> Mat;
        // if (sc)
        //     Mat = assembler_sc.LHS;
        // else
        Mat = assembler.LHS;


        RealType sigma_max, sigma_min;

        // Construct matrix operation object using the wrapper class SparseSymMatProd
        Spectra::SparseSymMatProd<RealType> op(Mat);
        // Construct eigen solver object, requesting the largest eigenvalue
        Spectra::SymEigsSolver< RealType, Spectra::LARGEST_MAGN,
                                Spectra::SparseSymMatProd<RealType> > max_eigs(&op, 1, 10);
        max_eigs.init();
        max_eigs.compute();
        if(max_eigs.info() == Spectra::SUCCESSFUL)
            sigma_max = max_eigs.eigenvalues()(0);


        // Construct eigen solver object, requesting the smallest eigenvalue
        Spectra::SymEigsSolver< RealType, Spectra::SMALLEST_MAGN,
                                Spectra::SparseSymMatProd<RealType> > min_eigs(&op, 1, 10);

        min_eigs.init();
        min_eigs.compute();
        if(min_eigs.info() == Spectra::SUCCESSFUL)
            sigma_min = min_eigs.eigenvalues()(0);

        // compute condition number
        RealType cond = sigma_max / sigma_min;
        TI.cond = cond;
        std::cout << "sigma_max = " << sigma_max << "   sigma_min = "
                  << sigma_min << "  cond = " << cond
                  << std::endl;
    }
    else
        TI.cond = 0.0;

    tc.toc();
    std::cout << bold << yellow << "Postprocessing: " << tc << " seconds" << reset << std::endl;


    return TI;
}




/////////////////////////   AUTOMATIC TESTS  //////////////////////////


void convergence_test(void)
{
    using T = double;

    std::vector<size_t> mesh_sizes, pol_orders;

    // meshes
    mesh_sizes.push_back(8);
    mesh_sizes.push_back(16);
    mesh_sizes.push_back(32);
    mesh_sizes.push_back(64);
    mesh_sizes.push_back(128);
    // mesh_sizes.push_back(256);

    // polynomial orders
    pol_orders.push_back(0);
    pol_orders.push_back(1);
    pol_orders.push_back(2);
    pol_orders.push_back(3);


    // export to files ...
    std::vector<std::string> files;
    files.push_back("./output/test_k0.txt");
    files.push_back("./output/test_k1.txt");
    files.push_back("./output/test_k2.txt");
    files.push_back("./output/test_k3.txt");

    for (std::vector<size_t>::iterator it = pol_orders.begin(); it != pol_orders.end(); it++)
    {
        size_t k = *it;

        std::cout << "start tests for k = " << k << std::endl;

        // init the file
        std::ofstream output;
        output.open (files.at(*it), std::ios::in | std::ios::trunc);
        if (!output.is_open())
            throw std::logic_error("file not open");

        // output << "N\th\tH1\tordre1\tL2\tordre2" << std::endl;
        output << "N\th\tH1\tordre1\tL2\tordre2\tLp\tordre3\tcond" << std::endl;

        // convergence tests
        T previous_H1 = 0.0;
        T previous_L2 = 0.0;
        T previous_p = 0.0;
        T previous_h = 0.0;
        for (std::vector<size_t>::iterator it_msh = mesh_sizes.begin();
             it_msh != mesh_sizes.end(); it_msh++)
        {
            size_t N = *it_msh;

            // init mesh (with agglomeration)
            mesh_init_params<T> mip;
            mip.Nx = N;
            mip.Ny = N;
            cuthho_poly_mesh<T> msh(mip);
            size_t int_refsteps = 1;
            T radius = 1.0/3.0;
            auto circle_level_set_function = circle_level_set<T>(radius, 0.5, 0.5);

            // auto level_set_function = flower_level_set<T>(0.31, 0.5, 0.5, 4, 0.04);
            // auto level_set_function = circle_level_set<T>(radius, 0.5, 0.5);
            // auto level_set_function = square_level_set<T>(1.05, -0.05, -0.05, 1.05);
            // auto level_set_function = square_level_set<T>(1.0, -0.0, -0.0, 1.0);
            auto level_set_function = square_level_set<T>(0.76, 0.24, 0.24, 0.76);
            detect_node_position(msh, level_set_function);
            detect_cut_faces(msh, level_set_function);
            if(1)  // AGGLOMERATION
            {
                detect_cut_cells(msh, level_set_function);
                detect_cell_agglo_set(msh, level_set_function);
                make_neighbors_info_cartesian(msh);
                refine_interface(msh, level_set_function, int_refsteps);
                make_agglomeration(msh, level_set_function);
            }
            else  // NODE DISPLACEMENT
            {
                move_nodes(msh, level_set_function);
                detect_cut_faces(msh, level_set_function); //do it again to update intersection points
                detect_cut_cells(msh, level_set_function);
                refine_interface(msh, level_set_function, int_refsteps);
            }

            // compute solution/errors
            stokes_test_info<T> TI;

            // auto TI = run_cuthho_interface(msh, level_set_function, k, 3, test_case);
            if(1) // sin(\pi x) * sin(\pi y)
            {
                auto test_case = make_test_case_stokes_1(msh, level_set_function);
                TI = run_cuthho_fictdom(msh, k, test_case);
            }

            // report info in the file
            T h = 1.0/N;
            if (it_msh == mesh_sizes.begin())
            {
                output << N << "\t" << h << "\t" << TI.H1_vel << "\t" << "."
                       << "\t" << TI.L2_vel << "\t" << "." << "\t" << TI.L2_p
                       << "\t" << "." << "\t" << TI.cond
                       << std::endl;
            }
            else
            {
                T orderH = log(previous_H1 / TI.H1_vel) / log(previous_h / h);
                T orderL = log(previous_L2 / TI.L2_vel) / log(previous_h / h);
                T orderp = log(previous_p / TI.L2_p) / log(previous_h / h);
                output << N << "\t" << h << "\t" << TI.H1_vel << "\t" << orderH
                       << "\t" << TI.L2_vel << "\t" << orderL << "\t" << TI.L2_p
                       << "\t" << orderp << "\t" << TI.cond
                       << std::endl;
            }
            previous_H1 = TI.H1_vel;
            previous_L2 = TI.L2_vel;
            previous_p = TI.L2_p;
            previous_h = h;
        }
        // close the file
        output.close();
    }

    // update the gnuplot curves
    system("gnuplot './output/gnuplot_script_stokes.txt'");

    // update the .pdf file
    system("pdflatex ./output/autom_tests_stokes.tex");

    // open the .pdf file
    system("xdg-open ./autom_tests_stokes.pdf");
}

//////////////////////////     MAIN        ////////////////////////////
#if 0
int main(int argc, char **argv)
{
    convergence_test();
    // tests_stabilization();
    // interface_residus();
    return 1;
}
#endif

#if 1
int main(int argc, char **argv)
{
    using RealType = double;

    size_t degree           = 0;
    size_t int_refsteps     = 4;

    bool dump_debug         = false;
    bool solve_interface    = false;
    bool solve_fictdom      = false;
    bool agglomeration      = false;

    mesh_init_params<RealType> mip;
    mip.Nx = 5;
    mip.Ny = 5;

    /* k <deg>:     method degree
     * M <num>:     number of cells in x direction
     * N <num>:     number of cells in y direction
     * r <num>:     number of interface refinement steps
     *
     * i:           solve interface problem
     * f:           solve fictitious domain problem
     *
     * D:           use node displacement to solve bad cuts (default)
     * A:           use agglomeration to solve bad cuts
     *
     * d:           dump debug data
     */

    int ch;
    while ( (ch = getopt(argc, argv, "k:M:N:r:ifDAd")) != -1 )
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

            case 'r':
                int_refsteps = atoi(optarg);
                break;

            case 'i':
                solve_interface = true;
                break;

            case 'f':
                solve_fictdom = true;
                break;

            case 'D':
                agglomeration = false;
                break;

            case 'A':
                agglomeration = true;
                break;

            case 'd':
                dump_debug = true;
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
    cuthho_poly_mesh<RealType> msh(mip);
    tc.toc();
    std::cout << bold << yellow << "Mesh generation: " << tc << " seconds" << reset << std::endl;
    /************** LEVEL SET FUNCTION **************/
    RealType radius = 1.0/3.0;
    auto level_set_function = circle_level_set<RealType>(radius, 0.5, 0.5);
    // auto level_set_function = line_level_set<RealType>(0.5);
    // auto level_set_function = flower_level_set<RealType>(0.31, 0.5, 0.5, 4, 0.04);
    /************** DO cutHHO MESH PROCESSING **************/

    tc.tic();
    detect_node_position(msh, level_set_function);
    detect_cut_faces(msh, level_set_function);

    if (agglomeration)
    {
        detect_cut_cells(msh, level_set_function);
        detect_cell_agglo_set(msh, level_set_function);
        make_neighbors_info_cartesian(msh);
        // make_neighbors_info(msh);
        refine_interface(msh, level_set_function, int_refsteps);
        make_agglomeration(msh, level_set_function);
    }
    else
    {
        move_nodes(msh, level_set_function);
        detect_cut_faces(msh, level_set_function); //do it again to update intersection points
        detect_cut_cells(msh, level_set_function);
        refine_interface(msh, level_set_function, int_refsteps);
    }


    tc.toc();
    std::cout << bold << yellow << "cutHHO-specific mesh preprocessing: " << tc << " seconds" << reset << std::endl;

    if (dump_debug)
    {
        dump_mesh(msh);
        output_mesh_info(msh, level_set_function);
    }

    output_mesh_info(msh, level_set_function);

    // jumps sin_sin -> exp_cos
    // auto test_case = make_test_case_laplacian_jumps_1(msh, level_set_function);
    // auto test_case = make_test_case_vector_laplacian_sin_sin_exp_cos(msh, level_set_function);
    // auto test_case = make_test_case_vector_laplacian_sin_sin(msh, level_set_function);
    // auto test_case = make_test_case_vector_laplacian_jumps_1(msh, level_set_function);
    auto test_case = make_test_case_stokes_1(msh, level_set_function);

    auto method = make_sym_gradrec_stokes_interface_method(msh, 1.0, 0.0, test_case);

    if (solve_interface)
        run_cuthho_interface(msh, degree, method, test_case);

    if (solve_fictdom)
        run_cuthho_fictdom(msh, degree, test_case);


    return 0;
}
#endif
