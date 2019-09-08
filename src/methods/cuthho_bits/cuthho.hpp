/*
 *       /\        Matteo Cicuttin (C) 2017,2018; Guillaume Delay 2018,2019
 *      /__\       matteo.cicuttin@enpc.fr        guillaume.delay@enpc.fr
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

#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "core/core"
#include "../hho"
#include "cuthho_mesh.hpp"



//////////////////////   MAKE_HHO_LAPLACIAN  ///////////////////////
template<typename T, size_t ET, typename Function>
std::pair<   Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>  >
make_hho_laplacian(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
                   const Function& level_set_function, hho_degree_info di,
                   element_location where)
{

    if ( !is_cut(msh, cl) )
        return make_hho_laplacian(msh, cl, di);

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    cell_basis<cuthho_mesh<T, ET>,T>     cb(msh, cl, recdeg);

    auto rbs = cell_basis<cuthho_mesh<T, ET>,T>::size(recdeg);
    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    Matrix<T, Dynamic, Dynamic> gr_lhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    Matrix<T, Dynamic, Dynamic> gr_rhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs, cbs + num_faces*fbs);

    /* Cell term (cut) */
    auto qps = integrate(msh, cl, 2*recdeg, where);
    for (auto& qp : qps)
    {
        auto dphi = cb.eval_gradients(qp.first);
        stiff += qp.second * dphi * dphi.transpose();
    }

    auto hT = diameter(msh, cl);

    /* Interface term */
    auto iqps = integrate_interface(msh, cl, 2*recdeg, where);
    for (auto& qp : iqps)
    {
        auto phi    = cb.eval_basis(qp.first);
        auto dphi   = cb.eval_gradients(qp.first);
        Matrix<T,2,1> n      = level_set_function.normal(qp.first);

        stiff -= qp.second * phi * (dphi * n).transpose();
        stiff -= qp.second * (dphi * n) * phi.transpose();
        stiff += qp.second * phi * phi.transpose() * cell_eta(msh, cl) / hT;
    }

    gr_lhs.block(0, 0, rbs, rbs) = stiff;
    gr_rhs.block(0, 0, rbs, cbs) = stiff.block(0, 0, rbs, cbs);

    auto ns = normals(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];
        auto n = ns[i];

        // face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg);
        cut_face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg, where);
        /* Terms on faces */
        auto qps = integrate(msh, fc, 2*recdeg, where);
        for (auto& qp : qps)
        {
            auto c_phi = cb.eval_basis(qp.first);
            auto f_phi = fb.eval_basis(qp.first);
            auto r_dphi_tmp = cb.eval_gradients(qp.first);
            auto r_dphi = r_dphi_tmp.block(0, 0, rbs, 2);
            gr_rhs.block(0, cbs+i*fbs, rbs, fbs) += qp.second * (r_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs, cbs) -= qp.second * (r_dphi * n) * c_phi.transpose();
        }
    }

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;
    return std::make_pair(oper, data);
}


template<typename T, size_t ET, typename Function>
std::pair<   Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>  >
make_hho_laplacian_interface(const cuthho_mesh<T, ET>& msh,
    const typename cuthho_mesh<T, ET>::cell_type& cl,
    const Function& level_set_function, hho_degree_info di, const params<T>& parms = params<T>())
{

    if ( !is_cut(msh, cl) )
        throw std::invalid_argument("The cell is not cut");

    auto recdeg = di.reconstruction_degree();
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    cell_basis<cuthho_mesh<T, ET>,T>     cb(msh, cl, recdeg);

    auto rbs = cell_basis<cuthho_mesh<T, ET>,T>::size(recdeg);
    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(2*rbs, 2*rbs);
    Matrix<T, Dynamic, Dynamic> gr_lhs = Matrix<T, Dynamic, Dynamic>::Zero(2*rbs, 2*rbs);
    Matrix<T, Dynamic, Dynamic> gr_rhs = Matrix<T, Dynamic, Dynamic>::Zero(2*rbs, 2*(cbs + num_faces*fbs));

    /* Cell term (cut) */

    auto qps_n = integrate(msh, cl, 2*recdeg, element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : qps_n)
    {
        auto dphi = cb.eval_gradients(qp.first);
        stiff.block(0,0,rbs,rbs) += parms.kappa_1 * qp.second * dphi * dphi.transpose();
    }

    auto qps_p = integrate(msh, cl, 2*recdeg, element_location::IN_POSITIVE_SIDE);
    for (auto& qp : qps_p)
    {
        auto dphi = cb.eval_gradients(qp.first);
        stiff.block(rbs,rbs,rbs,rbs) += parms.kappa_2 * qp.second * dphi * dphi.transpose();
    }

    auto hT = diameter(msh, cl);

    /* Interface term */
    auto iqps = integrate_interface(msh, cl, 2*recdeg, element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : iqps)
    {
        auto phi        = cb.eval_basis(qp.first);
        auto dphi       = cb.eval_gradients(qp.first);
        Matrix<T,2,1> n = level_set_function.normal(qp.first);

        Matrix<T, Dynamic, Dynamic> a = parms.kappa_1 * qp.second * phi * (dphi * n).transpose();
        Matrix<T, Dynamic, Dynamic> b = parms.kappa_1 * qp.second * (dphi * n) * phi.transpose();
        Matrix<T, Dynamic, Dynamic> c = parms.kappa_1 * qp.second * phi * phi.transpose() * cell_eta(msh, cl) / hT;

        stiff.block(  0,   0, rbs, rbs) -= a;
        stiff.block(rbs,   0, rbs, rbs) += a;

        stiff.block(  0,   0, rbs, rbs) -= b;
        stiff.block(  0, rbs, rbs, rbs) += b;

        stiff.block(  0,   0, rbs, rbs) += c;
        stiff.block(  0, rbs, rbs, rbs) -= c;
        stiff.block(rbs,   0, rbs, rbs) -= c;
        stiff.block(rbs, rbs, rbs, rbs) += c;

    }

    gr_lhs = stiff;
    gr_rhs.block(0,   0, 2*rbs, cbs) = stiff.block(0,   0, 2*rbs, cbs);
    gr_rhs.block(0, cbs, 2*rbs, cbs) = stiff.block(0, rbs, 2*rbs, cbs);

    auto ns = normals(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        auto fc = fcs[i];
        auto n = ns[i];

        // face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg);
        cut_face_basis<cuthho_mesh<T, ET>,T> fb_n(msh, fc, facdeg,
                                                  element_location::IN_NEGATIVE_SIDE);
        cut_face_basis<cuthho_mesh<T, ET>,T> fb_p(msh, fc, facdeg,
                                                  element_location::IN_POSITIVE_SIDE);
        /* Terms on faces */
        auto qps_n = integrate(msh, fc, 2*recdeg, element_location::IN_NEGATIVE_SIDE);
        for (auto& qp : qps_n)
        {
            auto c_phi = cb.eval_basis(qp.first);
            auto f_phi = fb_n.eval_basis(qp.first);
            auto r_dphi = cb.eval_gradients(qp.first);

            gr_rhs.block(0, 0, rbs, cbs) -= parms.kappa_1 * qp.second * (r_dphi * n) * c_phi.transpose();
            size_t col_ofs = 2*cbs + i*fbs;
            gr_rhs.block(0, col_ofs, rbs, fbs) += parms.kappa_1 * qp.second * (r_dphi * n) * f_phi.transpose();
        }

        auto qps_p = integrate(msh, fc, 2*recdeg, element_location::IN_POSITIVE_SIDE);
        for (auto& qp : qps_p)
        {
            auto c_phi = cb.eval_basis(qp.first);
            auto f_phi = fb_p.eval_basis(qp.first);
            auto r_dphi = cb.eval_gradients(qp.first);

            gr_rhs.block(rbs, cbs, rbs, cbs) -= parms.kappa_2 * qp.second * (r_dphi * n) * c_phi.transpose();
            size_t col_ofs = 2*cbs + fbs*fcs.size() + i*fbs;
            gr_rhs.block(rbs, col_ofs, rbs, fbs) += parms.kappa_2 * qp.second * (r_dphi * n) * f_phi.transpose();
        }
    }

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.ldlt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}


///////////////////////   STABILIZATION   //////////////////////////



// make_hho_cut_stabilization
// stabilization term on the faces
template<typename T, size_t ET>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>
make_hho_cut_stabilization(const cuthho_mesh<T, ET>& msh,
                           const typename cuthho_mesh<T, ET>::cell_type& cl,
                           const hho_degree_info& di, element_location where)
{
    if ( !is_cut(msh, cl) )
        return make_hho_naive_stabilization(msh, cl, di);

    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+num_faces*fbs, cbs+num_faces*fbs);
    Matrix<T, Dynamic, Dynamic> If = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);

    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, celdeg);

    auto hT = diameter(msh, cl);

    for (size_t i = 0; i < num_faces; i++)
    {
        auto fc = fcs[i];
        // face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg);
        cut_face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg, where);

        Matrix<T, Dynamic, Dynamic> oper = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs+num_faces*fbs);
        Matrix<T, Dynamic, Dynamic> mass = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> trace = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs);

        oper.block(0, cbs+i*fbs, fbs, fbs) = -If;

        auto qps = integrate(msh, fc, facdeg + celdeg, where);
        for (auto& qp : qps)
        {
            auto c_phi = cb.eval_basis(qp.first);
            auto f_phi = fb.eval_basis(qp.first);

            mass += qp.second * f_phi * f_phi.transpose();
            trace += qp.second * f_phi * c_phi.transpose();
        }

        if (qps.size() == 0) /* Avoid to invert a zero matrix */
            continue;


        oper.block(0, 0, fbs, cbs) = mass.ldlt().solve(trace);

        data += oper.transpose() * mass * oper * (1./hT);
    }

    return data;
}




// make_hho_cut_interface_penalty
// eta is the penalty (Nitsche's) parameter
// return eta h_T^{-1} (u_T , v_T)_{Gamma}
template<typename T, size_t ET>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>
make_hho_cut_interface_penalty(const cuthho_mesh<T, ET>& msh,
                               const typename cuthho_mesh<T, ET>::cell_type& cl,
                               const hho_degree_info& di, const T eta)
{
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);

    auto num_faces = faces(msh, cl).size();

    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, celdeg);

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+num_faces*fbs, cbs+num_faces*fbs);

    auto hT = diameter(msh, cl);

    auto iqps = integrate_interface(msh, cl, 2*celdeg, element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : iqps)
    {
        const auto c_phi  = cb.eval_basis(qp.first);

        data.block(0, 0, cbs, cbs) += qp.second * c_phi * c_phi.transpose() * eta / hT;
    }

    return data;
}




//// make_hho_stabilization_interface
// stabilization terms
// method = 1 : negative Nitsche's terms (for bold G -- bold G)
// method = 2 : kappa_2 + no Nitsche's terms (for tilde bold G -- tilde bold G)
// method = 3 : kappa_1 + no Nitsche's terms (for hat bold G -- bold G )
// method = 4 : positive Nitsche's terms (for hat bold G -- hat bold G)
template<typename T, size_t ET, typename Function>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>
make_hho_stabilization_interface(const cuthho_mesh<T, ET>& msh,
                                 const typename cuthho_mesh<T, ET>::cell_type& cl,
                                 const Function& level_set_function,
                                 const hho_degree_info& di, const size_t method,
                                 const params<T>& parms = params<T>())
{
    if ( !is_cut(msh, cl) )
        throw std::invalid_argument("The cell is not cut ...");

    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    Matrix<T, Dynamic, Dynamic> data
        = Matrix<T, Dynamic, Dynamic>::Zero(2*cbs+2*num_faces*fbs, 2*cbs+2*num_faces*fbs);
    // Matrix<T, Dynamic, Dynamic> If = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);

    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, celdeg);

    auto hT = diameter(msh, cl);


    const auto stab_n = make_hho_cut_stabilization(msh, cl, di,element_location::IN_NEGATIVE_SIDE);
    const auto stab_p = make_hho_cut_stabilization(msh, cl, di,element_location::IN_POSITIVE_SIDE);

    // cells--cells
    data.block(0, 0, cbs, cbs) += parms.kappa_1 * stab_n.block(0, 0, cbs, cbs);
    data.block(cbs, cbs, cbs, cbs) += parms.kappa_2 * stab_p.block(0, 0, cbs, cbs);
    // cells--faces
    data.block(0, 2*cbs, cbs, num_faces*fbs)
        += parms.kappa_1 * stab_n.block(0, cbs, cbs, num_faces*fbs);
    data.block(cbs, 2*cbs + num_faces*fbs, cbs, num_faces*fbs)
        += parms.kappa_2 * stab_p.block(0, cbs, cbs, num_faces*fbs);
    // faces--cells
    data.block(2*cbs, 0, num_faces*fbs, cbs)
        += parms.kappa_1 * stab_n.block(cbs, 0, num_faces*fbs, cbs);
    data.block(2*cbs + num_faces*fbs, cbs, num_faces*fbs, cbs)
        += parms.kappa_2 * stab_p.block(cbs, 0, num_faces*fbs, cbs);
    // faces--faces
    data.block(2*cbs, 2*cbs, num_faces*fbs, num_faces*fbs)
        += parms.kappa_1 * stab_n.block(cbs, cbs, num_faces*fbs, num_faces*fbs);
    data.block(2*cbs + num_faces*fbs, 2*cbs + num_faces*fbs, num_faces*fbs, num_faces*fbs)
        += parms.kappa_2 * stab_p.block(cbs, cbs, num_faces*fbs, num_faces*fbs);



    // complementary terms on the interface (cells--cells)
    Matrix<T, Dynamic, Dynamic> term_1 = make_hho_cut_interface_penalty(msh, cl, di, cell_eta(msh, cl)).block(0, 0, cbs, cbs);
    Matrix<T, Dynamic, Dynamic> term_2 = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);

    auto iqps = integrate_interface(msh, cl, 2*celdeg, element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : iqps)
    {
        const auto c_phi  = cb.eval_basis(qp.first);
        const auto dphi   = cb.eval_gradients(qp.first);
        const Matrix<T,2,1> n      = level_set_function.normal(qp.first);

        term_2 += qp.second * c_phi * (dphi * n).transpose();
    }

    if(method == 2)
    {
        data.block(0, cbs, cbs, cbs) -= parms.kappa_2 * term_1;
        data.block(cbs, 0, cbs, cbs) -= parms.kappa_2 * term_1;
        data.block(0, 0, cbs, cbs) += parms.kappa_2 * term_1;
        data.block(cbs, cbs, cbs, cbs) += parms.kappa_2 * term_1;
    }
    else
    {
        data.block(0, 0, cbs, cbs) += parms.kappa_1 * term_1;
        data.block(0, cbs, cbs, cbs) -= parms.kappa_1 * term_1;
        data.block(cbs, 0, cbs, cbs) -= parms.kappa_1 * term_1;
        data.block(cbs, cbs, cbs, cbs) += parms.kappa_1 * term_1;
    }

    if (method == 1)
    {
        data.block(0, 0, cbs, cbs) -= parms.kappa_1 * term_2;
        data.block(0, 0, cbs, cbs) -= parms.kappa_1 * term_2.transpose();
        data.block(cbs, 0, cbs, cbs) += parms.kappa_1 * term_2;
        data.block(0, cbs, cbs, cbs) += parms.kappa_1 * term_2.transpose();
    }
    if (method == 4)
    {
        data.block(0, cbs, cbs, cbs) += parms.kappa_2 * term_2;
        data.block(cbs, 0, cbs, cbs) += parms.kappa_2 * term_2.transpose();
        data.block(cbs, cbs, cbs, cbs) -= parms.kappa_2 * term_2;
        data.block(cbs, cbs, cbs, cbs) -= parms.kappa_2 * term_2.transpose();
    }

    return data;
}


///////////////////////////    GRADREC     //////////////////////////////

template<typename T, size_t ET, typename Function>
std::pair<   Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>  >
make_hho_gradrec_vector(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl, const Function& level_set_function, const hho_degree_info& di, element_location where, const size_t method = 1)
{

    if ( !is_cut(msh, cl) )
        return make_hho_gradrec_vector(msh, cl, di);

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto celdeg  = di.cell_degree();
    const auto facdeg  = di.face_degree();
    const auto graddeg = di.grad_degree();

    cell_basis<cuthho_mesh<T, ET>,T>            cb(msh, cl, celdeg);
    vector_cell_basis<cuthho_mesh<T, ET>,T>     gb(msh, cl, graddeg);


    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);
    auto gbs = vector_cell_basis<cuthho_mesh<T, ET>,T>::size(graddeg);

    const auto num_faces = faces(msh, cl).size();

    matrix_type        gr_lhs = matrix_type::Zero(gbs, gbs);
    matrix_type        gr_rhs = matrix_type::Zero(gbs, cbs + num_faces * fbs);

    const auto qps = integrate(msh, cl, celdeg - 1 + facdeg, where);
    for (auto& qp : qps)
    {
        const auto c_dphi = cb.eval_gradients(qp.first);
        const auto g_phi  = gb.eval_basis(qp.first);

        gr_lhs.block(0, 0, gbs, gbs) += qp.second * g_phi * g_phi.transpose();
        gr_rhs.block(0, 0, gbs, cbs) += qp.second * g_phi * c_dphi.transpose();
    }

    const auto fcs = faces(msh, cl);
    const auto ns = normals(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = ns[i];
        // face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg);
        cut_face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg, where);

        const auto qps_f = integrate(msh, fc, facdeg + std::max(facdeg, celdeg), where);
        for (auto& qp : qps_f)
        {
            const vector_type c_phi      = cb.eval_basis(qp.first);
            const vector_type f_phi      = fb.eval_basis(qp.first);
            const auto        g_phi      = gb.eval_basis(qp.first);
            const vector_type qp_g_phi_n = qp.second * g_phi * n;

            gr_rhs.block(0, cbs + i * fbs, gbs, fbs) += qp_g_phi_n * f_phi.transpose();
            gr_rhs.block(0, 0, gbs, cbs) -= qp_g_phi_n * c_phi.transpose();
        }
    }


    if(method == 1)
    {
        const auto iqps = integrate_interface(msh, cl, celdeg + graddeg, element_location::IN_NEGATIVE_SIDE);
        for (auto& qp : iqps)
        {
            const auto c_phi        = cb.eval_basis(qp.first);
            const auto g_phi        = gb.eval_basis(qp.first);

            Matrix<T,2,1> n = level_set_function.normal(qp.first);
            const vector_type qp_g_phi_n = qp.second * g_phi * n;

            gr_rhs.block(0 , 0, gbs, cbs) -= qp_g_phi_n * c_phi.transpose();
        }
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}




//// make_hho_gradrec_vector_interface
// return the gradient reconstruction for the interface pb
// method = 1 -> no terms on the interface (bold G)
// method = 2 -> terms on the interface scaled by 0.5 (tilde bold G)
// method = 3 -> terms on the interface scaled by 1 (hat bold G)
template<typename T, size_t ET, typename Function>
std::pair<   Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>  >
make_hho_gradrec_vector_interface(const cuthho_mesh<T, ET>& msh,
                                  const typename cuthho_mesh<T, ET>::cell_type& cl,
                                  const Function& level_set_function, const hho_degree_info& di,
                                  element_location where, size_t method)
{

    if ( !is_cut(msh, cl) )
        throw std::invalid_argument("The cell is not cut");

    typedef Matrix<T, Dynamic, Dynamic> matrix_type;
    typedef Matrix<T, Dynamic, 1>       vector_type;

    const auto celdeg  = di.cell_degree();
    const auto facdeg  = di.face_degree();
    const auto graddeg = di.grad_degree();

    cell_basis<cuthho_mesh<T, ET>,T>            cb(msh, cl, celdeg);
    vector_cell_basis<cuthho_mesh<T, ET>,T>     gb(msh, cl, graddeg);


    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(celdeg);
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);
    auto gbs = vector_cell_basis<cuthho_mesh<T, ET>,T>::size(graddeg);

    const auto num_faces = faces(msh, cl).size();

    matrix_type       rhs_tmp = matrix_type::Zero(gbs, cbs + num_faces * fbs);
    matrix_type        gr_lhs = matrix_type::Zero(gbs, gbs);
    matrix_type        gr_rhs = matrix_type::Zero(gbs, 2*cbs + 2*num_faces * fbs);

    const auto qps = integrate(msh, cl, celdeg - 1 + facdeg, where);
    for (auto& qp : qps)
    {
        const auto c_dphi = cb.eval_gradients(qp.first);
        const auto g_phi  = gb.eval_basis(qp.first);

        gr_lhs.block(0, 0, gbs, gbs) += qp.second * g_phi * g_phi.transpose();
        rhs_tmp.block(0, 0, gbs, cbs) += qp.second * g_phi * c_dphi.transpose();
    }

    const auto fcs = faces(msh, cl);
    const auto ns = normals(msh, cl);
    for (size_t i = 0; i < fcs.size(); i++)
    {
        const auto fc = fcs[i];
        const auto n  = ns[i];
        // face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg);
        cut_face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, facdeg, where);


        const auto qps_f = integrate(msh, fc, facdeg + std::max(facdeg, celdeg), where);
        for (auto& qp : qps_f)
        {
            const vector_type c_phi      = cb.eval_basis(qp.first);
            const vector_type f_phi      = fb.eval_basis(qp.first);
            const auto        g_phi      = gb.eval_basis(qp.first);
            const vector_type qp_g_phi_n = qp.second * g_phi * n;

            rhs_tmp.block(0, cbs + i * fbs, gbs, fbs) += qp_g_phi_n * f_phi.transpose();
            rhs_tmp.block(0, 0, gbs, cbs) -= qp_g_phi_n * c_phi.transpose();
        }
    }

    // term on the interface
    const auto iqps = integrate_interface(msh, cl, celdeg + graddeg, element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : iqps)
    {
        const auto c_phi        = cb.eval_basis(qp.first);
        const auto g_phi        = gb.eval_basis(qp.first);

        Matrix<T,2,1> n = level_set_function.normal(qp.first);
        const vector_type qp_g_phi_n = qp.second * g_phi * n;

        gr_rhs.block(0 , 0, gbs, cbs) -= qp_g_phi_n * c_phi.transpose();
        gr_rhs.block(0 , cbs, gbs, cbs) += qp_g_phi_n * c_phi.transpose();
    }
    if (method == 1) gr_rhs = 0.0 * gr_rhs;
    else if (method == 2) gr_rhs = 0.5 * gr_rhs;
    else if (method == 3) {}
    else throw std::invalid_argument("Method should be 1, 2 or 3");

    // other terms
    if(where == element_location::IN_NEGATIVE_SIDE)
    {
        gr_rhs.block(0, 0, gbs, cbs) += rhs_tmp.block(0, 0, gbs, cbs);
        gr_rhs.block(0, 2*cbs, gbs, num_faces*fbs)
            += rhs_tmp.block(0, cbs, gbs, num_faces*fbs);
    }
    else if( where == element_location::IN_POSITIVE_SIDE)
    {
        gr_rhs.block(0, cbs, gbs, cbs) += rhs_tmp.block(0, 0, gbs, cbs);
        gr_rhs.block(0, 2*cbs + num_faces*fbs, gbs, num_faces*fbs)
                     += rhs_tmp.block(0, cbs, gbs, num_faces*fbs);
    }

    matrix_type oper = gr_lhs.ldlt().solve(gr_rhs);
    matrix_type data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}



///////////////////////////  NITSCHE TERMS  /////////////////////////////


template<typename T, size_t ET, typename Function>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, Dynamic>
make_Nitsche(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
             const Function& level_set_function, hho_degree_info di)
{
    auto celdeg = di.cell_degree();
    auto facdeg = di.face_degree();

    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, celdeg);
    auto cbs = cb.size();
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(facdeg);

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    auto size_tot = cbs + num_faces * fbs;

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(size_tot, size_tot);

    if( !is_cut(msh, cl) )
        return ret;

    auto iqps = integrate_interface(msh, cl, 2*celdeg-1, element_location::IN_NEGATIVE_SIDE);
    for (auto& qp : iqps)
    {
        const auto c_phi  = cb.eval_basis(qp.first);
        const auto c_dphi  = cb.eval_gradients(qp.first);
        const Matrix<T,2,1> n = level_set_function.normal(qp.first);
        const auto c_dphi_n = c_dphi * n;

        ret.block(0, 0, cbs, cbs) -= qp.second * c_phi * c_dphi_n.transpose();
        ret.block(0, 0, cbs, cbs) -= qp.second * c_dphi_n * c_phi.transpose();
    }

    return ret;
}




///////////////////////////   RHS VECTOR   //////////////////////////////

// make_rhs for the file cuthho_square.cpp
template<typename T, size_t ET, typename F1, typename F2, typename F3>
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

        auto hT = diameter(msh, cl);

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

        auto qps = integrate(msh, cl, 2*degree, where);
        for (auto& qp : qps)
        {
            auto phi = cb.eval_basis(qp.first);
            ret += qp.second * phi * f(qp.first);
        }


        auto qpsi = integrate_interface(msh, cl, degree, where);
        for (auto& qp : qpsi)
        {
            auto phi = cb.eval_basis(qp.first);
            auto dphi = cb.eval_gradients(qp.first);
            auto n = level_set_function.normal(qp.first);

            ret += qp.second * bcs(qp.first) * ( phi * cell_eta(msh, cl)/hT - dphi*n);
        }


        return ret;
    }
    else
    {
        auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(degree);
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);
        return ret;
    }
}



// make_rhs for the file cuthho_meth2_square.cpp
template<typename T, size_t ET, typename F1, typename F2, typename F3>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, 1>
make_rhs(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
         size_t degree, const F1& f, const element_location where, const F2& level_set_function, const F3& bcs, Matrix<T, Dynamic, Dynamic> GR, size_t method = 1)
{
    if ( location(msh, cl) == where )
        return make_rhs(msh, cl, degree, f);
    else if ( location(msh, cl) == element_location::ON_INTERFACE )
    {
        cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, degree);
        auto cbs = cb.size();

        vector_cell_basis<cuthho_mesh<T, ET>,T> gb(msh, cl, degree-1);
        auto gbs = gb.size();


        auto hT = diameter(msh, cl);

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(GR.cols());
        Matrix<T, Dynamic, 1> source_vect = Matrix<T, Dynamic, 1>::Zero(gbs);
        Matrix<T, Dynamic, 1> grad_term = Matrix<T, Dynamic, 1>::Zero(GR.cols());

        auto qps = integrate(msh, cl, 2*degree, where);
        for (auto& qp : qps)
        {
            auto phi = cb.eval_basis(qp.first);
            ret.block(0, 0, cbs, 1) += qp.second * phi * f(qp.first);
        }


        auto qpsi = integrate_interface(msh, cl, 2*degree, element_location::IN_NEGATIVE_SIDE);
        for (auto& qp : qpsi)
        {
            const auto phi = cb.eval_basis(qp.first);

            ret.block(0, 0, cbs, 1)
                += qp.second * bcs(qp.first) * phi * cell_eta(msh, cl)/hT;
        }

        auto qpsi2 = integrate_interface(msh, cl, 2*degree - 1, element_location::IN_NEGATIVE_SIDE);
        if(method == 1)
        {
            for (auto& qp : qpsi2)
            {
                const auto n = level_set_function.normal(qp.first);
                const auto g_phi  = gb.eval_basis(qp.first);

                source_vect += qp.second * bcs(qp.first) * g_phi * n;
            }

            grad_term += GR.transpose() * source_vect;

            ret -= grad_term;
        }
        else if (method == 2)
        {
            for (auto& qp : qpsi2)
            {
                const auto n = level_set_function.normal(qp.first);
                const auto c_dphi  = cb.eval_gradients(qp.first);
                const auto c_dphi_n  = c_dphi * n;

                ret.block(0, 0, cbs, 1) -= qp.second * bcs(qp.first) * c_dphi_n;
            }
        }
        return ret;
    }
    else
    {
        auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(degree);
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);
        return ret;
    }
}





// make_Dirichlet_jump -> for the file cuthho_square.cpp
template<typename T, size_t ET, typename F1, typename F2>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, 1>
make_Dirichlet_jump(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
                   size_t degree, const element_location where, const F1& level_set_function, 
                   const F2& dir_jump)
{
    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, degree);
    auto cbs = cb.size();
    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

    if( location(msh, cl) != element_location::ON_INTERFACE )
        return ret;

    auto hT = diameter(msh, cl);

    if(where == element_location::IN_NEGATIVE_SIDE) {
        auto qpsi = integrate_interface(msh, cl, 2*degree, element_location::IN_NEGATIVE_SIDE );

        for (auto& qp : qpsi)
        {
            auto phi = cb.eval_basis(qp.first);
            auto dphi = cb.eval_gradients(qp.first);
            auto n = level_set_function.normal(qp.first);

            ret += qp.second * dir_jump(qp.first) * ( phi * cell_eta(msh, cl)/hT - dphi*n);
        }
    }
    else if(where == element_location::IN_POSITIVE_SIDE) {
        auto qpsi = integrate_interface(msh, cl, 2*degree, element_location::IN_NEGATIVE_SIDE );

        for (auto& qp : qpsi)
        {
            auto phi = cb.eval_basis(qp.first);
            ret -= qp.second * dir_jump(qp.first) * phi * cell_eta(msh, cl)/hT;
        }
    }
    return ret;
}



//// make_Dirichlet_jump -> for the file cuthho_meth2_square.cpp
// method = 1 : Nitsche in negative part (for bold G -- bold G)
// method = 2 : kappa_2 + no Nitsche's term (for tilde bold G -- tilde bold G)
// method = 3 : kappa_1 + no Nitsche's term (for hat bold G -- bold G)
// method = 4 : Nitsche in positive part (for hat bold G -- hat bold G)
template<typename T, size_t ET, typename F1, typename F2>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, 1>
make_Dirichlet_jump(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
                    size_t degree, const element_location where, const F1& level_set_function,
                    const F2& dir_jump, const size_t method, const params<T>& parms = params<T>())
{
    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, degree);
    auto cbs = cb.size();
    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

    if( location(msh, cl) != element_location::ON_INTERFACE )
        return ret;

    Matrix<T, Dynamic, 1> term1 = Matrix<T, Dynamic, 1>::Zero(cbs);

    auto hT = diameter(msh, cl);

    auto qpsi = integrate_interface(msh, cl, 2*degree, element_location::IN_NEGATIVE_SIDE );

    for (auto& qp : qpsi)
    {
        auto phi = cb.eval_basis(qp.first);

        term1 += qp.second * dir_jump(qp.first) * phi * cell_eta(msh, cl)/hT;
    }


    if (method == 2 && where == element_location::IN_NEGATIVE_SIDE)
        return parms.kappa_2 * term1;

    if (method == 2 && where == element_location::IN_POSITIVE_SIDE)
        return - parms.kappa_2 * term1;

    if (method == 3 && where == element_location::IN_NEGATIVE_SIDE)
        return parms.kappa_1 * term1;

    if (method == 3 && where == element_location::IN_POSITIVE_SIDE)
        return - parms.kappa_1 * term1;

    Matrix<T, Dynamic, 1> term2 = Matrix<T, Dynamic, 1>::Zero(cbs);
    for (auto& qp : qpsi)
    {
        auto phi = cb.eval_basis(qp.first);
        auto dphi = cb.eval_gradients(qp.first);
        auto n = level_set_function.normal(qp.first);

        term2 += qp.second * dir_jump(qp.first) * dphi * n ;
    }


    if (method == 1 && where == element_location::IN_NEGATIVE_SIDE)
        return parms.kappa_1 * ( term1 - term2 );

    if (method == 1 && where == element_location::IN_POSITIVE_SIDE)
        return - parms.kappa_1 * term1;

    if (method == 4 && where == element_location::IN_NEGATIVE_SIDE)
        return parms.kappa_1 * term1;

    if (method == 4 && where == element_location::IN_POSITIVE_SIDE)
        return parms.kappa_2 * term2 - parms.kappa_1 * term1;
}



template<typename T, size_t ET, typename F1>
Matrix<typename cuthho_mesh<T, ET>::coordinate_type, Dynamic, 1>
make_flux_jump(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
               size_t degree, const element_location where, const F1& flux_jump)
{
    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, degree);
    auto cbs = cb.size();
    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

    if( location(msh, cl) != element_location::ON_INTERFACE )
        return ret;

    auto qpsi = integrate_interface(msh, cl, 2*degree, element_location::IN_NEGATIVE_SIDE);

    for (auto& qp : qpsi)
    {
        auto phi = cb.eval_basis(qp.first);
        ret += qp.second * flux_jump(qp.first) * phi;
    }

    return ret;
}
