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

#pragma once

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "basic_geom.hpp"
#include "bases.hpp"
#include "quadratures.hpp"

template<typename Mesh>
std::pair<   Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>,
             Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>  >
make_hho_laplacian(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
{
    using T = typename Mesh::coordinate_type;

    cell_basis<Mesh,T>     cb(msh, cl, degree+1);

    auto rbs = cell_basis<Mesh,T>::size(degree+1);
    auto cbs = cell_basis<Mesh,T>::size(degree);
    auto fbs = face_basis<Mesh,T>::size(degree);

    Matrix<T, Dynamic, Dynamic> stiff = Matrix<T, Dynamic, Dynamic>::Zero(rbs, rbs);
    Matrix<T, Dynamic, Dynamic> gr_lhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, rbs-1);
    Matrix<T, Dynamic, Dynamic> gr_rhs = Matrix<T, Dynamic, Dynamic>::Zero(rbs-1, cbs + 4*fbs);

    auto qps = integrate(msh, cl, 2*degree+2);

    for (auto& qp : qps)
    {
        auto dphi = cb.eval_gradients(qp.first);
        stiff += qp.second * dphi * dphi.transpose();
    }

    gr_lhs = stiff.block(1, 1, rbs-1, rbs-1);
    gr_rhs.block(0, 0, rbs-1, cbs) = stiff.block(1, 0, rbs-1, cbs);

    auto fcs = faces(msh, cl);
    auto ns = normals(msh, cl);

    for (size_t i = 0; i < 4; i++)
    {
        auto fc = fcs[i];
        auto n = ns[i];
        face_basis<Mesh,T> fb(msh, fc, degree);
        auto qps_f = integrate(msh, fc, 2*degree);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 1> c_phi_tmp = cb.eval_basis(qp.first);
            Matrix<T, Dynamic, 1> c_phi = c_phi_tmp.head(cbs);
            Matrix<T, Dynamic, 2> c_dphi_tmp = cb.eval_gradients(qp.first);
            Matrix<T, Dynamic, 2> c_dphi = c_dphi_tmp.block(1, 0, rbs-1, 2);
            Matrix<T, Dynamic, 1> f_phi = fb.eval_basis(qp.first);
            gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.second * (c_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs-1, cbs) -= qp.second * (c_dphi * n) * c_phi.transpose();
        }
    }

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}


template<typename Mesh>
Matrix<typename Mesh::coordinate_type, Dynamic, Dynamic>
make_hho_naive_stabilization(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
{
    using T = typename Mesh::coordinate_type;

    auto cbs = cell_basis<Mesh,T>::size(degree);
    auto fbs = face_basis<Mesh,T>::size(degree);

    auto fcs = faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+4*fbs, cbs+4*fbs);
    Matrix<T, Dynamic, Dynamic> If = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);

    cell_basis<Mesh,T> cb(msh, cl, degree);

    for (size_t i = 0; i < 4; i++)
    {
        auto fc = fcs[i];
        face_basis<Mesh,T> fb(msh, fc, degree);

        Matrix<T, Dynamic, Dynamic> oper = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs+4*fbs);
        Matrix<T, Dynamic, Dynamic> mass = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);
        Matrix<T, Dynamic, Dynamic> trace = Matrix<T, Dynamic, Dynamic>::Zero(fbs, cbs);

        oper.block(0, cbs+i*fbs, fbs, fbs) = -If;

        auto qps = integrate(msh, fc, 2*degree);
        for (auto& qp : qps)
        {
            auto c_phi = cb.eval_basis(qp.first);
            auto f_phi = fb.eval_basis(qp.first);

            mass += qp.second * f_phi * f_phi.transpose();
            trace += qp.second * f_phi * c_phi.transpose();
        }

        oper.block(0, 0, fbs, cbs) = mass.llt().solve(trace);

        /* Don't divide by h with this stabilization. It breaks convergence! */
        data += oper.transpose() * mass * oper;
    }

    return data;
}
