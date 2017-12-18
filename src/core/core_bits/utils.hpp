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

#include "bases.hpp"
#include "quadratures.hpp"

template<typename Mesh, typename T = typename Mesh::coordinate_type>
Matrix<T, Dynamic, Dynamic>
make_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
{
    cell_basis<Mesh,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);

    auto qps = integrate(msh, cl, 2*degree);

    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        ret += qp.second * phi * phi.transpose();
    }

    return ret;
}

template<typename Mesh, typename T = typename Mesh::coordinate_type>
Matrix<T, Dynamic, Dynamic>
make_mass_matrix(const Mesh& msh, const typename Mesh::face_type& fc, size_t degree)
{
    face_basis<Mesh,T> fb(msh, fc, degree);
    auto fbs = fb.size();

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);

    auto qps = integrate(msh, fc, 2*degree);

    for (auto& qp : qps)
    {
        auto phi = fb.eval_basis(qp.first);
        ret += qp.second * phi * phi.transpose();
    }

    return ret;
}

template<typename Mesh, typename Function>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_rhs(const Mesh& msh, const typename Mesh::cell_type& cl,
         size_t degree, const Function& f)
{
    using T = typename Mesh::coordinate_type;

    cell_basis<Mesh,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

    auto qps = integrate(msh, cl, 2*degree);

    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        ret += qp.second * phi * f(qp.first);
    }

    return ret;
}

template<typename Mesh, typename Function>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
make_rhs(const Mesh& msh, const typename Mesh::face_type& fc,
         size_t degree, const Function& f)
{
    using T = typename Mesh::coordinate_type;

    face_basis<Mesh,T> fb(msh, fc, degree);
    auto fbs = fb.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(fbs);

    auto qps = integrate(msh, fc, 2*degree);

    for (auto& qp : qps)
    {
        auto phi = fb.eval_basis(qp.first);
        ret += qp.second * phi * f(qp.first);
    }

    return ret;
}

template<typename Mesh, typename Function>
Matrix<typename Mesh::coordinate_type, Dynamic, 1>
project_function(const Mesh& msh, const typename Mesh::cell_type& cl,
                 size_t degree, const Function& f)
{
    using T = typename Mesh::coordinate_type;

    auto cbs = cell_basis<Mesh,T>::size(degree);
    auto fbs = face_basis<Mesh,T>::size(degree);

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs+4*fbs);

    Matrix<T, Dynamic, Dynamic> cell_mm = make_mass_matrix(msh, cl, degree);
    Matrix<T, Dynamic, 1> cell_rhs = make_rhs(msh, cl, degree, f);
    ret.block(0, 0, cbs, 1) = cell_mm.llt().solve(cell_rhs);

    auto fcs = faces(msh, cl);
    for (size_t i = 0; i < 4; i++)
    {
        auto fc = fcs[i];
        Matrix<T, Dynamic, Dynamic> face_mm = make_mass_matrix(msh, fc, degree);
        Matrix<T, Dynamic, 1> face_rhs = make_rhs(msh, fc, degree, f);
        ret.block(cbs+i*fbs, 0, fbs, 1) = face_mm.llt().solve(face_rhs);
    }

    return ret;
}
