/*
 *       /\        Matteo Cicuttin (C) 2017,2018
 *      /__\       matteo.cicuttin@enpc.fr
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

template<typename T, size_t ET>
Matrix<T, Dynamic, Dynamic>
make_mass_matrix(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
                 size_t degree, element_location where)
{
    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);

    auto qps = integrate(msh, cl, 2*degree, where);
    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        ret += qp.second * phi * phi.transpose();
    }

    return ret;
}

template<typename T, size_t ET>
Matrix<T, Dynamic, Dynamic>
make_mass_matrix(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::face_type& fc,
                 size_t degree, element_location where)
{
    face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, degree);
    auto fbs = fb.size();

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(fbs, fbs);

    auto qps = integrate(msh, fc, 2*degree, where);
    for (auto& qp : qps)
    {
        auto phi = fb.eval_basis(qp.first);
        ret += qp.second * phi * phi.transpose();
    }

    return ret;
}

template<typename T, size_t ET, typename Function>
Matrix<T, Dynamic, 1>
make_rhs(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
         size_t degree, element_location where, const Function& f)
{
    cell_basis<cuthho_mesh<T, ET>,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);

    auto qps = integrate(msh, cl, 2*degree, where);

    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        ret += qp.second * phi * f(qp.first);
    }

    return ret;
}

template<typename T, size_t ET, typename Function>
Matrix<T, Dynamic, 1>
make_rhs(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::face_type& fc,
         size_t degree, element_location where, const Function& f)
{
    face_basis<cuthho_mesh<T, ET>,T> fb(msh, fc, degree);
    auto fbs = fb.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(fbs);

    auto qps = integrate(msh, fc, 2*degree, where);

    for (auto& qp : qps)
    {
        auto phi = fb.eval_basis(qp.first);
        ret += qp.second * phi * f(qp.first);
    }

    return ret;
}

template<typename T, size_t ET, typename Function>
Matrix<T, Dynamic, 1>
project_function(const cuthho_mesh<T, ET>& msh, const typename cuthho_mesh<T, ET>::cell_type& cl,
                 hho_degree_info hdi, element_location where, const Function& f)
{
    auto cbs = cell_basis<cuthho_mesh<T, ET>,T>::size(hdi.cell_degree());
    auto fbs = face_basis<cuthho_mesh<T, ET>,T>::size(hdi.face_degree());

    auto fcs = faces(msh, cl);
    auto num_faces = fcs.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs+num_faces*fbs);

    if ( location(msh, cl) != element_location::ON_INTERFACE &&
         location(msh, cl) != where )
        return ret;

    Matrix<T, Dynamic, Dynamic> cell_mm = make_mass_matrix(msh, cl, hdi.cell_degree(), where);
    Matrix<T, Dynamic, 1> cell_rhs = make_rhs(msh, cl, hdi.cell_degree(), where, f);
    ret.block(0, 0, cbs, 1) = cell_mm.llt().solve(cell_rhs);

    for (size_t i = 0; i < num_faces; i++)
    {
        auto fc = fcs[i];

        if ( location(msh, fc) != element_location::ON_INTERFACE &&
             location(msh, fc) != where )
        {
            ret.block(cbs+i*fbs, 0, fbs, 1) = Matrix<T, Dynamic, 1>::Zero(fbs);
        }
        else
        {
            Matrix<T, Dynamic, Dynamic> face_mm = make_mass_matrix(msh, fc, hdi.face_degree(), where);
            Matrix<T, Dynamic, 1> face_rhs = make_rhs(msh, fc, hdi.face_degree(), where, f);
            ret.block(cbs+i*fbs, 0, fbs, 1) = face_mm.llt().solve(face_rhs);
        }
    }

    return ret;
}
