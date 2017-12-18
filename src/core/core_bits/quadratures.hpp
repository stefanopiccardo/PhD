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

#include <iostream>
#include <vector>
#include <cassert>

#include "point.hpp"

template<typename T>
std::vector<std::pair<point<T,1>, T>>
edge_quadrature(size_t doe)
{
    std::vector<std::pair<point<T,1>, T>> ret;

    using namespace Eigen;


    if (doe%2 == 0)
        doe++;

    size_t num_nodes = (doe+1)/2;

    //size_t num_nodes = doe+1;

    if (num_nodes == 1)
    {
        auto qp = std::make_pair(point<T,1>({0.5}), 1.0);
        ret.push_back(qp);
        return ret;
    }

    Matrix<T, Dynamic, Dynamic> M = Matrix<T, Dynamic, Dynamic>::Zero(num_nodes, num_nodes);
    for (size_t i = 1; i < num_nodes; i++)
    {
        T p = 4.0 - 1.0/(i*i);
        M(i, i-1) = sqrt(1.0/p);
    }

    SelfAdjointEigenSolver<Matrix<T, Dynamic, Dynamic>> solver;
    solver.compute(M);

    Matrix<T, Dynamic, Dynamic> weights = solver.eigenvectors().row(0);
                                weights = weights.array().square();
    Matrix<T, Dynamic, Dynamic> nodes = solver.eigenvalues();

    assert( weights.size() == nodes.size() );

    for (size_t i = 0; i < nodes.size(); i++)
    {
        auto qp = std::make_pair(point<T,1>({(nodes(i) + 1.)/2.}), weights(i));
        ret.push_back(qp);
    }

    return ret;
}



template<typename Mesh>
std::vector<std::pair<point<typename Mesh::coordinate_type,2>, typename Mesh::coordinate_type>>
integrate(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
{
    using T = typename Mesh::coordinate_type;

    typedef typename Mesh::point_type    point_type;

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

template<typename Mesh>
std::vector<std::pair<point<typename Mesh::coordinate_type,2>, typename Mesh::coordinate_type>>
integrate(const Mesh& msh, const typename Mesh::face_type& fc, size_t degree)
{
    using T = typename Mesh::coordinate_type;

    typedef typename Mesh::point_type    point_type;

    auto qps = edge_quadrature<T>(degree);
    auto pts = points(msh, fc);

    auto scale = (pts[1] - pts[0]);
    auto meas = scale.to_vector().norm();

    std::vector<std::pair<point<T,2>, T>>   ret;

    for (auto itor = qps.begin(); itor != qps.end(); itor++)
    {
        auto qp = *itor;
        auto p = qp.first.x() * scale + pts[0];
        auto w = qp.second * meas;

        ret.push_back( std::make_pair(p, w) );
    }

    return ret;
}


template<typename T, size_t N>
std::ostream&
operator<<(std::ostream& os, const std::pair<point<T,N>, T>& pw_pair)
{
    os << "[ " << pw_pair.first << ", " << pw_pair.second << " ]";
    return os;
}
