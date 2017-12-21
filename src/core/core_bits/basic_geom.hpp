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
#include <stdexcept>

#include "point.hpp"

template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    auto itor = std::lower_bound(msh.cells.begin(), msh.cells.end(), cl);
    if ( itor == msh.cells.end() )
        throw std::logic_error("Cell not found: this is likely a bug.");

    return std::distance(msh.cells.begin(), itor);
}

template<typename Mesh>
size_t
offset(const Mesh& msh, const typename Mesh::face_type& fc)
{
    auto itor = std::lower_bound(msh.faces.begin(), msh.faces.end(), fc);
    if ( itor == msh.faces.end() )
        throw std::logic_error("Face not found: this is likely a bug.");

    return std::distance(msh.faces.begin(), itor);
}

template<typename Mesh>
std::array< typename Mesh::node_type, 4 >
nodes(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    std::array< typename Mesh::node_type, 4 > ret;

    auto ptid2node = [&](size_t ptid) -> auto {
        return msh.nodes.at(ptid);
    };

    std::transform( cl.ptids.begin(), cl.ptids.end(), ret.begin(), ptid2node );

    return ret;
}

template<typename Mesh>
std::array< typename Mesh::point_type, 4 >
points(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    std::array< typename Mesh::point_type, 4 > ret;

    auto ptid2pt = [&](size_t ptid) -> auto {
        return msh.points.at(ptid);
    };

    std::transform( cl.ptids.begin(), cl.ptids.end(), ret.begin(), ptid2pt );

    return ret;
}

template<typename Mesh>
std::array< typename Mesh::point_type, 2 >
points(const Mesh& msh, const typename Mesh::face_type& fc)
{
    std::array< typename Mesh::point_type, 2 > ret;

    auto ptid2pt = [&](size_t ptid) -> auto {
        return msh.points.at(ptid);
    };

    std::transform( fc.ptids.begin(), fc.ptids.end(), ret.begin(), ptid2pt );

    return ret;
}

template<typename Mesh>
typename Mesh::point_type
points(const Mesh& msh, const typename Mesh::node_type& node)
{
    return msh.points.at( node.ptid );
}

template<typename Mesh>
std::array< typename Mesh::face_type, 4 >
faces(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    typedef typename Mesh::face_type face_type;
    std::array< typename Mesh::face_type, 4 > ret;

    face_type f;

    f.ptids[0] = cl.ptids[0];
    f.ptids[1] = cl.ptids[1];
    auto itor = std::lower_bound(msh.faces.begin(), msh.faces.end(), f);
    if (itor == msh.faces.end())
        throw std::logic_error("Face not found, this is likely a bug.");
    ret[0] = *itor;

    f.ptids[0] = cl.ptids[1];
    f.ptids[1] = cl.ptids[2];
    itor = std::lower_bound(msh.faces.begin(), msh.faces.end(), f);
    if (itor == msh.faces.end())
        throw std::logic_error("Face not found, this is likely a bug.");
    ret[1] = *itor;

    f.ptids[0] = cl.ptids[3];
    f.ptids[1] = cl.ptids[2];
    itor = std::lower_bound(msh.faces.begin(), msh.faces.end(), f);
    if (itor == msh.faces.end())
        throw std::logic_error("Face not found, this is likely a bug.");
    ret[2] = *itor;

    f.ptids[0] = cl.ptids[0];
    f.ptids[1] = cl.ptids[3];
    itor = std::lower_bound(msh.faces.begin(), msh.faces.end(), f);
    if (itor == msh.faces.end())
        throw std::logic_error("Face not found, this is likely a bug.");
    ret[3] = *itor;

    return ret;
}

template<typename Mesh>
typename Mesh::point_type
barycenter(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    auto pts = points(msh, cl);
    return (pts[0] + pts[1] + pts[2] + pts[3])/4.0;
}

template<typename Mesh>
typename Mesh::point_type
barycenter(const Mesh& msh, const typename Mesh::face_type& fc)
{
    auto pts = points(msh, fc);
    return (pts[0] + pts[1])/2.0;
}

template<typename Mesh>
typename Mesh::coordinate_type
diameter(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    auto p0 = msh.points.at( cl.ptids[0] );
    auto p1 = msh.points.at( cl.ptids[2] );

    return (p1-p0).to_vector().norm();
}

template<typename Mesh>
typename Mesh::coordinate_type
diameter(const Mesh& msh, const typename Mesh::face_type& fc)
{
    auto p0 = msh.points.at( fc.ptids[0] );
    auto p1 = msh.points.at( fc.ptids[1] );

    return (p1-p0).to_vector().norm();
}

template<typename Mesh>
typename Mesh::coordinate_type
measure(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    auto p0 = msh.points.at( cl.ptids[0] );
    auto p2 = msh.points.at( cl.ptids[2] );

    return (p2.x() - p0.x()) * (p2.y() - p0.y());
}

template<typename Mesh>
typename Mesh::coordinate_type
measure(const Mesh& msh, const typename Mesh::face_type& fc)
{
    auto p0 = msh.points.at( fc.ptids[0] );
    auto p1 = msh.points.at( fc.ptids[1] );

    return (p1-p0).to_vector().norm();
}

template<typename Mesh>
std::array< Matrix<typename Mesh::coordinate_type,2,1>, 4 >
normals(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    using T = typename Mesh::coordinate_type;
    std::array< Matrix<T,2,1>, 4 >  ret;

    ret[0](0) = 0;
    ret[0](1) = -1;

    ret[1](0) = 1;
    ret[1](1) = 0;

    ret[2](0) = 0;
    ret[2](1) = 1;

    ret[3](0) = -1;
    ret[3](1) = 0;

    return ret;
}

template<typename Mesh>
std::vector< typename Mesh::point_type >
make_test_points(const Mesh& msh, const typename Mesh::cell_type& cl)
{
    const size_t N = 10;

    std::vector< typename Mesh::point_type > ret;

    auto cell_pts = points(msh, cl);

    auto min_x = cell_pts[0].x();
    auto min_y = cell_pts[0].y();

    auto hx = (cell_pts[2].x() - cell_pts[0].x())/N;
    auto hy = (cell_pts[2].y() - cell_pts[0].y())/N;

    for (size_t j = 0; j < N+1; j++)
    {
        for(size_t i = 0; i < N+1; i++)
        {
            auto px = min_x + i*hx;
            auto py = min_y + j*hy;
            typename Mesh::point_type pt(px, py);
            ret.push_back(pt);
        }
    }

    return ret;
}

template<typename Mesh>
std::vector< typename Mesh::point_type >
make_test_points(const Mesh& msh, const typename Mesh::face_type& fc)
{
    const size_t N = 10;

    std::vector< typename Mesh::point_type > ret;

    auto face_pts = points(msh, fc);

    auto min = face_pts[0];

    auto h = (face_pts[1] - face_pts[0])/N;

    for (size_t i = 0; i < N+1; i++)
    {
        auto p = min + i*h;
        ret.push_back(p);
    }

    return ret;
}
