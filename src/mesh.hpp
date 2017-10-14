/*
 *       /\        Matteo Cicuttin (C) 2017
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    cutHHO prototype code.
 *  /_\/_\/_\/_\
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it.
 *
 */

#pragma once

#include <vector>
#include <fstream>
#include <algorithm>

#include "src/point.hpp"


struct node
{
    size_t                  point_id;
};

struct edge
{
    std::array<size_t, 2>   point_ids;
    bool                    is_boundary;
    size_t                  boundary_id;

    bool operator<(const edge& other) const
    {
        return point_ids < other.point_ids;
    }

    bool operator==(const edge& other) const
    {
        return point_ids == other.point_ids;
    }
};

struct polygon
{
    std::vector<size_t>     point_ids;
    size_t                  domain_id;
    bool                    cut_by_interface;

    polygon()
    {
        point_ids.clear();
        domain_id = 0;
        cut_by_interface = false;
    }

    bool operator<(const polygon& other) const
    {
        return point_ids < other.point_ids;
    }

    bool operator==(const polygon& other) const
    {
        return point_ids == other.point_ids;
    }
};


template<typename T>
struct mesh
{
    typedef point<T,2>          point_type;
    typedef node                node_type;
    typedef edge                face_type;
    typedef polygon             cell_type;

    std::vector<point<T,2>>     points;
    std::vector<node_type>      nodes;
    std::vector<face_type>      faces;
    std::vector<cell_type>      cells;
};

template<typename T>
typename mesh<T>::point_type
points(const mesh<T>& msh, const typename mesh<T>::node_type& node)
{
    return msh.points.at( node.point_id );
}

template<typename T>
std::array< typename mesh<T>::point_type, 2 >
points(const mesh<T>& msh, const typename mesh<T>::face_type& face)
{
    std::array< typename mesh<T>::point_type, 2 >    ret;

    auto id_to_point = [&](size_t ptid) -> typename mesh<T>::point_type {
        return msh.points.at( ptid );
    };

    std::transform(face.point_ids.begin(), face.point_ids.end(), ret.begin());
    return ret;
}

template<typename T>
std::vector<typename mesh<T>::point_type>
points(const mesh<T>& msh, const typename mesh<T>::cell_type& cell)
{
    std::vector<typename mesh<T>::point_type>    ret;

    ret.resize( cell.point_ids.size() );

    auto id_to_point = [&](size_t ptid) -> typename mesh<T>::point_type {
        return msh.points.at( ptid );
    };

    std::transform(cell.point_ids.begin(), cell.point_ids.end(), ret.begin(), id_to_point);
    return ret;
}

template<typename T, typename Function>
void
detect_cut_cells(mesh<T>& msh, const Function& level_set_function)
{
    for (auto& cl : msh.cells)
    {
        auto pts = points(msh, cl);
        std::vector<int> signs;
        signs.resize(pts.size());

        auto compute_sign = [&](typename mesh<T>::point_type pt) -> int {
            return level_set_function(pt) < 0 ? -1 : +1;
        };

        std::transform(pts.begin(), pts.end(), signs.begin(), compute_sign);

        cl.cut_by_interface = true;

        if ( std::count(signs.begin(), signs.end(), 1) == signs.size() )
            cl.cut_by_interface = false;

        if ( std::count(signs.begin(), signs.end(), -1) == signs.size() )
            cl.cut_by_interface = false;
    }
}


template<typename T>
bool load_netgen_2d(const std::string& filename, mesh<T>& msh)
{
    std::ifstream ifs(filename);
    if (!ifs.is_open())
    {
        std::cout << "Problem opening input mesh" << std::endl;
        return false;
    }

    size_t count;

    ifs >> count;
    msh.points.reserve(count);

    for (size_t i = 0; i < count; i++)
    {
        T x, y;
        ifs >> x >> y;
        msh.points.push_back( typename mesh<T>::point_type(x,y) );
    }

    std::cout << "Points: " << msh.points.size() << std::endl;

    ifs >> count;
    msh.cells.reserve(count);

    for (size_t i = 0; i < count; i++)
    {
        size_t dom, p0, p1, p2;
        ifs >> dom >> p0 >> p1 >> p2;

        if (dom == 0 || p0 == 0 || p1 == 0 || p2 == 0 )
        {
            std::cout << "Indices in netgen file should be 1-based, found a 0.";
            std::cout << std::endl;
            return false;
        }

        polygon p;
        p.domain_id = dom-1;
        p.point_ids.push_back(p0-1);
        p.point_ids.push_back(p1-1);
        p.point_ids.push_back(p2-1);
        msh.cells.push_back(p);
    }

    std::cout << "Cells: " << msh.cells.size() << std::endl;
    std::sort(msh.cells.begin(), msh.cells.end());

    msh.faces.reserve( 3 * msh.cells.size() );

    for (auto& cl : msh.cells)
    {
        assert(cl.point_ids.size() == 3);

        for (size_t i = 0; i < 3; i++)
        {
            auto p0 = cl.point_ids[i];
            auto p1 = cl.point_ids[i%3];

            if (p1 < p0)
                std::swap(p0, p1);

            typename mesh<T>::face_type face;
            face.point_ids[0] = p0;
            face.point_ids[1] = p1;

            msh.faces.push_back(face);
        }
    }

    std::sort(msh.faces.begin(), msh.faces.end());
    msh.faces.erase(msh.faces.begin(), std::unique(msh.faces.begin(), msh.faces.end()) );

    std::cout << "Faces: " << msh.faces.size() << std::endl;

    for (size_t i = 0; i < msh.points.size(); i++)
    {
        typename mesh<T>::node_type node;
        node.point_id = i;
        msh.nodes.push_back(node);
    }

    std::cout << "Nodes: " << msh.nodes.size() << std::endl;

    ifs.close();
    return true;
}
