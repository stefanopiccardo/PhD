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
#include "cuthho_mesh.hpp"

template<typename T>
bool
is_cut(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::cell_type& cl)
{
    assert(cl.user_data.location != element_location::UNDEF);
    return cl.user_data.location == element_location::ON_INTERFACE;
}

template<typename T>
bool
is_cut(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::face_type& fc)
{
    assert(fc.user_data.location != element_location::UNDEF);
    return fc.user_data.location == element_location::ON_INTERFACE;
}

template<typename T>
element_location
location(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::cell_type& cl)
{
    assert(cl.user_data.location != element_location::UNDEF);
    return cl.user_data.location;
}

template<typename T>
element_location
location(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::face_type& fc)
{
    assert(fc.user_data.location != element_location::UNDEF);
    return fc.user_data.location;
}

template<typename T>
element_location
location(const cuthho_mesh<T>& msh, const typename cuthho_mesh<T>::node_type& n)
{
    assert(n.user_data.location != element_location::UNDEF);
    return n.user_data.location;
}


template<typename T, typename Function>
point<T, 2>
find_zero_crossing(const point<T,2>& p0, const point<T,2>& p1, const Function& level_set_function,
                   const T& threshold)
{
    /* !!! We assume that the level set function *has* a zero crossing
     * between p0 and p1 !!! */
    auto pa = p0;
    auto pb = p1;
    auto pm = (pa+pb)/2.0;
    auto pm_prev = pm;

    T x_diff_sq, y_diff_sq;

    /* A threshold of 1/10000 the diameter of the element is considered
     * acceptable. Since with 24 iterations we reduce the error by 16384
     * and the worst case is that the two points are at the opposite sides
     * of the element, we put 30 as limit. */
    size_t max_iter = 30;

    do {
        auto la = level_set_function(pa);
        auto lb = level_set_function(pb);
        auto lm = level_set_function(pm);

        if ( (lb >= 0 && lm >= 0) || (lb < 0 && lm < 0) )
        {   /* intersection is between pa and pm */
            pm_prev = pm;
            pb = pm;
            pm = (pa+pb)/2.0;
        }
        else
        {   /* intersection is between pm and pb */
            pm_prev = pm;
            pa = pm;
            pm = (pa+pb)/2.0;
        }

        x_diff_sq = (pm_prev.x() - pm.x()) * (pm_prev.x() - pm.x());
        y_diff_sq = (pm_prev.y() - pm.y()) * (pm_prev.y() - pm.y());

    } while ( (sqrt(x_diff_sq + y_diff_sq) > threshold) && max_iter-- );

    return pm;

    /* Affine zero crossing was like that: */
    //auto t = l0/(l0-l1);
    //auto ip = (pts[1] - pts[0]) * t + pts[0];
}

template<typename T, typename Function>
void
detect_node_position(cuthho_mesh<T>& msh, const Function& level_set_function)
{
    for (auto& n : msh.nodes)
    {
        auto pt = points(msh, n);
        if ( level_set_function(pt) < 0 )
            n.user_data.location = element_location::IN_NEGATIVE_SIDE;
        else
            n.user_data.location = element_location::IN_POSITIVE_SIDE;
    }
}

template<typename T, typename Function>
void
detect_cut_faces(cuthho_mesh<T>& msh, const Function& level_set_function)
{
    for (auto& fc : msh.faces)
    {
        auto pts = points(msh, fc);
        auto l0 = level_set_function(pts[0]);
        auto l1 = level_set_function(pts[1]);
        if (l0 >= 0 && l1 >= 0)
        {
            fc.user_data.location = element_location::IN_POSITIVE_SIDE;
            continue;
        }

        if (l0 < 0 && l1 < 0)
        {
            fc.user_data.location = element_location::IN_NEGATIVE_SIDE;
            continue;
        }

        auto threshold = diameter(msh, fc) / 10000.0;
        auto pm = find_zero_crossing(pts[0], pts[1], level_set_function, threshold);

        /* If node 0 is in the negative region, mark it as node inside, otherwise mark node 1 */
        fc.user_data.node_inside = ( l0 < 0 ) ? 0 : 1;
        fc.user_data.location = element_location::ON_INTERFACE;
        fc.user_data.intersection_point = pm;
    }
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
        auto fcs = faces(msh, cl);
        std::array< std::pair<size_t, point_type>, 2 >  cut_faces;

        size_t k = 0;
        for (size_t i = 0; i < fcs.size(); i++)
        {
            if ( is_cut(msh, fcs[i]) )
                cut_faces.at(k++) = std::make_pair(i, fcs[i].user_data.intersection_point);
        }

        /* If a face is cut, the cells that own the face are cut. Is this
         * unconditionally true? It should...fortunately this isn't avionics
         * software */

        if (k == 0)
        {
            auto is_positive = [&](const point_type& pt) -> bool {
                return level_set_function(pt) > 0;
            };

            auto pts = points(msh, cl);
            if ( std::all_of(pts.begin(), pts.end(), is_positive) )
                cl.user_data.location = element_location::IN_POSITIVE_SIDE;
            else
                cl.user_data.location = element_location::IN_NEGATIVE_SIDE;
        }

        if (k == 2)
        {
            cl.user_data.location = element_location::ON_INTERFACE;
            auto p0 = cut_faces[0].second;
            auto p1 = cut_faces[1].second;
            auto pt = p1 - p0;
            auto pn = p0 + point<T,2>(-pt.y(), pt.x());

            if ( level_set_function(pn) >= 0 )
            {
                cl.user_data.p0 = p1;
                cl.user_data.p1 = p0;
            }
            else
            {
                cl.user_data.p0 = p0;
                cl.user_data.p1 = p1;
            }

            cl.user_data.interface.push_back(cl.user_data.p0);
            cl.user_data.interface.push_back(cl.user_data.p1);
        }

        if ( k != 0 && k != 2 )
            throw std::logic_error("invalid number of cuts in cell");

        cell_i++;
    }
}

template<typename T, typename Function>
void
refine_interface(cuthho_mesh<T>& msh, typename cuthho_mesh<T>::cell_type& cl,
                 const Function& level_set_function, size_t min, size_t max)
{
    if ( (max-min) < 2 )
        return;

    typedef typename cuthho_mesh<T>::point_type     point_type;

    size_t mid = (max+min)/2;
    auto p0 = cl.user_data.interface.at(min);
    auto p1 = cl.user_data.interface.at(max);
    auto pm = (p0+p1)/2.0;
    auto pt = p1 - p0;
    auto pn = point_type(-pt.y(), pt.x());
    auto ps1 = pm + pn;
    auto ps2 = pm - pn;

    auto lm = level_set_function(pm);
    auto ls1 = level_set_function(ps1);
    auto ls2 = level_set_function(ps2);

    point_type ip;

    if ( !((lm >= 0 && ls1 >= 0) || (lm < 0 && ls1 < 0)) )
    {
        auto threshold = diameter(msh, cl) / 10000.0;
        ip = find_zero_crossing(pm, ps1, level_set_function, threshold);
    }
    else if ( !((lm >= 0 && ls2 >= 0) || (lm < 0 && ls2 < 0)) )
    {
        auto threshold = diameter(msh, cl) / 10000.0;
        ip = find_zero_crossing(pm, ps2, level_set_function, threshold);
    }
    else
        throw std::logic_error("interface not found in search range");

    cl.user_data.interface.at(mid) = ip;

    refine_interface(msh, cl, level_set_function, min, mid);
    refine_interface(msh, cl, level_set_function, mid, max);
}

template<typename T, typename Function>
void
refine_interface(cuthho_mesh<T>& msh, const Function& level_set_function, size_t levels)
{
    if (levels == 0)
        return;

    size_t interface_points = iexp_pow(2, levels);

    for (auto& cl : msh.cells)
    {
        if ( !is_cut(msh, cl) )
            continue;

        cl.user_data.interface.resize(interface_points+1);
        cl.user_data.interface.at(0)                = cl.user_data.p0;
        cl.user_data.interface.at(interface_points) = cl.user_data.p1;

        refine_interface(msh, cl, level_set_function, 0, interface_points);
    }
}

template<typename T>
std::vector< typename cuthho_mesh<T>::point_type >
collect_triangulation_points(const cuthho_mesh<T>& msh,
                             const typename cuthho_mesh<T>::cell_type& cl,
                             const element_location& where)
{
    typedef typename cuthho_mesh<T>::point_type     point_type;
    typedef typename cuthho_mesh<T>::node_type      node_type;

    assert( is_cut(msh, cl) );
    auto ns = nodes(msh, cl);

    std::vector< point_type > ret;

    auto node2pt = [&](const node_type& n) -> auto {
        return msh.points.at(n.ptid);
    };

    auto insert_interface = [&](void) -> void {
        if (where == element_location::IN_NEGATIVE_SIDE)
            ret.insert(ret.end(), cl.user_data.interface.begin(), cl.user_data.interface.end());
        else if (where == element_location::IN_POSITIVE_SIDE)
            ret.insert(ret.end(), cl.user_data.interface.rbegin(), cl.user_data.interface.rend());
        else
            throw std::logic_error("If you've got here there is some issue...");
    };

    bool case1 = location(msh, ns.front()) == where && location(msh, ns.back()) != where;
    bool case2 = location(msh, ns.front()) != where && location(msh, ns.back()) == where;
    bool case3 = location(msh, ns.front()) != where && location(msh, ns.back()) != where;
    //bool case4 = location(msh, ns.front()) == where && location(msh, ns.back()) == where;

    if ( case1 || case2 || case3 )
    {
        for (size_t i = 0; i < ns.size(); i++)
            if ( location(msh, ns[i]) == where )
                ret.push_back( node2pt(ns[i]) );

        insert_interface();
    }
    else // case4, the only possible remaining
    {
        size_t i = 0;
        while ( i < ns.size() && location(msh, ns.at(i)) == where )
            ret.push_back( node2pt(ns.at(i++)) );
        insert_interface();
        while ( i < ns.size() && location(msh, ns.at(i)) != where )
            i++;
        while ( i < ns.size() && location(msh, ns.at(i)) == where )
            ret.push_back( node2pt(ns.at(i++)) );
    }

    return ret;
}

template<typename T>
auto
barycenter(const std::vector< point<T, 2> >& pts)
{
    point<T, 2> ret;

    T den = 0.0;

    for (size_t i = 2; i < pts.size(); i++)
    {
        auto pprev  = pts[i-1]-pts[0];
        auto pcur   = pts[i]-pts[0];
        auto d      = det(pprev, pcur) / 2.0;
        ret         = ret + (pprev + pcur) * d;
        den         += d;
    }

    return pts[0] + ret/(den*3);
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
        else if ( is_cut(msh, fc) )
            ofs << "line([" << pts[0].x() << ", " << pts[1].x() << "], [" << pts[0].y() << ", " << pts[1].y() << "], 'Color', 'g');" << std::endl;
        else
            ofs << "line([" << pts[0].x() << ", " << pts[1].x() << "], [" << pts[0].y() << ", " << pts[1].y() << "], 'Color', 'k');" << std::endl;

        auto bar = barycenter(msh, fc);
        ofs << "text(" << bar.x() << ", " << bar.y() << ", '" << i << "');" << std::endl;

        i++;
    }

    for (auto& cl : msh.cells)
    {
        if ( is_cut(msh, cl) )
        {
            auto p0 = cl.user_data.p0;
            auto p1 = cl.user_data.p1;
            auto q = p1 - p0;
            ofs << "quiver(" << p0.x() << ", " << p0.y() << ", " << q.x() << ", " << q.y() << ", 0)" << std::endl;

            for (size_t i = 1; i < cl.user_data.interface.size(); i++)
            {
                auto s = cl.user_data.interface.at(i-1);
                auto l = cl.user_data.interface.at(i) - s;
                ofs << "quiver(" << s.x() << ", " << s.y() << ", " << l.x() << ", " << l.y() << ", 0)" << std::endl;
            }

            for (auto& ip : cl.user_data.interface)
                ofs << "plot(" << ip.x() << ", " << ip.y() << ", '*k');" << std::endl;

            auto tpn = collect_triangulation_points(msh, cl, element_location::IN_NEGATIVE_SIDE);
            auto bn = barycenter(tpn);
            ofs << "plot(" << bn.x() << ", " << bn.y() << ", 'dr');" << std::endl;

            auto tpp = collect_triangulation_points(msh, cl, element_location::IN_POSITIVE_SIDE);
            auto bp = barycenter(tpp);
            ofs << "plot(" << bp.x() << ", " << bp.y() << ", 'db');" << std::endl;

        }
    }

    ofs << "t = linspace(0,2*pi,1000);plot(sqrt(0.15)*cos(t)+0.5, sqrt(0.15)*sin(t)+0.5)" << std::endl;

    ofs.close();
}
