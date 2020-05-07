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

#include <unsupported/Eigen/MatrixFunctions> // ADD BY STEFANO

using namespace Eigen;

#include "core/core"
#include "core/solvers"
#include "dataio/silo_io.hpp"

#include "methods/hho"
#include "methods/cuthho"



///////////////////////   STEFANO FUNCTIONS  ///////////////////////////
/// NEW VERSIONE OF PREVIOUS FUNCTIONS REWRITTEN FOR DISCRETE LEVEL SET FUNCTIONS

/// Each mesh has memory of the original nodes: useful to link agglomerated - original mesh
template<typename Mesh>
void offset_definition( Mesh& msh)
{
    //size_t counter = 0;
    for (auto& cl:msh.cells) {
        auto offset_orig = offset(msh,cl);
        cl.user_data.offset_subcells.push_back(offset_orig);
       // counter++;
       // std::cout<<"Initialisation cell num "<< cl.user_data.offset_subcells[offset_orig]<<std::endl;
    }
}


/// Useful to plot level set pre FE transport problem
/// in cuthho_export.hpp
template<typename Mesh, typename Function>
void
output_mesh_info2_pre_FEM(const Mesh& msh, const Function& level_set_function)
{
    using RealType = typename Mesh::coordinate_type;

    /************** OPEN SILO DATABASE **************/
    silo_database silo;
    silo.create("cuthho_meshinfo.silo");
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

    /************** MAKE A SILO VARIABLE FOR CELL HIGHLIGHT **************/
    std::vector<RealType> highlight_markers;
    for (auto& cl : msh.cells)
    {
        if ( cl.user_data.highlight )
            highlight_markers.push_back(1.0);
        else
            highlight_markers.push_back(0.0);

    }
    silo.add_variable("mesh", "highlighted_cells", highlight_markers.data(), highlight_markers.size(), zonal_variable_t);

    /************** MAKE A SILO VARIABLE FOR LEVEL SET FUNCTION **************/
    std::vector<RealType> level_set_vals;
    // for (auto& pt : msh.points)
    //    level_set_vals.push_back( level_set_function(pt) );
    for (auto& n : msh.nodes)
        level_set_vals.push_back( level_set_function(n) );
    
    silo.add_variable("mesh", "level_set", level_set_vals.data(), level_set_vals.size(), nodal_variable_t);

    /************** MAKE A SILO VARIABLE FOR NODE POSITIONING **************/
    std::vector<RealType> node_pos;
    for (auto& n : msh.nodes)
        node_pos.push_back( location(msh, n) == element_location::IN_POSITIVE_SIDE ? +1.0 : -1.0 );
    silo.add_variable("mesh", "node_pos", node_pos.data(), node_pos.size(), nodal_variable_t);

    std::vector<RealType> cell_set;
    for (auto& cl : msh.cells)
    {
        RealType r;

        switch ( cl.user_data.agglo_set )
        {
            case cell_agglo_set::UNDEF:
                r = 0.0;
                break;

            case cell_agglo_set::T_OK:
                r = 1.0;
                break;

            case cell_agglo_set::T_KO_NEG:
                r = 2.0;
                break;

            case cell_agglo_set::T_KO_POS:
                r = 3.0;
                break;

        }

        cell_set.push_back( r );
    }
    silo.add_variable("mesh", "agglo_set", cell_set.data(), cell_set.size(), zonal_variable_t);

    silo.close();

    /*************  MAKE AN OUTPUT FOR THE INTERSECTION POINTS *************/
    std::vector<RealType> int_pts_x;
    std::vector<RealType> int_pts_y;

    for (auto& fc : msh.faces)
    {
        if( fc.user_data.location != element_location::ON_INTERFACE ) continue;

        RealType x = fc.user_data.intersection_point.x();
        RealType y = fc.user_data.intersection_point.y();

        int_pts_x.push_back(x);
        int_pts_y.push_back(y);
    }

    std::ofstream points_file_pre_FEM("int_points_pre_FEM.3D", std::ios::out | std::ios::trunc);

    if(points_file_pre_FEM)
    {
        // instructions
        points_file_pre_FEM << "X   Y   Z   val" << std::endl;

        for( size_t i = 0; i<int_pts_x.size(); i++)
        {
            points_file_pre_FEM << int_pts_x[i] << "   " <<  int_pts_y[i]
                        << "   0.0     0.0" << std::endl;
        }

        points_file_pre_FEM.close();
    }

    else
        std::cerr << "points_file_pre_FEM has not been opened" << std::endl;


    /*************  MAKE AN OUTPUT FOR THE INTERFACE *************/
    std::vector<RealType> int_x;
    std::vector<RealType> int_y;

    for (auto& cl : msh.cells)
    {
        if( cl.user_data.location != element_location::ON_INTERFACE ) continue;

        for(size_t i = 0; i < cl.user_data.interface.size(); i++)
        {
            RealType x = cl.user_data.interface.at(i).x();
            RealType y = cl.user_data.interface.at(i).y();

            int_x.push_back(x);
            int_y.push_back(y);
        }
    }
    std::ofstream interface_file_preFEM("interface_pre_FEM.3D", std::ios::out | std::ios::trunc);

    if(interface_file_preFEM)
    {
        // instructions
        interface_file_preFEM << "X   Y   Z   val" << std::endl;

        for( size_t i = 0; i<int_x.size(); i++)
        {
            interface_file_preFEM << int_x[i] << "   " <<  int_y[i]
                        << "   0.0     0.0" << std::endl;
        }

        interface_file_preFEM.close();
    }

    else
        std::cerr << "interface_file_preFEM has not been opened" << std::endl;
}


/// Useful to plot level set post FE transport problem
/// in cuthho_export.hpp
template<typename Mesh, typename Function>
void
output_mesh_info2(const Mesh& msh, const Function& level_set_function)
{
    using RealType = typename Mesh::coordinate_type;

    /************** OPEN SILO DATABASE **************/
    silo_database silo;
    silo.create("cuthho_meshinfo.silo");
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

    /************** MAKE A SILO VARIABLE FOR CELL HIGHLIGHT **************/
    std::vector<RealType> highlight_markers;
    for (auto& cl : msh.cells)
    {
        if ( cl.user_data.highlight )
            highlight_markers.push_back(1.0);
        else
            highlight_markers.push_back(0.0);

    }
    silo.add_variable("mesh", "highlighted_cells", highlight_markers.data(), highlight_markers.size(), zonal_variable_t);

    /************** MAKE A SILO VARIABLE FOR LEVEL SET FUNCTION **************/
    std::vector<RealType> level_set_vals;
    // for (auto& pt : msh.points)
    //    level_set_vals.push_back( level_set_function(pt) );
    for (auto& n : msh.nodes)
        level_set_vals.push_back( level_set_function(n) );
    
    silo.add_variable("mesh", "level_set", level_set_vals.data(), level_set_vals.size(), nodal_variable_t);

    /************** MAKE A SILO VARIABLE FOR NODE POSITIONING **************/
    std::vector<RealType> node_pos;
    for (auto& n : msh.nodes)
        node_pos.push_back( location(msh, n) == element_location::IN_POSITIVE_SIDE ? +1.0 : -1.0 );
    silo.add_variable("mesh", "node_pos", node_pos.data(), node_pos.size(), nodal_variable_t);

    std::vector<RealType> cell_set;
    for (auto& cl : msh.cells)
    {
        RealType r;

        switch ( cl.user_data.agglo_set )
        {
            case cell_agglo_set::UNDEF:
                r = 0.0;
                break;

            case cell_agglo_set::T_OK:
                r = 1.0;
                break;

            case cell_agglo_set::T_KO_NEG:
                r = 2.0;
                break;

            case cell_agglo_set::T_KO_POS:
                r = 3.0;
                break;

        }

        cell_set.push_back( r );
    }
    silo.add_variable("mesh", "agglo_set", cell_set.data(), cell_set.size(), zonal_variable_t);

    silo.close();

    /*************  MAKE AN OUTPUT FOR THE INTERSECTION POINTS *************/
    std::vector<RealType> int_pts_x;
    std::vector<RealType> int_pts_y;

    for (auto& fc : msh.faces)
    {
        if( fc.user_data.location != element_location::ON_INTERFACE ) continue;

        RealType x = fc.user_data.intersection_point.x();
        RealType y = fc.user_data.intersection_point.y();

        int_pts_x.push_back(x);
        int_pts_y.push_back(y);
    }

    std::ofstream points_file("int_points.3D", std::ios::out | std::ios::trunc);

    if(points_file)
    {
        // instructions
        points_file << "X   Y   Z   val" << std::endl;

        for( size_t i = 0; i<int_pts_x.size(); i++)
        {
            points_file << int_pts_x[i] << "   " <<  int_pts_y[i]
                        << "   0.0     0.0" << std::endl;
        }

        points_file.close();
    }

    else
        std::cerr << "Points_file has not been opened" << std::endl;


    /*************  MAKE AN OUTPUT FOR THE INTERFACE *************/
    std::vector<RealType> int_x;
    std::vector<RealType> int_y;

    for (auto& cl : msh.cells)
    {
        if( cl.user_data.location != element_location::ON_INTERFACE ) continue;

        for(size_t i = 0; i < cl.user_data.interface.size(); i++)
        {
            RealType x = cl.user_data.interface.at(i).x();
            RealType y = cl.user_data.interface.at(i).y();

            int_x.push_back(x);
            int_y.push_back(y);
        }
    }
    std::ofstream interface_file("interface.3D", std::ios::out | std::ios::trunc);

    if(interface_file)
    {
        // instructions
        interface_file << "X   Y   Z   val" << std::endl;

        for( size_t i = 0; i<int_x.size(); i++)
        {
            interface_file << int_x[i] << "   " <<  int_y[i]
                        << "   0.0     0.0" << std::endl;
        }

        interface_file.close();
    }

    else
        std::cerr << "Interface_file has not been opened" << std::endl;
}

/// New version of find_zero_crossing for discret level set functions
/// in cuthho_geom.hpp
template<typename T, typename Function, typename Mesh>
point<T, 2>
find_zero_crossing_on_face(const point<T,2>& p0, const point<T,2>& p1, const Function& level_set_function,
                   const T& threshold, const Mesh & msh, const typename Mesh::face_type& fc)
{
    /* !!! We assume that the level set function *has* a zero crossing
     * between p0 and p1 !!! */
    
    // I SHOULD CHECK THAT pm IS ALWAYS ON THE FACE p0-p1 ???????
    auto pa = p0;
    auto pb = p1;
    auto pm = (pa+pb)/2.0;
    auto pm_prev = pm;

    T x_diff_sq, y_diff_sq;

    /* A threshold of 1/10000 the diameter of the element is considered
     * acceptable. Since with 24 iterations we reduce the error by 16384
     * and the worst case is that the two points are at the opposite sides
     * of the element, we put 30 as limit. */
    size_t max_iter = 50;

    do {
        auto la = level_set_function(pa,msh,fc);
        auto lb = level_set_function(pb,msh,fc);
        auto lm = level_set_function(pm,msh,fc);

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


/// New version of find_zero_crossing for discret level set functions
/// in cuthho_geom.hpp
template<typename T, typename Function, typename Mesh >
point<T, 2>
find_zero_crossing_in_cell(const point<T,2>& p0, const point<T,2>& p1, const Function& level_set_function,
                   const T& threshold, const Mesh & msh, const typename Mesh::cell_type& cl)
{
    /* !!! We assume that the level set function *has* a zero crossing
     * between p0 and p1 !!! */
    
    // I SHOULD CHECK THAT pm IS ALWAYS IN THE CELL ???????
    auto pa = p0;
    auto pb = p1;
    auto pm = (pa+pb)/2.0;
    auto pm_prev = pm;

    T x_diff_sq, y_diff_sq;

    /* A threshold of 1/10000 the diameter of the element is considered
     * acceptable. Since with 24 iterations we reduce the error by 16384
     * and the worst case is that the two points are at the opposite sides
     * of the element, we put 30 as limit. */
    size_t max_iter = 50; // ERA 50, METTO 100

    do {
        auto la = level_set_function(pa,msh,cl);
        auto lb = level_set_function(pb,msh,cl);
        auto lm = level_set_function(pm,msh,cl);

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




/// New version of detect_node_position for discret level functions
/// in cuthho_geom.hpp
template<typename T, size_t ET, typename Function>
void
detect_node_position2(cuthho_mesh<T, ET>& msh, const Function& level_set_function)
{
    for (auto& n : msh.nodes)
    {
        //auto pt = points(msh, n); //deleted by Stefano
        //if ( level_set_function(pt) < 0 ) //deleted by Stefano
        if ( level_set_function(n) < 0 ) // add by Stefano
            n.user_data.location = element_location::IN_NEGATIVE_SIDE;
        else
            n.user_data.location = element_location::IN_POSITIVE_SIDE;
    }
}

/// New version of detect_cut_faces for discret level functions
/// in cuthho_geom.hpp
template<typename T, size_t ET, typename Function>
void
detect_cut_faces2(cuthho_mesh<T, ET>& msh, const Function& level_set_function)
{
    for (auto& fc : msh.faces)
    {
        auto pts = points(msh, fc);
        //auto l0 = level_set_function(pts[0]);      //deleted by Stefano
        //auto l1 = level_set_function(pts[1]);       //deleted by Stefano
        
        auto l0 = level_set_function(pts[0],msh,fc);      // add by Stefano
        auto l1 = level_set_function(pts[1],msh,fc);       // add by Stefano
        
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

        auto threshold = diameter(msh, fc) / 1e20;
        //auto pm = find_zero_crossing(pts[0], pts[1], level_set_function, threshold);
        auto pm = find_zero_crossing_on_face(pts[0], pts[1], level_set_function, threshold,msh,fc);

        /* If node 0 is in the negative region, mark it as node inside, otherwise mark node 1 */
        fc.user_data.node_inside = ( l0 < 0 ) ? 0 : 1;
        fc.user_data.location = element_location::ON_INTERFACE;
        fc.user_data.intersection_point = pm;
    }
}


/// New version of detect_cut_cells for discret level functions
/// in cuthho_geom.hpp
template<typename T, size_t ET, typename Function>
void
detect_cut_cells2(cuthho_mesh<T, ET>& msh, const Function& level_set_function)
{
    std::cout<<"I AM IN DETECT CUT CELL2!!!!"<<std::endl;
    typedef typename cuthho_mesh<T, ET>::face_type  face_type;
    typedef typename cuthho_mesh<T, ET>::point_type point_type;
    typedef typename cuthho_mesh<T, ET>::cell_type cell_type;

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
     
            
            
            
            /*
            auto pts = points(msh, cl);
            auto pt = pts.begin();
            size_t counter = 0;
            while( ( pt!= pts.end() ) && ( level_set_function(*pt,msh,cl) > 0 ) )
            {
                counter++;
                pt++;
                
            }
             
            if ( counter == pts.size() )
                cl.user_data.location = element_location::IN_POSITIVE_SIDE;
            else
                cl.user_data.location = element_location::IN_NEGATIVE_SIDE;
            */
             
        }

        if (k == 2)
        {
            cl.user_data.location = element_location::ON_INTERFACE;
            auto p0 = cut_faces[0].second;
            auto p1 = cut_faces[1].second;
            auto pt = p1 - p0;
            auto pn = p0 + point<T,2>(-pt.y(), pt.x());

            if ( level_set_function(pn,msh,cl) >= 0 )
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


/// New version of refine_interface for discret level functions

template<typename T, size_t ET, typename Function>
void
refine_interface2(cuthho_mesh<T, ET>& msh, typename cuthho_mesh<T, ET>::cell_type& cl,
                 const Function& level_set_function, size_t min, size_t max)
{
    if ( (max-min) < 2 )
        return;

    typedef typename cuthho_mesh<T, ET>::point_type     point_type;

    size_t mid = (max+min)/2;
    auto p0 = cl.user_data.interface.at(min);
    auto p1 = cl.user_data.interface.at(max);
    auto pm = (p0+p1)/2.0;
    auto pt = p1 - p0;
    auto pn = point_type(-pt.y(), pt.x());
    auto ps1 = pm + pn;
    auto ps2 = pm - pn;

    auto lm = level_set_function(pm,msh,cl);
    auto ls1 = level_set_function(ps1,msh,cl);
    auto ls2 = level_set_function(ps2,msh,cl);

    point_type ip;
   // std::cout<<"the node of interface are "<<p0<<" and "<<p1<<". I search pm= "<<pm<<" in which phi = "<<lm<<" and ps1 e ps2 "<<ps1<<" and "<<ps2<<"equal to "<<ls1<<" , "<<ls2<<std::endl;
    if ( !((lm >= 0 && ls1 >= 0) || (lm < 0 && ls1 < 0)) )
    {
        auto threshold = diameter(msh, cl) / 1e20;
        ip = find_zero_crossing_in_cell(pm, ps1, level_set_function, threshold,msh,cl);
    }
    else if ( !((lm >= 0 && ls2 >= 0) || (lm < 0 && ls2 < 0)) )
    {
        auto threshold = diameter(msh, cl) / 1e20;
        ip = find_zero_crossing_in_cell(pm, ps2, level_set_function, threshold,msh,cl);
    }
    else
        throw std::logic_error("interface not found in search range");

    cl.user_data.interface.at(mid) = ip;

    refine_interface2(msh, cl, level_set_function, min, mid);
    refine_interface2(msh, cl, level_set_function, mid, max);
}

template<typename T, size_t ET, typename Function>
void
refine_interface2(cuthho_mesh<T, ET>& msh, const Function& level_set_function, size_t levels)
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

        refine_interface2(msh, cl, level_set_function, 0, interface_points);
    }
}







/*****************************************************************************
 *   Test stuff ADDED by StePicca
 *****************************************************************************/

//template<typename T, typename Mesh>
//std::vector<size_t>
//subcell_finder<T,Mesh>(const Mesh& , const point<T,2>& , const typename Mesh::cell_type& , const mesh_init_params<T>&);

template<typename T, typename Mesh>
bool
pt_in_cell(const Mesh& msh, const point<T,2>& , const typename Mesh::cell_type& );

template<typename T,typename Mesh>
std::vector< point<T,2> >
equidistriduted_nodes(const Mesh& , const typename Mesh::cell_type& , size_t );

template<typename T,typename Mesh>
std::vector< point<T,2> >
equidistriduted_nodes_subcell(const Mesh& ,
          const typename Mesh::cell_type& ,
                              size_t , const std::vector<size_t>& ) ;


template<typename T>
std::vector<point<T,1> >
reference_nodes(size_t );

template< typename FonctionD , typename Mesh , typename FonctionA >
void
testing_level_set(const Mesh msh , const FonctionD& ,  const FonctionA& );

template< typename FonctionD , typename Mesh >
void
test_new_method(const Mesh , const FonctionD&  , const typename Mesh::cell_type & );


// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh , typename FonctionA >
void
gradient_checking1(const Mesh msh , const FonctionD& level_set_disc , const FonctionA& level_set_anal , const typename Mesh::cell_type & cl)
{
    typedef typename Mesh::point_type       point_type;
   
    double derD1x , derD1y , derAx , derAy , derD2x , derD2y ;
    Eigen::Matrix<double,2,1> derD1 ,derD2 , derA;
    point<double,2> node;

    
    auto pts = points(msh,cl);
    for(auto & node:pts)
    {
        derD1 = level_set_disc.gradient(node);
        derD2 = level_set_disc.gradient(node,msh,cl);
        derA = level_set_anal.gradient(node);
    
            derD1x = derD1(0);
            derD1y = derD1(1);
            derD2x = derD2(0);
            derD2y = derD2(1);
            derAx = derA(0);
            derAy = derA(1);

        
        /*
            if((derD1x-derD2x)>1e-2)
            {
            std::cout<<"Differnce between two x-evaluation system "<<(derD1x-derD2x)<<std::endl;
            }
        */
        
        if((derAx-derD2x)>1e-2)
        {
            std::cout<<"Differnce between analytic and NEW X-evaluation system "<<(derAx-derD2x)<<std::endl;
        }
   
        if((derAx-derD1x)>1e-2)
        {
            std::cout<<"Differnce between analytic and OLD X-evaluation system "<<(derAx-derD1x)<<std::endl;
        }
        
        
        /*
        if((derD1y-derD2y)>1e-2)
        {
            std::cout<<"Differnce between two y-evaluation system "<<(derD1y-derD2y)<<std::endl;
        }
         */
        
        
        if((derAy-derD2y)>1e-2)
        {
            std::cout<<"Differnce between analytic and NEW Y-evaluation system "<<(derAy-derD2y)<<std::endl;
        }
        
        if((derAy-derD1y)>1e-2)
        {
            std::cout<<"Differnce between analytic and OLD Y-evaluation system "<<(derAy-derD1y)<<std::endl;
        }
             
            
    }
        


}

// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh , typename FonctionA >
void
gradient_checking(const Mesh msh , const FonctionD& level_set_disc , const FonctionA& level_set_anal , const typename Mesh::cell_type & cl)
{
    typedef typename Mesh::point_type       point_type;
   
    double derD1x , derD1y , derAx , derAy , derD2x , derD2y ;
    Eigen::Matrix<double,2,1> derD1 ,derD2 , derA;
    point<double,2> node;

    auto pts = points(msh,cl);
    for(auto & node:pts)
    {
        derD1 = level_set_disc.gradient(node);
       // derD2 = level_set_disc.gradient(node,msh,cl);
        derA = level_set_anal.gradient(node);
    
            derD1x = derD1(0);
            derD1y = derD1(1);
         //   derD2x = derD2(0);
         //   derD2y = derD2(1);
            derAx = derA(0);
            derAy = derA(1);

        
        /*
            if((derD1x-derD2x)>1e-2)
            {
            std::cout<<"Differnce between two x-evaluation system "<<(derD1x-derD2x)<<std::endl;
            }
        */
        /*
        if((derAx-derD2x)>1e-2)
        {
            std::cout<<"Differnce between analytic and NEW X-evaluation system "<<(derAx-derD2x)<<std::endl;
        }
        */
        if((derAx-derD1x)>1e-2)
        {
            std::cout<<"Differnce between analytic and OLD X-evaluation system "<<(derAx-derD1x)<<std::endl;
        }
        
        
        /*
        if((derD1y-derD2y)>1e-2)
        {
            std::cout<<"Differnce between two y-evaluation system "<<(derD1y-derD2y)<<std::endl;
        }
         */
        
        /*
        if((derAy-derD2y)>1e-2)
        {
            std::cout<<"Differnce between analytic and NEW Y-evaluation system "<<(derAy-derD2y)<<std::endl;
        }
        */
        if((derAy-derD1y)>1e-2)
        {
            std::cout<<"Differnce between analytic and OLD Y-evaluation system "<<(derAy-derD1y)<<std::endl;
        }
             
            
    }
        


}


// Qualitative testing of the discrete level set function wrt the analytical one


template< typename FonctionD , typename Mesh , typename FonctionA >
void
testing_velocity(const Mesh msh , const FonctionD& vel_disc , const FonctionA& vel_anal)
{
    typedef typename Mesh::point_type       point_type;
    postprocess_output<double> postoutput1;
    
    Eigen::Matrix<double,2,1> valueA;
    point<double,2> node;
    size_t N, M;
    N = 40;
    M = 40;
    auto test_discx  = std::make_shared< gnuplot_output_object<double> >("vel_discX.dat");
    auto test_analx = std::make_shared< gnuplot_output_object<double> >("vel_analX.dat");
    auto test_discy  = std::make_shared< gnuplot_output_object<double> >("vel_discY.dat");
    auto test_analy = std::make_shared< gnuplot_output_object<double> >("vel_analY.dat");
    
    for(auto& cl:msh.cells)
    {
        auto pts = points(msh,cl);
        for(auto&pt : pts)
        {
            auto valueD = vel_disc(pt,msh,cl);
            valueA = vel_anal(pt);
            
            test_discx->add_data(pt,valueD.first);
            test_discy->add_data(pt,valueD.second);
            test_analx->add_data(pt,valueA(0));
            test_analy->add_data(pt,valueA(1));
            
            
        }
    }
   
    postoutput1.add_object(test_discx);
    postoutput1.add_object(test_analx);
    postoutput1.add_object(test_discy);
    postoutput1.add_object(test_analy);

    
    postoutput1.write();
    
}





template< typename FonctionD , typename Mesh , typename FonctionA >
void
testing_level_set(const Mesh msh , const FonctionD& level_set_disc , const FonctionA& level_set_anal)
{
    typedef typename Mesh::point_type       point_type;
    postprocess_output<double> postoutput1;
    double valueD , valueA , derDx , derDy , derAx , derAy;
    Eigen::Matrix<double,2,1> derD , derA;
    point<double,2> node;
    size_t N, M;
    N = 40;
    M = 40;
    auto test_disc  = std::make_shared< gnuplot_output_object<double> >("testing_interface_disc.dat");
    auto test_anal = std::make_shared< gnuplot_output_object<double> >("testing_interface_anal.dat");
    
    auto test_disc_gradX  = std::make_shared< gnuplot_output_object<double> >("testing_der_discX.dat");
    auto test_anal_gradX = std::make_shared< gnuplot_output_object<double> >("testing_der_analX.dat");
    
    auto test_disc_gradY  = std::make_shared< gnuplot_output_object<double> >("testing_der_discY.dat");
    auto test_anal_gradY = std::make_shared< gnuplot_output_object<double> >("testing_der_analY.dat");
    
    for(size_t i = 0 ; i<= N ; i++ )
    {
        for (size_t j = 0 ; j<= M ; j++ )
        {
            double px = i*(1.0/N);
            double py = j*(1.0/M);
            node = point_type(px,py);
            /*
            if (std::abs( level_set_disc(node)) <1e-2 ) {
                valueD = 1;
            }
            else
                valueD = 0;
            */
            valueD = level_set_disc(node);
            /*
            if (std::abs( level_set_anal(node)) <1e-2 ) {
                valueA = 1;
            }
            else
                valueA = 0;
            */
            valueA = level_set_anal(node);
            
            derD = level_set_disc.gradient(node);
            derA = level_set_anal.gradient(node);
        
            
            derDx = derD(0);
            derDy = derD(1);
            derAx = derA(0);
            derAy = derA(1);
            
            test_disc->add_data(node,valueD);
            test_anal->add_data(node,valueA);
            
            test_disc_gradX->add_data(node,derDx);
            test_anal_gradX->add_data(node,derAx);
            test_disc_gradY->add_data(node,derDy);
            test_anal_gradY->add_data(node,derAy);
            
        }
        
    }
    postoutput1.add_object(test_disc);
    postoutput1.add_object(test_anal);
    
    postoutput1.add_object(test_disc_gradX);
    postoutput1.add_object(test_anal_gradX);
    postoutput1.add_object(test_disc_gradY);
    postoutput1.add_object(test_anal_gradY);
    
    postoutput1.write();
    
}

// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh  >
void
testing_level_set2(const Mesh msh , const FonctionD& level_set_disc ,const FonctionD& level_set_fin)
{
    typedef typename Mesh::point_type       point_type;
    postprocess_output<double> postoutput1;
    double valueD , valueA , derDx , derDy , derAx , derAy;
    Eigen::Matrix<double,2,1> derD , derA;
    point<double,2> node;
    size_t N, M;
    N = 40;
    M = 40;
    auto test_disc  = std::make_shared< gnuplot_output_object<double> >("testing_interface_fin_disc.dat");
    auto test_anal = std::make_shared< gnuplot_output_object<double> >("testing_interface_fin_anal.dat");
    /*
    auto test_disc_gradX  = std::make_shared< gnuplot_output_object<double> >("testing_der_discX.dat");
    auto test_anal_gradX = std::make_shared< gnuplot_output_object<double> >("testing_der_analX.dat");
    
    auto test_disc_gradY  = std::make_shared< gnuplot_output_object<double> >("testing_der_discY.dat");
    auto test_anal_gradY = std::make_shared< gnuplot_output_object<double> >("testing_der_analY.dat");
    */
    for(size_t i = 0 ; i<= N ; i++ )
    {
        for (size_t j = 0 ; j<= M ; j++ )
        {
            double px = i*(1.0/N);
            double py = j*(1.0/M);
            node = point_type(px,py);
            /*
            if (std::abs( level_set_disc(node)) <1e-2 ) {
                valueD = 1;
            }
            else
                valueD = 0;
            */
            valueD = level_set_disc(node);
            /*
            if (std::abs( level_set_anal(node)) <1e-2 ) {
                valueA = 1;
            }
            else
                valueA = 0;
            */
            valueA = level_set_fin(node);
            /*
            derD = level_set_disc.gradient(node);
            derA = level_set_fin.gradient(node);
        
            
            derDx = derD(0);
            derDy = derD(1);
            derAx = derA(0);
            derAy = derA(1);
            */
            test_disc->add_data(node,valueD);
            test_anal->add_data(node,valueA);
            /*
            test_disc_gradX->add_data(node,derDx);
            test_anal_gradX->add_data(node,derAx);
            test_disc_gradY->add_data(node,derDy);
            test_anal_gradY->add_data(node,derAy);
            */
        }
        
    }
    postoutput1.add_object(test_disc);
    postoutput1.add_object(test_anal);
    /*
    postoutput1.add_object(test_disc_gradX);
    postoutput1.add_object(test_anal_gradX);
    postoutput1.add_object(test_disc_gradY);
    postoutput1.add_object(test_anal_gradY);
    */
    postoutput1.write();
    
}

// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh  >
void
testing_level_set_mesh2(const Mesh msh , const FonctionD& level_set_disc ,const FonctionD& level_set_fin)
{
    typedef typename Mesh::point_type       point_type;
    postprocess_output<double> postoutput1;
    double valueD , valueA , derDx , derDy , derAx , derAy;
    Eigen::Matrix<double,2,1> derD , derA;
    point<double,2> node;
    size_t N, M;
    N = 40;
    M = 40;
    auto test_disc  = std::make_shared< gnuplot_output_object<double> >("testing_interface_fin_disc_mesh2.dat");
    auto test_anal = std::make_shared< gnuplot_output_object<double> >("testing_interface_fin_anal_mesh2.dat");
    /*
    auto test_disc_gradX  = std::make_shared< gnuplot_output_object<double> >("testing_der_discX.dat");
    auto test_anal_gradX = std::make_shared< gnuplot_output_object<double> >("testing_der_analX.dat");
    
    auto test_disc_gradY  = std::make_shared< gnuplot_output_object<double> >("testing_der_discY.dat");
    auto test_anal_gradY = std::make_shared< gnuplot_output_object<double> >("testing_der_analY.dat");
    */
    for(size_t i = 0 ; i<= N ; i++ )
    {
        for (size_t j = 0 ; j<= M ; j++ )
        {
            double px = i*(1.0/N);
            double py = j*(1.0/M);
            node = point_type(px,py);
            /*
            if (std::abs( level_set_disc(node)) <1e-2 ) {
                valueD = 1;
            }
            else
                valueD = 0;
            */
            valueD = level_set_disc(node);
            /*
            if (std::abs( level_set_anal(node)) <1e-2 ) {
                valueA = 1;
            }
            else
                valueA = 0;
            */
            valueA = level_set_fin(node);
            /*
            derD = level_set_disc.gradient(node);
            derA = level_set_fin.gradient(node);
        
            
            derDx = derD(0);
            derDy = derD(1);
            derAx = derA(0);
            derAy = derA(1);
            */
            test_disc->add_data(node,valueD);
            test_anal->add_data(node,valueA);
            /*
            test_disc_gradX->add_data(node,derDx);
            test_anal_gradX->add_data(node,derAx);
            test_disc_gradY->add_data(node,derDy);
            test_anal_gradY->add_data(node,derAy);
            */
        }
        
    }
    postoutput1.add_object(test_disc);
    postoutput1.add_object(test_anal);
    /*
    postoutput1.add_object(test_disc_gradX);
    postoutput1.add_object(test_anal_gradX);
    postoutput1.add_object(test_disc_gradY);
    postoutput1.add_object(test_anal_gradY);
    */
    postoutput1.write();
    
}

// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh  >
void
testing_level_set_mesh1(const Mesh msh , const FonctionD& level_set_disc ,const FonctionD& level_set_fin)
{
    typedef typename Mesh::point_type       point_type;
    postprocess_output<double> postoutput1;
    double valueD , valueA , derDx , derDy , derAx , derAy;
    Eigen::Matrix<double,2,1> derD , derA;
    point<double,2> node;
    size_t N, M;
    N = 40;
    M = 40;
    auto test_disc  = std::make_shared< gnuplot_output_object<double> >("testing_interface_fin_disc_mesh1.dat");
    auto test_anal = std::make_shared< gnuplot_output_object<double> >("testing_interface_fin_anal_mesh1.dat");
    /*
    auto test_disc_gradX  = std::make_shared< gnuplot_output_object<double> >("testing_der_discX.dat");
    auto test_anal_gradX = std::make_shared< gnuplot_output_object<double> >("testing_der_analX.dat");
    
    auto test_disc_gradY  = std::make_shared< gnuplot_output_object<double> >("testing_der_discY.dat");
    auto test_anal_gradY = std::make_shared< gnuplot_output_object<double> >("testing_der_analY.dat");
    */
    for(size_t i = 0 ; i<= N ; i++ )
    {
        for (size_t j = 0 ; j<= M ; j++ )
        {
            double px = i*(1.0/N);
            double py = j*(1.0/M);
            node = point_type(px,py);
            /*
            if (std::abs( level_set_disc(node)) <1e-2 ) {
                valueD = 1;
            }
            else
                valueD = 0;
            */
            valueD = level_set_disc(node);
            /*
            if (std::abs( level_set_anal(node)) <1e-2 ) {
                valueA = 1;
            }
            else
                valueA = 0;
            */
            valueA = level_set_fin(node);
            /*
            derD = level_set_disc.gradient(node);
            derA = level_set_fin.gradient(node);
        
            
            derDx = derD(0);
            derDy = derD(1);
            derAx = derA(0);
            derAy = derA(1);
            */
            test_disc->add_data(node,valueD);
            test_anal->add_data(node,valueA);
            /*
            test_disc_gradX->add_data(node,derDx);
            test_anal_gradX->add_data(node,derAx);
            test_disc_gradY->add_data(node,derDy);
            test_anal_gradY->add_data(node,derAy);
            */
        }
        
    }
    postoutput1.add_object(test_disc);
    postoutput1.add_object(test_anal);
    /*
    postoutput1.add_object(test_disc_gradX);
    postoutput1.add_object(test_anal_gradX);
    postoutput1.add_object(test_disc_gradY);
    postoutput1.add_object(test_anal_gradY);
    */
    postoutput1.write();
    
}

// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh  >
void
testing_level_set_mesh0(const Mesh msh , const FonctionD& level_set_disc ,const FonctionD& level_set_fin)
{
    typedef typename Mesh::point_type       point_type;
    postprocess_output<double> postoutput1;
    double valueD , valueA , derDx , derDy , derAx , derAy;
    Eigen::Matrix<double,2,1> derD , derA;
    point<double,2> node;
    size_t N, M;
    N = 40;
    M = 40;
    auto test_disc  = std::make_shared< gnuplot_output_object<double> >("testing_interface_fin_disc_mesh0.dat");
    auto test_anal = std::make_shared< gnuplot_output_object<double> >("testing_interface_fin_anal_mesh0.dat");
    /*
    auto test_disc_gradX  = std::make_shared< gnuplot_output_object<double> >("testing_der_discX.dat");
    auto test_anal_gradX = std::make_shared< gnuplot_output_object<double> >("testing_der_analX.dat");
    
    auto test_disc_gradY  = std::make_shared< gnuplot_output_object<double> >("testing_der_discY.dat");
    auto test_anal_gradY = std::make_shared< gnuplot_output_object<double> >("testing_der_analY.dat");
    */
    for(size_t i = 0 ; i<= N ; i++ )
    {
        for (size_t j = 0 ; j<= M ; j++ )
        {
            double px = i*(1.0/N);
            double py = j*(1.0/M);
            node = point_type(px,py);
            /*
            if (std::abs( level_set_disc(node)) <1e-2 ) {
                valueD = 1;
            }
            else
                valueD = 0;
            */
            valueD = level_set_disc(node);
            /*
            if (std::abs( level_set_anal(node)) <1e-2 ) {
                valueA = 1;
            }
            else
                valueA = 0;
            */
            valueA = level_set_fin(node);
            /*
            derD = level_set_disc.gradient(node);
            derA = level_set_fin.gradient(node);
        
            
            derDx = derD(0);
            derDy = derD(1);
            derAx = derA(0);
            derAy = derA(1);
            */
            test_disc->add_data(node,valueD);
            test_anal->add_data(node,valueA);
            /*
            test_disc_gradX->add_data(node,derDx);
            test_anal_gradX->add_data(node,derAx);
            test_disc_gradY->add_data(node,derDy);
            test_anal_gradY->add_data(node,derAy);
            */
        }
        
    }
    postoutput1.add_object(test_disc);
    postoutput1.add_object(test_anal);
    /*
    postoutput1.add_object(test_disc_gradX);
    postoutput1.add_object(test_anal_gradX);
    postoutput1.add_object(test_disc_gradY);
    postoutput1.add_object(test_anal_gradY);
    */
    postoutput1.write();
    
}
// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh , typename FonctionA >
void
test_new_method(const Mesh msh , const FonctionD& level_set_disc , const FonctionA& level_set_anal , const typename Mesh::cell_type & cl)
{
    typedef typename Mesh::point_type       point_type;
    double valueD1 , valueD2 ,valueA ;
    
    auto pts =points(msh,cl);
    
    
    
    
    timecounter tc1;
    tc1.tic();
    for(auto& pt : pts)
    {
        valueD1 = level_set_disc(pt,msh,cl);
       // valueA = level_set_anal(pt);
       // valueD2 = level_set_disc(pt);
       // std::cout<<"Differnce between two evaluation system "<<(valueD1-valueD2)<<std::endl;
      //  std::cout<<"Im in test_new_method";
    }
    
   
    
    
    tc1.toc();
    std::cout << bold << yellow << "Time for the new evaluation: " << tc1 << " seconds" << reset << std::endl;
    
   tc1.tic();
       for(auto& pt : pts)
       {
           valueD2 = level_set_disc(pt);
         //  std::cout<<"Im in test_new_method";
       }
    tc1.toc();
    std::cout << bold << yellow << "Time for the old evaluation: " << tc1 << " seconds" << reset << std::endl;
    
    tc1.tic();
         for(auto& pt : pts)
         {
             valueA = level_set_anal(pt);
         }
      tc1.toc();
      std::cout << bold << yellow << "Time for the analytic evaluation: " << tc1 << " seconds" << reset << std::endl;
    
    
}

// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh >
void
time_NEWdiscrete_testing(const Mesh msh , const FonctionD& level_set_disc , const typename Mesh::cell_type & cl)
{
    typedef typename Mesh::point_type       point_type;
    double valueD1 , valueD2 ,valueA ;
    
    auto pts =points(msh,cl);
    
  
    for(auto& pt : pts)
    {
        valueD1 = level_set_disc(pt,msh,cl);
       
    }
    
    
}

// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh >
void
time_OLDdiscrete_testing(const Mesh msh , const FonctionD& level_set_disc , const typename Mesh::cell_type & cl)
{
    typedef typename Mesh::point_type       point_type;
    double valueD1 , valueD2 ,valueA ;
    
    auto pts =points(msh,cl);
    
  
    for(auto& pt : pts)
    {
        valueD1 = level_set_disc(pt);
       
    }
    
    
}

template< typename FonctionA , typename Mesh >
void
time_analytic_testing(const Mesh msh , const FonctionA& level_set_anal , const typename Mesh::cell_type & cl)
{
    typedef typename Mesh::point_type       point_type;
    double valueD1 , valueD2 ,valueA ;
    
    auto pts =points(msh,cl);
    
  
    for(auto& pt : pts)
    {
        valueA = level_set_anal(pt);
       
    }
    
    
}


template< typename FonctionD , typename Mesh >
void
time_face_testing(const Mesh msh , const FonctionD& level_set_disc , const typename Mesh::face_type & fc)
{
    
    typedef typename Mesh::point_type       point_type;
    double valueD1 , valueD2 ,valueA ;
       
    auto pts =points(msh,fc);
       
     
    for(auto& pt : pts)
    {
        valueD1 = level_set_disc(pt,msh,fc);
          
    }
       
    
    
}

template< typename FonctionA , typename Mesh >
void
time_faceANALITIC_testing(const Mesh msh , const FonctionA& level_set_anal , const typename Mesh::face_type & fc)
{
    
    typedef typename Mesh::point_type       point_type;
    double valueD1 , valueD2 ,valueA ;
       
    auto pts =points(msh,fc);
       
     
    for(auto& pt : pts)
    {
        valueA = level_set_anal(pt);
          
    }
       
    
    
}



// Qualitative testing of the discrete level set function wrt the analytical one
template< typename FonctionD , typename Mesh , typename FonctionA >
void
test_new_method(const Mesh msh , const FonctionD& level_set_disc , const FonctionA& level_set_anal , const typename Mesh::cell_type & cl , const typename Mesh::face_type & fc)
{
    typedef typename Mesh::point_type       point_type;
    double valueD1 , valueD2 ,valueA , valueD3;
    
    auto pts =points(msh,fc);
    
    
    
    
    timecounter tc1;
    tc1.tic();
    for(auto& pt : pts)
    {
        valueD1 = level_set_disc(pt,msh,fc);
        valueA = level_set_anal(pt);
        //valueD2 = level_set_disc(pt,msh,cl);
        valueD3 = level_set_disc(pt);
        
        std::cout<<"Differnce between FACE and OLD evaluation system "<<(valueD1-valueD3)<<std::endl;
        std::cout<<"Error between analytic and face evaluation "<<(valueD1-valueA)<<std::endl;
    }
    

    
}



template<typename T>
std::vector<point<T,1> >
reference_nodes_ordered(size_t degree)
{
    auto comp_degree = degree + 1;
    size_t reqd_nodes = comp_degree;

    std::vector<point<T,1> > ret;
    ret.reserve(reqd_nodes);

    point<T,1>  qp , qp_1;
    T           a1, a2;
    T           delta_x;
    switch(reqd_nodes)
    {
        case 1:
            qp = point<T,1>({0.0});
            ret.push_back(qp);
            return ret;

        case 2:
            qp = point<T,1>({ 1.0 });
            ret.push_back( -qp );
            ret.push_back( qp );
            return ret;

        case 3:
            qp = point<T,1>({ 1.0 });
            ret.push_back( -qp );
            ret.push_back( qp );
            qp = point<T,1>({0.0});
            ret.push_back( qp );
            return ret;

        case 4:
            a1 = 1.0/3.0;
            qp = point<T,1>({ 1.0 });
            ret.push_back( -qp );
            ret.push_back( qp );
            qp = point<T,1>({ a1 });
            ret.push_back( -qp );
            ret.push_back( qp );
            return ret;

        case 5:
            // Be carefull in what order data is inserted in ret!
            // In Gauss Legendre the first one was 0.0, now is the last one
            a2 = 0.5;
            a1 = 1.0;
            qp = point<T,1>({ a1 });
            ret.push_back( -qp );
            ret.push_back( qp );

            qp = point<T,1>({ a2 });
            ret.push_back( -qp );
            qp_1 = point<T,1>({ 0.0 });
            ret.push_back( qp_1 );
            ret.push_back( qp );

         //   qp = point<T,1>({ 0.0 });
         //   ret.push_back( qp );
            return ret;

        default:

            delta_x = 2.0/degree;
            a1 = 1.0;
            while (a1>1e-10) {
                qp = point<T,1>({ a1 });
                ret.push_back( -qp );
                ret.push_back( qp );
                a1-=delta_x;

            }
            if(a1<1e-10 && a1>-1e-10)
            {
                qp = point<T,1>({0.0});
                ret.push_back( qp );
            }
            std::sort(ret.begin()+2, ret.end() ,[](point<T,1>  a, point<T,1>  b) {
                return a.x() < b.x();
            } );
            return ret;
    }
    return ret;
}


// Lagrangian basis b_kl(x,y) = b_k(x)*b_l(y) over a set of equidistributed 2-dimensional nodes (3D CASE NOT YET IMPLEMENTED)

template<typename T,typename Mesh>
std::vector< point<T,2> >
equidistriduted_nodes_ordered_bis(const Mesh& msh,
          const typename Mesh::cell_type& cl,
          size_t degree)
{
    typedef typename Mesh::point_type    point_type;

    auto qps = reference_nodes_ordered<T>(degree); //Ordering:  0 2 3 4 5 6 ... 1
   
    // 3 -  12 - 11 - 10 - 2
    // 13 - 22 - 21 - 20 - 9
    // 14 - 23 - 24 - 19 - 8
    // 15 - 16 - 17 - 18 - 7
    // 0 -  4 -  5  - 6  - 1

    auto pts = points(msh, cl);
    
    auto v0 = pts[1] - pts[0];
    auto v1 = pts[2] - pts[1];
    auto v2 = pts[3] - pts[2];
    auto v3 = pts[3] - pts[0];

    std::vector< point<T,2> > ret((degree+1)*(degree+1));

    auto P = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].x() * (1-xi)*(1-eta) +
               0.25 * pts[1].x() * (1+xi)*(1-eta) +
               0.25 * pts[2].x() * (1+xi)*(1+eta) +
               0.25 * pts[3].x() * (1-xi)*(1+eta);
    };

    auto Q = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].y() * (1-xi)*(1-eta) +
               0.25 * pts[1].y() * (1+xi)*(1-eta) +
               0.25 * pts[2].y() * (1+xi)*(1+eta) +
               0.25 * pts[3].y() * (1-xi)*(1+eta);
    };

    /// ADDING VERTICES:
    
    // (-1,-1)
    auto qp_x = qps[0];
    auto qp_y = qps[0];
    auto xi = qp_x.x();
    auto eta = qp_y.x();
    auto px = P(xi, eta);
    auto py = Q(xi, eta);
    ret[0] = ( point_type(px, py) );
    // (1,-1)
    qp_x = qps[1];
    qp_y = qps[0];
    xi = qp_x.x();
    eta = qp_y.x();
    px = P(xi, eta);
    py = Q(xi, eta);
    ret[1] = ( point_type(px, py) );
    // (1,1)
    qp_x = qps[1];
    qp_y = qps[1];
    xi = qp_x.x();
    eta = qp_y.x();
    px = P(xi, eta);
    py = Q(xi, eta);
    ret[2] = ( point_type(px, py) );
    // (-1,1)
    qp_x = qps[0];
    qp_y = qps[1];
    xi = qp_x.x();
    eta = qp_y.x();
    px = P(xi, eta);
    py = Q(xi, eta);
    ret[3] = ( point_type(px, py) );
    
    /// Counter for each side of the 2D - square : starting to count from a vertice, the position where I save the point in ret is (degree -1)*j , with j the j-esima face.
    
    int  count0 = 4 , count1 = 4 + degree - 1 , count2 = 4 + 2*(degree - 1), count3 = 4 + 3*(degree - 1)  ;
    
    int count_bis0 = 4*degree ; // counter initialisation (USELESS)
    int j_max = floor((degree-1)/2) ; // number of internal layour of points
  //  std::cout<<"j max "<<j_max<<std::endl;
    int pos_right = 100; // initial point from the right from bottom -> up
    int pos_left = 100; // initial point from the left from left -> right
   // size_t i_min = 0;
 //   std::cout<<"inizia ciclo"<<std::endl;
    for (int j = 0 ; j <= j_max ; j++ ) { // for each layout of point
       // bool i_loop = FALSE;
        for(int i = std::max(0,j-1) ; i < degree -1 - j ; i++ ) // I move from into the points over a side of each layout
        {
         
            if( i == std::max(0,j-1) && j > 0) // vertices
            {
                // different pos_left depending on the layout. Especially this rules the y starting point
                if(j == 0)
                    pos_left = 0;
                else if(j == 1)
                    pos_left = 2;
                else
                    pos_left = 2+(j-1);
                           
                //   std::cout<<"pos_left "<<pos_left<<std::endl;
                qp_x = qps[2+i];
                qp_y = qps[pos_left]; //qps[0 + 2*j];
                //   std::cout<<"qp_y "<<qp_y<<std::endl;
                xi = qp_x.x();
                eta = qp_y.x();
                px = P(xi, eta);
                py = Q(xi, eta);
                // Here change counters. No more count0, count1 etc.. Just count_bis0, to have the vertices enumerate one after the other.
                count_bis0 = count0; // counter_bis0 re-initialisation for each loop (first loop)
                //std::cout<<"count bis 0 "<<count_bis0<<std::endl;
                ret[count_bis0] = ( point_type(px, py) ); // adding point bottom layout
                //std::cout<<"count0 is "<<count0<<" , pt0"<<point_type(px, py)<<std::endl;
                //std::cout<<ret[count0]<<std::endl;
                count_bis0++;
                
                if(j==0)
                    pos_right = 1;
                else
                    pos_right = degree + 1 - j;
                // size_t pos = (1 + j*(degree-1));
                // if(pos>degree)
                //     pos -= j*degree ;
                           
                qp_x =  qps[pos_right];
                qp_y = qps[2 + i];
                xi = qp_x.x();
                eta = qp_y.x();
                px = P(xi, eta);
                py = Q(xi, eta);
                ret[count_bis0] = ( point_type(px, py) ); // adding point right layout in count_bis0 pos.
                //std::cout<<"count1 is "<<count1<<" , pt1"<<point_type(px, py)<<std::endl;
                //std::cout<<"count_bis0 is "<<count_bis0<<std::endl;
                //std::cout<<ret[count_bis0]<<std::endl;
                count_bis0++;
                
                qp_x = qps[degree - i];
                qp_y =  qps[pos_right] ;
                xi = qp_x.x();
                eta = qp_y.x();
                px = P(xi, eta);
                py = Q(xi, eta);
                ret[count_bis0] = ( point_type(px, py) ); // adding point top layout in count_bis0 pos.
                //std::cout<<"count_bis0 is "<<count_bis0<<" and pt "<<ret[count_bis0]<<std::endl;
                //std::cout<<"count2 is "<<count2<<" , pt2"<<point_type(px, py)<<std::endl;
                
                count_bis0++;
                qp_x = qps[pos_left];
                qp_y = qps[degree - i];
                xi = qp_x.x();
                eta = qp_y.x();
                px = P(xi, eta);
                py = Q(xi, eta);
                ret[count_bis0] = ( point_type(px, py) ); // adding point left layout in count_bis0 pos.
                //std::cout<<"count_bis0 is "<<count_bis0<<" and pt "<<ret[count_bis0]<<std::endl;
                //std::cout<<"count3 is "<<count3<<" , pt3"<<point_type(px, py)<<std::endl;
                count_bis0++;
                
                // Uploading counters according to count_bis0. I start enumerate "face" unknown for the current layout (IN ELSE HERE BELOW) from the last position count_bis0.
                count0 = count_bis0 ;
                count1 = count0 + std::ceil((degree - 2.0*j)/2.0); //  count0 + (degree - 2); it was
                count2 = count1 + std::ceil((degree - 2.0*j)/2.0);
                count3 = count2 + std::ceil((degree - 2.0*j)/2.0);

            }
            else // NOT vertices -> node in the sides of each layout
            {
                
                //   std::cout<<"i "<<i<<" j "<<j<<std::endl;
                // vertical position where starting for each bottom side layout.
                if(j == 0)
                    pos_left = 0;
                else if(j == 1)
                    pos_left = 2;
                else
                    pos_left = 2+(j-1);
                
                //   std::cout<<"pos_left "<<pos_left<<std::endl;
                qp_x = qps[2+i];
                qp_y = qps[pos_left]; //qps[0 + 2*j];
                //   std::cout<<"qp_y "<<qp_y<<std::endl;
                xi = qp_x.x();
                eta = qp_y.x();
                px = P(xi, eta);
                py = Q(xi, eta);
                ret[count0] = ( point_type(px, py) ); // from left, node of each bottom side layout
                //std::cout<<"count0 is "<<count0<<" , pt0"<<point_type(px, py)<<std::endl;
                //std::cout<<ret[count0]<<std::endl;
                count0++;
                // x-position where to start to increase to get points for each right side layout of points
                if(j==0)
                    pos_right = 1;
                else
                    pos_right = degree + 1 - j;
                // size_t pos = (1 + j*(degree-1));
                // if(pos>degree)
                //     pos -= j*degree ;
            
                qp_x =  qps[pos_right];
                qp_y = qps[2 + i];
                xi = qp_x.x();
                eta = qp_y.x();
                px = P(xi, eta);
                py = Q(xi, eta);
                ret[count1] = ( point_type(px, py) ); // using count1 that allows correct enumeration
                //std::cout<<"count1 is "<<count1<<" , pt1"<<point_type(px, py)<<std::endl;
                //std::cout<<ret[count1]<<std::endl;
                count1++; // count1 bigger than count0. Count1 just enumerate right faces' layout
                
                qp_x = qps[degree - i];
                qp_y =  qps[pos_right] ;
                xi = qp_x.x();
                eta = qp_y.x();
                px = P(xi, eta);
                py = Q(xi, eta);
                ret[count2] = ( point_type(px, py) );
                //std::cout<<"count2 is "<<count2<<" , pt2"<<point_type(px, py)<<std::endl;
                count2++; // count2 just enumerate top faces' layout
                
                qp_x = qps[pos_left];
                qp_y = qps[degree - i];
                xi = qp_x.x();
                eta = qp_y.x();
                px = P(xi, eta);
                py = Q(xi, eta);
                ret[count3] = ( point_type(px, py) );
                //std::cout<<"count3 is "<<count3<<" , pt3"<<point_type(px, py)<<std::endl;
                count3++; // count3 just enumerate left faces' layout
            
                
            }
        }
        // Uploading for the next layout.
        count0 = count3 ;
        count1 = count0 + (degree - 2*(j+1)); //  count0 + (degree - 2); it was
        count2 = count1 + (degree - 2*(j+1));
        count3 = count2 + (degree - 2*(j+1));
        //}
    }
    
    /// Middle point --> the internal node is treated a part of the others. Just even degrees have it.
    if( degree % 2 == 0)
    {
        qp_x = qps[degree - floor((degree-1)/2)];
        qp_y = qps[degree - floor((degree-1)/2)];
        xi = qp_x.x();
        eta = qp_y.x();
        px = P(xi, eta);
        py = Q(xi, eta);
        ret[count0] = ( point_type(px, py) );
        //std::cout<<"counto MIDDLE is "<<count0<<" , pt3"<<point_type(px, py)<<std::endl;
    }
    return ret;
}

template<typename Mesh, typename T >
class cell_basis_Lagrangian_ordered
{
    typedef typename Mesh::coordinate_type  coordinate_type;
    typedef typename Mesh::point_type       point_type;

    point_type          cell_bar;
    coordinate_type     cell_h;
    size_t              basis_degree, basis_size;
    std::vector<point<T, 2> >          nodes;
    std::vector<size_t>         row_indeces , col_indeces ;
    
public:
    cell_basis_Lagrangian_ordered(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
    {
        //nodes           = equidistriduted_nodes<T,Mesh>(msh, cl, degree);
        // nodes           = equidistriduted_nodes_ordered<T,Mesh>(msh, cl, degree);
        nodes           = equidistriduted_nodes_ordered_bis<T,Mesh>(msh, cl, degree);
        basis_degree    = degree;
        basis_size      = (basis_degree+1)*(basis_degree+1);
        set_col_row_indeces();
    
    }
    
    void set_col_row_indeces()
    {
        //std::cout<<"SONO QUA NEL SET INDECES "<<row_indeces.size()<<" and "<<col_indeces.size()<<std::endl;
        // It could be cool to just save once and not call each time. To do that or call the  variable Static(?) or give it just as external fx and initiliaze once.
        if(row_indeces.size()==0 && col_indeces.size()==0 )
        {
            row_indeces.push_back(0);
            row_indeces.push_back(1);
            col_indeces.push_back(0);
            col_indeces.push_back(3);
            for(size_t i = 0 ; i<basis_degree-1 ; i++){
                row_indeces.push_back(4+i);
                col_indeces.push_back(3*basis_degree+1);
            }
            
        }
    }
        
    
    /*
    cell_basis_Lagrangian(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, const std::vector<size_t>& indeces)
    {
        
        nodes = equidistriduted_nodes_subcell<T,Mesh>(msh, cl, degree, indeces);
        basis_degree    = degree;
        basis_size      = (basis_degree+1)*(basis_degree+1);
        
    }
    */

    Matrix<T, Dynamic, 1>
    eval_basis(const point_type& pt)
    {
        Matrix<T, Dynamic, 1> rety = Matrix<T, Dynamic, 1>::Zero(basis_degree+1);
        Matrix<T, Dynamic, 1> retx = Matrix<T, Dynamic, 1>::Zero(basis_degree+1);
        Matrix<T, Dynamic, 1> ret =  Matrix<T, Dynamic, 1>::Zero(basis_size);

    
        /*
        for(size_t i = 0 ; i<basis_degree+1 ; i++){
            std::cout<<row_indeces[i]<<" , ";
        }
        std::cout<<std::endl;
        for(size_t i = 0 ; i<basis_degree+1 ; i++){
            std::cout<<col_indeces[i]<<" , ";
        }
        std::cout<<std::endl;
        */
        // the indeces of first row and column are already saved. Just need a loop over them and then do the right tensor product maintaining the ordering of the nodes.
        size_t ct = 0;
        for(auto& k : row_indeces)
        {
            T bk = 1;
            for ( auto& j : row_indeces ){
                if(j!=k){
            bk *= ( ( pt.x() - (nodes.at(j)).x() )/ ( (nodes.at(k)).x() - (nodes.at(j)).x() ) );
                }
            }
            
            retx(ct) = bk;
            ct++;
        }
        
        ct = 0;
        for(auto& l : col_indeces)
        {
            T bl = 1;
            for ( auto& j : col_indeces ){
                if(j!=l){
            bl *= ( ( pt.y() - (nodes.at(j)).y() )/ ( (nodes.at(l)).y() - (nodes.at(j)).y() ) );
                }
            }
            rety(ct) = bl;
            ct++;
        }
        
    
    
        // VERTICES FUNCTIONS:
        if( basis_degree == 0){
            ret(0) = rety(0)*retx(0);
            return ret ;
        }
        
        
        ret(0) = rety(0)*retx(0);
        ret(1) = rety(0)*retx(1);
        ret(2) = rety(1)*retx(1);
        ret(3) = rety(1)*retx(0);
        
        if( basis_degree == 1)
            return ret ;
        
        int  count0 = 4 , count1 = 4 + basis_degree - 1 , count2 = 4 + 2*(basis_degree - 1), count3 = 4 + 3*(basis_degree - 1)  ;
        int count_bis0 = 4*basis_degree ; // counter initialisation (USELESS)
        int j_max = floor((basis_degree-1)/2) ; // number of internal layour of points
        int pos_right = 100; // initial point from the right from bottom -> up
        int pos_left = 100; // initial point from the left from left -> right
     
        for (int j = 0 ; j <= j_max ; j++ ) { // for each layout of point
            // bool i_loop = FALSE;
            for(int i = std::max(0,j-1) ; i < basis_degree -1 - j ; i++ ) // I move from left into the points over a side of each layout
            {
                if( i == std::max(0,j-1) && j > 0) // vertices
                {
                    // different pos_left depending on the layout. Especially this rules the y starting point
                    if(j == 0)
                        pos_left = 0;
                    else if(j == 1)
                        pos_left = 2;
                    else
                        pos_left = 2+(j-1);
                                             
                    // Here change counters. No more count0, count1 etc.. Just count_bis0, to have the vertices enumerate one after the other.
                    count_bis0 = count0; // counter_bis0 re-initialisation for each loop (first loop)
                    //std::cout<<"count bis 0 "<<count_bis0<<std::endl;
                    ret(count_bis0) = ( retx(2+i)*rety(pos_left) ); // adding point bottom layout
                    count_bis0++;
                    if(j==0)
                        pos_right = 1;
                    else
                        pos_right = basis_degree + 1 - j;
    
                    ret(count_bis0) = ( retx(pos_right)*rety(2 + i) ); // adding point right layout in count_bis0 pos.
                    count_bis0++;

                    ret(count_bis0) = ( retx(basis_degree - i)*rety(pos_right) ); // adding point top layout in count_bis0 pos.
                                 
                    count_bis0++;
                                 
                    ret(count_bis0) = ( retx(pos_left)*rety(basis_degree - i) );
                                 
                    count_bis0++;
                                  
                    // Uploading counters according to count_bis0. I start enumerate "face" unknown for the current layout (IN ELSE HERE BELOW) from the last position count_bis0.
                    count0 = count_bis0 ;
                    count1 = count0 + std::ceil((basis_degree - 2.0*j)/2.0); //  count0 + (degree - 2); it was
                    count2 = count1 + std::ceil((basis_degree - 2.0*j)/2.0);
                    count3 = count2 + std::ceil((basis_degree - 2.0*j)/2.0);

                }
                else // NOT vertices -> node in the sides of each layout
                {
                    //   std::cout<<"i "<<i<<" j "<<j<<std::endl;
                    // vertical position where starting for each bottom side layout.
                    if(j == 0)
                        pos_left = 0;
                    else if(j == 1)
                        pos_left = 2;
                    else
                        pos_left = 2+(j-1);
                                
                    ret(count0) = ( retx(2+i)*rety(pos_left) );
                    count0++;
                    // x-position where to start to increase to get points for each right side layout of points
                    if(j==0)
                        pos_right = 1;
                    else
                        pos_right = basis_degree + 1 - j;
                  
                    ret(count1) = ( retx(pos_right)*rety(2 + i) );
                    count1++; // count1 bigger than count0. Count1 just enumerate right faces' layout
              
                    ret(count2) = ( retx(basis_degree - i)*rety(pos_right) );
                    count2++; // count2 just enumerate top faces' layout
                                  
                    ret(count3) = ( retx(pos_left)*rety(basis_degree - i) );
                    count3++; // count3 just enumerate left faces' layout
                              
                                  
                }
            }
            // Uploading for the next layout.
            count0 = count3 ;
            count1 = count0 + (basis_degree - 2*(j+1));
            count2 = count1 + (basis_degree - 2*(j+1));
            count3 = count2 + (basis_degree - 2*(j+1));
                     
        }
                      
        /// Middle point --> the internal node is treated a part of the others. Just even degrees have it.
        if( basis_degree % 2 == 0)
            ret(count0) = (retx(basis_degree - floor((basis_degree-1)/2))*rety(basis_degree - floor((basis_degree-1)/2)));
                          
        return ret;
        

    }


    
    Matrix<T, Dynamic, 2>
    eval_gradients(const point_type& pt)
    {
        // Modified Yves Daoust Algorithm (https://math.stackexchange.com/questions/809927/first-derivative-of-lagrange-polynomial)
        
        Matrix<T, Dynamic, 2> ret = Matrix<T, Dynamic, 2>::Zero(basis_size, 2);
       
        Matrix<T, Dynamic, 1> rety = Matrix<T, Dynamic, 1>::Zero(basis_size);
        Matrix<T, Dynamic, 1> retx = Matrix<T, Dynamic, 1>::Zero(basis_size);
        Matrix<T, Dynamic, 1> sy = Matrix<T, Dynamic, 1>::Zero(basis_size);
        Matrix<T, Dynamic, 1> sx = Matrix<T, Dynamic, 1>::Zero(basis_size);


        // for each l, b_l(y)' = {sum(tmpy!=l)[prod(jy!=l,jy!=tmpy)[x-x_jy]]}/prod(tmpy!=l)[x_l-x_tmpy]
        
        for ( size_t l = 0; l < basis_size ; l++ )
        {
            size_t col = l%(basis_degree+1);
            T bl = 1.0 , bl_der = 1.0 ;
            T sumy = 0.0;
            for (size_t tmpy = col; tmpy <= col+(basis_degree+1)*basis_degree; tmpy+=(basis_degree+1))
            {
                 T sumyy = 1.0 ;
                if (tmpy!=l)
                {

                    bl *= ( ( pt.y() - (nodes.at(tmpy)).y() )/ ( (nodes.at(l)).y() - (nodes.at(tmpy)).y() ) );

                    bl_der *= ( 1.0/ ( (nodes.at(l)).y() - (nodes.at(tmpy)).y() ) );
                    for ( size_t jy = col; jy <= col+(basis_degree+1)*basis_degree; jy+=(basis_degree+1) )
                   {
                        if (jy!=tmpy && jy!=l)
                        {
                            sumyy *= ( pt.y()-(nodes.at(jy)).y() );
                        }
                    }
                    sumy +=sumyy;
                }
            }
            rety(l) = bl;
            sy(l) = bl_der*sumy;
        }
        
        // For the x-derivative of b_k(x), same procedure of b_l(y)'
         
        for (size_t k = 0 ; k < basis_size ; k++ )
        {
            size_t row=floor( k/(basis_degree+1) );
            T bk = 1.0 , bk_der = 1.0 ;
            T sumx = 0.0;
            for (size_t tmpx = (basis_degree+1)*row; tmpx <= (basis_degree+1)*(row+1)-1; tmpx++)
            {
                T sumxx = 1.0 ;
                
                if (tmpx!=k) {
                    
                    bk *= ( ( pt.x() - (nodes.at(tmpx)).x() )/ ( (nodes.at(k)).x() - (nodes.at(tmpx)).x() ) );
                    bk_der *= ( 1.0/ ( (nodes.at(k)).x() - (nodes.at(tmpx)).x() ) );
                    for (size_t jx = (basis_degree+1)*row; jx <= (basis_degree+1)*(row+1)-1; jx++)
                    {
                        if (jx!=tmpx && jx!=k)
                        {
                            sumxx *= ( pt.x()-(nodes.at(jx)).x() );
                        }
                    }
                    sumx += sumxx;
                    
                }
            }
            retx(k) = bk;
            sx(k) = bk_der*sumx;
        }
        
        for (size_t i = 0; i<basis_size; i++)
        {
            ret(i,0) = rety(i)*sx(i);
            ret(i,1) = retx(i)*sy(i);
            
        }
        return ret;
        
    }

    size_t size() const
    {
        return basis_size;
    }

    size_t degree() const
    {
        return basis_degree;
    }

    static size_t size(size_t degree)
    {
        return (degree+1)*(degree+1);
    }
};


int binomial_coeff_fx(int n , int k) {

        int C[n + 1][k + 1];
        int i, j;
      
        // Caculate value of Binomial Coefficient
        // in bottom up manner
        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= std::min(i, k); j++)
            {
                // Base Cases
                if (j == 0 || j == i)
                    C[i][j] = 1;
      
                // Calculate value using previously
                // stored values
                else
                    C[i][j] = C[i - 1][j - 1] +
                              C[i - 1][j];
            }
        }
      
        return C[n][k];
}



template<typename Mesh, typename VT>
class cell_basis_Bernstein
{
    typedef typename Mesh::coordinate_type  coordinate_type;
    typedef typename Mesh::point_type       point_type;

    point_type          cell_bar;
    coordinate_type     cell_h;
    size_t              basis_degree, basis_size;
    coordinate_type min_x , max_x , min_y , max_y ;
    coordinate_type scaling_x , scaling_y ;

#ifdef POWER_CACHE
    std::vector<coordinate_type>  power_cache , power_cache_bis;
    std::vector<size_t>  binomial_coeff , binomial_coeff_bis ;
#endif

public:
    cell_basis_Bernstein(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
    {
        //cell_bar        = barycenter(msh, cl);
        //cell_h          = diameter(msh, cl);
        basis_degree    = degree;
        basis_size      = (basis_degree+1)*(basis_degree+1);
        min_x = points(msh,cl)[0].x();
        max_x = points(msh,cl)[1].x();
        min_y = points(msh,cl)[0].y();
        max_y = points(msh,cl)[2].y();
        scaling_x = 1.0/( pow( (max_x - min_x), basis_degree) );
        scaling_y = 1.0/( pow( (max_y - min_y), basis_degree) );
       
    }

    Matrix<VT, Dynamic, 1>
    eval_basis(const point_type& pt)
    {
        Matrix<VT, Dynamic, 1> ret = Matrix<VT, Dynamic, 1>::Zero(basis_size);

        auto bx0 = pt.x() - min_x ;
        auto bx1 = max_x - pt.x() ;
        auto by0 = pt.y() - min_y ;
        auto by1 = max_y - pt.y() ;
        
#ifdef POWER_CACHE
        if ( power_cache.size() != (basis_degree+1)*2 )
          power_cache.resize( (basis_degree+1)*2);
        
        if ( binomial_coeff.size() != (basis_degree+1) )
        binomial_coeff.resize( basis_degree+1 );

        // Creation of the exponential terms
        for (size_t i = 0 ; i <= basis_degree ; i++)
        {
            size_t j = basis_degree - i;
            //if(i == 0)
             //   binomial_coeff[0] = 1.0 ;
            //else if(i == basis_degree)
            //    binomial_coeff[basis_degree+1] = 1.0 ;
           // else
            
            binomial_coeff[i] = binomial_coeff_fx(basis_degree,i);
           
            
            power_cache[2*i]    = binomial_coeff[i]*iexp_pow(bx0, i)*iexp_pow(bx1, j);
            power_cache[2*i+1]  = binomial_coeff[i]*iexp_pow(by0, i)*iexp_pow(by1, j);
          
        }
        
#endif
        size_t pos = 0;
    
        // Case degree FEM = 0
        if (basis_degree == 0){
            VT one = 1.0;
            ret(pos++) = one;
        }
        // Case degree FEM = 1 -> i.e. just vertices
        // Cell ordering:
        // 3 -- 2
        // 0 -- 1
        
        else if(basis_degree == 1){
            for (int pow_y = 0; pow_y <= basis_degree; pow_y +=basis_degree){
                if( pow_y == 0 ) // ordering from sx to dx 0->1
                {
                for (int pow_x = 0; pow_x <= basis_degree; pow_x+=basis_degree){
#ifdef POWER_CACHE
                    auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
#else
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                                            
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                        
                    auto bv = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) * scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
#endif
                    ret(pos++) = bv;
                }
                }
                else{ // ordering from dx to sx 2->3
                    for (int pow_x = basis_degree; pow_x >= 0; pow_x-=basis_degree){
#ifdef POWER_CACHE
                    auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
#else
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                                                                
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                                            
                    auto bv = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) * scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
#endif
                    ret(pos++) = bv;
                }
                }
            }
        }
        
        else
        {
            // case degree FEM >= 2
            int N = basis_degree ;
            int starter = 0;
            int internal_bases = basis_degree + 1 ; // counting the number of layout of internal nodes: initially there are basis_degree +1 nodes
           
            while(internal_bases > 1) // if there is a deeper layout of nodes
            {
                // Vertices nodes of a quadrangular cell. For each loop of internal_bases I add the vertice of "inner quadrangualar nodes" ( it's not properly true speaking of nodes, being a modal approach; better say the picks of the shape functions).
                for (int pow_y = starter; pow_y <= N; pow_y +=(N-starter)){
                  
                    if(pow_y == starter ) // if bottom line -> from sx to dx
                    {
                       // if not external layout of nodes ( vertices level ) -> I don't take extremes since internal layout has 2 nodes less.
                        for (int pow_x = starter; pow_x <= N; pow_x+=(N-starter)){
#ifdef POWER_CACHE
                        auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
#else
                        size_t j_x = basis_degree - pow_x;
                        size_t j_y = basis_degree - pow_y;
                                                                
                        auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                        auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                                            
                        auto bv = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) * scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
#endif
                        ret(pos++) = bv;
                        // std::cout<<"pos vertices bottom is "<<pos<<std::endl;
                        }
                    }
                    else{ // if top line -> from dx to sx
                        for (int pow_x = N; pow_x >= starter; pow_x-=(N-starter)){
#ifdef POWER_CACHE
                        auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
#else
                        size_t j_x = basis_degree - pow_x;
                        size_t j_y = basis_degree - pow_y;
                                                                                    
                        auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                        auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                                                                
                        auto bv = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) * scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
#endif
                        ret(pos++) = bv;
                       // std::cout<<"pos vertices upper is "<<pos<<std::endl;
                        }
                        
                    }
                    
                }
                    
                    
                // After the vertices elements, step by faces bases.

                //size_t pos = 4;
                // face 0 (BOTTOM)
                for (size_t pow_x = starter + 1 ; pow_x <= N - 1 ; pow_x++)
                {
                    size_t pow_y = starter;
#ifdef POWER_CACHE
                    auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
#else
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                            
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                           
                    auto bv = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) * scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
#endif
                    ret(pos++) = bv;
                   //  std::cout<<"pos face 0 is "<<pos<<std::endl;
                }
       
                // face 1 (RIGHT)
                for (size_t pow_y = starter + 1; pow_y <= N - 1 ; pow_y++)
                {
                    size_t pow_x = N;
        #ifdef POWER_CACHE
                    auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
        #else
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                                    
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                   
                    auto bv = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) * scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
        #endif
                    ret(pos++) = bv;
                  //   std::cout<<"pos face 1 is "<<pos<<std::endl;
                }
        
                // face 2 (TOP)
                for (size_t pow_x = N - 1; pow_x >= starter + 1 ; pow_x--)
                {
                    size_t pow_y = N;
        #ifdef POWER_CACHE
                    auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
        #else
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                                    
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                   
                    auto bv = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) * scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
        #endif
                    ret(pos++) = bv;
                  //   std::cout<<"pos face 2 is "<<pos<<std::endl;
                }
        
              
                // face 3 (LEFT)
                
                for (size_t pow_y = N - 1; pow_y >= starter + 1 ; pow_y--)
                {
                    size_t pow_x = starter;
        #ifdef POWER_CACHE
                    auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
        #else
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                                    
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                   
                    auto bv = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) * scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
        #endif
                    ret(pos++) = bv;
                  //   std::cout<<"pos face 3 is "<<pos<<std::endl;
                }
        
                N--; // for each layout of "nodes", I have un unknown less from right
                starter++; // for each layout of "nodes", I have un unknown less from left
                internal_bases -= 2 ; // each layout I have 2 nodes less than the previous
                
            }
            // for B_k, with k even, odd number of bases: there is a central one.
            if( basis_degree % 2 == 0 )
            {
             //   std::cout<<"N is "<<N<< " and starter is "<<starter<<std::endl;
                assert( N == starter );
                size_t pow_x = starter ;
                size_t pow_y = starter ;
#ifdef POWER_CACHE
                auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
#else
                size_t j_x = basis_degree - pow_x;
                size_t j_y = basis_degree - pow_y;
                                
                auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                            
                auto bv = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) * scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
#endif
                ret(pos++) = bv;
               //  std::cout<<"pos %2 is "<<pos<<std::endl;
            }
        }
       // std::cout<<"pos is"<<pos<<" and basis size is "<< basis_size<<std::endl;
        assert(pos == basis_size);

        return ret;
    }

    // Same approach of eval
    Matrix<VT, Dynamic, 2>
    eval_gradients(const point_type& pt)
    {
        Matrix<VT, Dynamic, 2> ret = Matrix<VT, Dynamic, 2>::Zero(basis_size, 2);

        auto bx0 = pt.x() - min_x ;
        auto bx1 = max_x - pt.x() ;
        auto by0 = pt.y() - min_y ;
        auto by1 = max_y - pt.y() ;
        //The stuff below should be in the constructor, but gradient is used less times than eval.
        auto N = basis_degree - 1 ;
        auto coeff_dx = basis_degree/(max_x - min_x); // n/(b-a)
        auto coeff_dy = basis_degree/(max_y - min_y);
        auto scaling_x_bis = 1.0/(pow(max_x - min_x,N) ); // 1/(b-a)^(n-1)
        auto scaling_y_bis = 1.0/(pow(max_y - min_y,N) );
        //std::cout<<"scaling_x_bis "<<scaling_x_bis<<std::endl;
        
#ifdef POWER_CACHE
        if ( power_cache.size() != (basis_degree+1)*2 )
          power_cache.resize( (basis_degree+1)*2) ;
        
        if ( power_cache_bis.size() != (N+1)*2 )
            power_cache_bis.resize( (N+1)*2) ;
        
        if ( binomial_coeff.size() != (basis_degree+1) )
            binomial_coeff.resize( basis_degree+1 ) ;
            
        if ( binomial_coeff_bis.size() != (N+1) )
            binomial_coeff_bis.resize( N+1 ) ;
        
       // Construction of the exponenatial term for bernstein basis B^N and B^(N-1) (useful for derivative)
        for (size_t i = 0; i <= basis_degree ; i++)
        {
            size_t j = basis_degree - i;
            binomial_coeff[i] = binomial_coeff_fx( basis_degree , i );
            power_cache[2*i]    = scaling_x *binomial_coeff[i]*iexp_pow(bx0, i)*iexp_pow(bx1, j);
            power_cache[2*i+1]  = scaling_y *binomial_coeff[i]*iexp_pow(by0, i)*iexp_pow(by1, j);
             
            if( i < basis_degree )
            {
                size_t j_bis = N - i;
                binomial_coeff_bis[i] = binomial_coeff_fx( N , i );
                power_cache_bis[2*i] = scaling_x_bis * binomial_coeff_bis[i]*iexp_pow(bx0, i)*iexp_pow(bx1, j_bis);
                power_cache_bis[2*i+1] = scaling_y_bis * binomial_coeff_bis[i]*iexp_pow(by0, i)*iexp_pow(by1, j_bis);
               
            }
            
        }
#endif

        size_t pos = 0;
        // Case FEM degree = 0
        if (basis_degree == 0){
            VT zero = 0.0;
            ret(pos,0) =  zero;
            ret(pos,1) = zero;
            pos++;
        }
        // Case FEM degree = 1
        else if(basis_degree == 1){
            VT coeffx = 1.0/(max_x - min_x);
            VT coeffy = 1.0/(max_y - min_y);
            //std::cout<<"min x "<<min_x<<"max x "<<max_x<<"min y "<<min_y<<"max y "<<max_y<<std::endl;
            auto px = 0.0 , py = 0.0 ;
            int j = 1; // useful for the sign of the derivative in y
            // Pick just vertices
            for (int pow_y = 0; pow_y <= basis_degree; pow_y +=basis_degree){
                if( pow_y == 0 ) // bottom
                {
                    int i = 1;  // useful for the sign of the derivative in x
                    for (int pow_x = 0; pow_x <= basis_degree; pow_x+=basis_degree){
        #ifdef POWER_CACHE
                    auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
        #else
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                                                        
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                                    
                    px = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) ;
                    py = scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
        #endif
                    
                   
                    ret(pos,0) =  py * coeffx * pow(-1,i);
                    ret(pos,1) = px * coeffy * pow(-1,j);
                    //std::cout<<"i is "<<i<< " and pow (-1,i) = "<<pow(-1,i)<<std::endl;
                    pos++;
                    i++;
                    }
                  
                }
                else{
                    int i = 2;
                    for (int pow_x = basis_degree; pow_x >= 0; pow_x-=basis_degree){
        #ifdef POWER_CACHE
                    auto bv = scaling_x * power_cache[2*pow_x] * scaling_y * power_cache[2*pow_y+1];
        #else
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                                                                            
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                    px = scaling_x * coeff_n_x * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x) ;
                    py = scaling_y * coeff_n_y * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y) ;
        #endif
                   
                    ret(pos,0) =  py * coeffx * pow(-1,i);
                    ret(pos,1) = px * coeffy * pow(-1,j);
                  //  std::cout<<"i is "<<i<< " and pow (-1,i) = "<<pow(-1,i)<<std::endl;
                    pos++;
                    i--;
                    }
                }
                j++;
            }
           
        }
        
        else{
            // if degree FEM >= 2
            //std::cout<<"Bernstein Gradient, basis >=2 "<<std::endl;
            int N_partial = basis_degree ;
            int starter = 0;
            int internal_bases = basis_degree + 1 ;
            
            while(internal_bases > 1) // for each layour of internal node
            {
                for (int pow_y = starter; pow_y <= N_partial ; pow_y +=(N_partial-starter)){
                    
                    if(pow_y == starter ) // Bottom side
                    {
                        
                        for (int pow_x = starter; pow_x <= N_partial; pow_x+=(N_partial-starter)){
                            
                            VT ix_bis = pow_x-1; // element i-1 for derivative
                            VT iy_bis = pow_y-1; // element i-1 for derivative
                            
#ifdef POWER_CACHE
                            std::cout<<"Bernstein Gradient, basis >=2: POWER_CACHE"<<std::endl;
                            auto px = power_cache[2*pow_x];
                            auto py = power_cache[2*pow_y+1];
                            //std::cout<<"px POWER_CACHE "<<px<<std::endl;
                            if ( pow_x == 0 )
                                auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                            else if ( pow_x == basis_degree )
                                auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            else
                                auto dx = coeff_dx*( power_cache_bis[2*ix_bis] - power_cache_bis[2*pow_x] );
                            
                            if ( pow_y == 0 )
                                auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                            else if ( pow_y == basis_degree )
                                auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                            else
                                auto dy = coeff_dy*( power_cache_bis[2*iy_bis+1] - power_cache_bis[2*pow_y+1] );
                            //std::cout<<"dx POWER_CACHE "<<dx<<std::endl;
                            
#else
                            
                            //std::cout<<"Bernstein Gradient, basis >=2: NOT CACHE"<<std::endl;
                            size_t j_x = basis_degree - pow_x;
                            size_t j_y = basis_degree - pow_y;
                            //std::cout<<"j_x PRE "<<j_x<<std::endl;
                            auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                            auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                            //std::cout<<"coeff_n_x PRE "<<coeff_n_x<<std::endl;
                            //std::cout<<"scaling_x PRE "<<scaling_x<<std::endl;
                            //std::cout<<"pow_x PRE "<<pow_x<<std::endl;
                            //std::cout<<"bx0 PRE "<<bx0<<std::endl;
                            //std::cout<<"bx1 PRE "<<bx1<<std::endl;
                            
                            auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                            auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                            //std::cout<<"px PRE "<<px<<std::endl;
                            
                            // DERIVATIVES
                            //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                            auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 ;
                            VT dx = 0.0 , dy = 0.0 ;
                            
                            if ( pow_x == 0 ){
                                //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                dx = -coeff_dx * px_bis0 ;
                                //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                                
                            }
                            else if ( pow_x == basis_degree ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis1 = N - ix_bis;
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_dx * px_bis1 ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else{
                                size_t j_x_bis0 = N - pow_x;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                
                                size_t j_x_bis1 = N - ix_bis;
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_dx * (px_bis1 -px_bis0);
                                //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                            }
                            
                            if ( pow_y == 0 ){
                                //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                               // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                                dy = -coeff_dy *  py_bis0 ;
                                //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                                
                            }
                            else if ( pow_y == basis_degree ){
                                //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                                size_t j_y_bis1 = N - iy_bis;
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_dy * py_bis1  ;
                                //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                            }
                            else{
                                size_t j_y_bis0 = N - pow_y;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                
                                size_t j_y_bis1 = N - iy_bis;
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_dy * ( py_bis1 - py_bis0 );
                            //auto dy=coeff_dy*(power_cache_bis[2*iy_bis+1]-power_cache_bis[2*pow_y+1]);
                            }
                            
                            
#endif
                            //std::cout<<"pos "<<pos<<std::endl;
                            //std::cout<<"dx "<<dx<<std::endl;
                            //std::cout<<"dy "<<dy<<std::endl;
                            //std::cout<<"px "<<px<<std::endl;
                            //std::cout<<"py "<<py<<std::endl;
                            ret(pos,0) = dx*py;
                            ret(pos,1) = px*dy;
                            pos++;
                            
                            
                        }
                    }
                    else{ // Top side
                        for (int pow_x = N_partial; pow_x >= starter; pow_x-=(N_partial-starter)){
                            VT ix_bis = pow_x-1;
                            VT iy_bis = pow_y-1;
                            
#ifdef POWER_CACHE
                            auto px = power_cache[2*pow_x];
                            auto py = power_cache[2*pow_y+1];
                            if ( pow_x == 0 )
                                auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                            else if ( pow_x == basis_degree )
                                auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            else
                                auto dx = coeff_dx*( power_cache_bis[2*ix_bis] - power_cache_bis[2*pow_x] );
                            
                            if ( pow_y == 0 )
                                auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                            else if ( pow_y == basis_degree )
                                auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                            else
                                auto dy = coeff_dy*( power_cache_bis[2*iy_bis+1] - power_cache_bis[2*pow_y+1] );
#else
                            
                            size_t j_x = basis_degree - pow_x;
                            size_t j_y = basis_degree - pow_y;
                            
                            auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                            auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                            
                            auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                            auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                            
                            // DERIVATIVES
                            //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                            auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 ;
                            VT dx = 0.0 , dy = 0.0 ;
                            
                            if ( pow_x == 0 ){
                                //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                dx = -coeff_dx * px_bis0 ;
                                //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                                
                            }
                            else if ( pow_x == basis_degree ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis1 = N - ix_bis;
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_dx * px_bis1 ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else{
                                size_t j_x_bis0 = N - pow_x;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                
                                size_t j_x_bis1 = N - ix_bis;
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_dx * (px_bis1 -px_bis0);
                                //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                            }
                            
                            if ( pow_y == 0 ){
                                //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                               // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                                dy = -coeff_dy *  py_bis0 ;
                                //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                                
                            }
                            else if ( pow_y == basis_degree ){
                                //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                                size_t j_y_bis1 = N - iy_bis;
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_dy * py_bis1  ;
                                //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                            }
                            else{
                                size_t j_y_bis0 = N - pow_y;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                
                                size_t j_y_bis1 = N - iy_bis;
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_dy * ( py_bis1 - py_bis0 );
                            //auto dy=coeff_dy*(power_cache_bis[2*iy_bis+1]-power_cache_bis[2*pow_y+1]);
                            }
                            
                            
#endif
                            
                            ret(pos,0) = dx*py;
                            ret(pos,1) = px*dy;
                            pos++;
                            
                            
                        } // loop for top side vertices
                        
                        
                        
                        
                    } // else bottom-top side
                    
                } // end loop for vertices
                
                
                // face 0 (BOTTOM)
                for (size_t pow_x = starter + 1 ; pow_x <= N_partial - 1 ; pow_x++)
                {
                    size_t pow_y = starter;
                    
                    VT ix_bis = pow_x-1;
                    VT iy_bis = pow_y-1;
                    
#ifdef POWER_CACHE
                    auto px = power_cache[2*pow_x];
                    auto py = power_cache[2*pow_y+1];
                    if ( pow_x == 0 )
                        auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                    else if ( pow_x == basis_degree )
                        auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                    else
                        auto dx = coeff_dx*( power_cache_bis[2*ix_bis] - power_cache_bis[2*pow_x] );
                    
                    if ( pow_y == 0 )
                        auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                    else if ( pow_y == basis_degree )
                        auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                    else
                        auto dy = coeff_dy*( power_cache_bis[2*iy_bis+1] - power_cache_bis[2*pow_y+1] );
#else
                    
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                    
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                    
                    auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                    auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                    //std::cout<<"j_x Face0 "<<j_x<<std::endl;
                    //std::cout<<"coeff_n_x Face0 "<<coeff_n_x<<std::endl;
                    //std::cout<<"scaling_x Face0 "<<scaling_x<<std::endl;
                    //std::cout<<"pow_x Face0 "<<pow_x<<std::endl;
                    //std::cout<<"bx0 Face0 "<<bx0<<std::endl;
                    //std::cout<<"bx1 Face0 "<<bx1<<std::endl;
                    
                    
                    // DERIVATIVES
                    //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                    auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 ;
                    VT dx = 0.0 , dy = 0.0 ;
                    
                    if ( pow_x == 0 ){
                        //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                        size_t j_x_bis0 = N - pow_x;
                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                        dx = -coeff_dx * px_bis0 ;
                        //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                        
                    }
                    else if ( pow_x == basis_degree ){
                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                        size_t j_x_bis1 = N - ix_bis;
                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                        dx = coeff_dx * px_bis1 ;
                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                    }
                    else{
                        size_t j_x_bis0 = N - pow_x;
                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                        
                        size_t j_x_bis1 = N - ix_bis;
                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                        dx = coeff_dx * (px_bis1 -px_bis0);
                        //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                    }
                    
                    if ( pow_y == 0 ){
                        //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                        size_t j_y_bis0 = N - pow_y;
                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                       // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                        dy = -coeff_dy *  py_bis0 ;
                        //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                        
                    }
                    else if ( pow_y == basis_degree ){
                        //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                        size_t j_y_bis1 = N - iy_bis;
                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                        dy = coeff_dy * py_bis1  ;
                        //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                    }
                    else{
                        size_t j_y_bis0 = N - pow_y;
                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                        
                        size_t j_y_bis1 = N - iy_bis;
                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                        dy = coeff_dy * ( py_bis1 - py_bis0 );
                    //auto dy=coeff_dy*(power_cache_bis[2*iy_bis+1]-power_cache_bis[2*pow_y+1]);
                    }
                    
#endif
                    
                    ret(pos,0) = dx*py;
                    ret(pos,1) = px*dy;
                    pos++;
                    
                }
                
                // face 1 (RIGHT)
                for (size_t pow_y = starter + 1; pow_y <= N_partial - 1 ; pow_y++)
                {
                    size_t pow_x = N_partial;
                    VT ix_bis = pow_x-1;
                    VT iy_bis = pow_y-1;
                    
#ifdef POWER_CACHE
                    auto px = power_cache[2*pow_x];
                    auto py = power_cache[2*pow_y+1];
                    if ( pow_x == 0 )
                        auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                    else if ( pow_x == basis_degree )
                        auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                    else
                        auto dx = coeff_dx*( power_cache_bis[2*ix_bis] - power_cache_bis[2*pow_x] );
                    
                    if ( pow_y == 0 )
                        auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                    else if ( pow_y == basis_degree )
                        auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                    else
                        auto dy = coeff_dy*( power_cache_bis[2*iy_bis+1] - power_cache_bis[2*pow_y+1] );
#else
                    
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                    
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                    
                    auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                    auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                    
                    // DERIVATIVES
                    //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                    auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 ;
                    VT dx = 0.0 , dy = 0.0 ;
                    
                    if ( pow_x == 0 ){
                        //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                        size_t j_x_bis0 = N - pow_x;
                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                        dx = -coeff_dx * px_bis0 ;
                        //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                        
                    }
                    else if ( pow_x == basis_degree ){
                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                        size_t j_x_bis1 = N - ix_bis;
                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                        dx = coeff_dx * px_bis1 ;
                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                    }
                    else{
                        size_t j_x_bis0 = N - pow_x;
                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                        
                        size_t j_x_bis1 = N - ix_bis;
                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                        dx = coeff_dx * (px_bis1 -px_bis0);
                        //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                    }
                    
                    if ( pow_y == 0 ){
                        //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                        size_t j_y_bis0 = N - pow_y;
                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                       // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                        dy = -coeff_dy *  py_bis0 ;
                        //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                        
                    }
                    else if ( pow_y == basis_degree ){
                        //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                        size_t j_y_bis1 = N - iy_bis;
                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                        dy = coeff_dy * py_bis1  ;
                        //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                    }
                    else{
                        size_t j_y_bis0 = N - pow_y;
                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                        
                        size_t j_y_bis1 = N - iy_bis;
                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                        dy = coeff_dy * ( py_bis1 - py_bis0 );
                    //auto dy=coeff_dy*(power_cache_bis[2*iy_bis+1]-power_cache_bis[2*pow_y+1]);
                    }
                    
                    
#endif
                    
                    ret(pos,0) = dx*py;
                    ret(pos,1) = px*dy;
                    pos++;
                }
                
                // face 2 (TOP)
                for (size_t pow_x = N_partial - 1; pow_x >= starter + 1 ; pow_x--)
                {
                    size_t pow_y = N_partial;
                    VT ix_bis = pow_x-1;
                    VT iy_bis = pow_y-1;
                    
#ifdef POWER_CACHE
                    auto px = power_cache[2*pow_x];
                    auto py = power_cache[2*pow_y+1];
                    if ( pow_x == 0 )
                        auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                    else if ( pow_x == basis_degree )
                        auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                    else
                        auto dx = coeff_dx*( power_cache_bis[2*ix_bis] - power_cache_bis[2*pow_x] );
                    
                    if ( pow_y == 0 )
                        auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                    else if ( pow_y == basis_degree )
                        auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                    else
                        auto dy = coeff_dy*( power_cache_bis[2*iy_bis+1] - power_cache_bis[2*pow_y+1] );
#else
                    
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                    
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                    
                    auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                    auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                    
                    // DERIVATIVES
                    //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                    auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 ;
                    VT dx = 0.0 , dy = 0.0 ;
                    
                    if ( pow_x == 0 ){
                        //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                        size_t j_x_bis0 = N - pow_x;
                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                        dx = -coeff_dx * px_bis0 ;
                        //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                        
                    }
                    else if ( pow_x == basis_degree ){
                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                        size_t j_x_bis1 = N - ix_bis;
                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                        dx = coeff_dx * px_bis1 ;
                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                    }
                    else{
                        size_t j_x_bis0 = N - pow_x;
                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                        
                        size_t j_x_bis1 = N - ix_bis;
                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                        dx = coeff_dx * (px_bis1 -px_bis0);
                        //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                    }
                    
                    if ( pow_y == 0 ){
                        //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                        size_t j_y_bis0 = N - pow_y;
                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                       // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                        dy = -coeff_dy *  py_bis0 ;
                        //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                        
                    }
                    else if ( pow_y == basis_degree ){
                        //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                        size_t j_y_bis1 = N - iy_bis;
                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                        dy = coeff_dy * py_bis1  ;
                        //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                    }
                    else{
                        size_t j_y_bis0 = N - pow_y;
                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                        
                        size_t j_y_bis1 = N - iy_bis;
                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                        dy = coeff_dy * ( py_bis1 - py_bis0 );
                    //auto dy=coeff_dy*(power_cache_bis[2*iy_bis+1]-power_cache_bis[2*pow_y+1]);
                    }
                    
                    
#endif
                    
                    ret(pos,0) = dx*py;
                    ret(pos,1) = px*dy;
                    pos++;
                }
                
                
                // face 3 (LEFT)
                
                for (size_t pow_y = N_partial - 1; pow_y >= starter + 1 ; pow_y--)
                {
                    size_t pow_x = starter;
                    VT ix_bis = pow_x-1;
                    VT iy_bis = pow_y-1;
                    
#ifdef POWER_CACHE
                    auto px = power_cache[2*pow_x];
                    auto py = power_cache[2*pow_y+1];
                    if ( pow_x == 0 )
                        auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                    else if ( pow_x == basis_degree )
                        auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                    else
                        auto dx = coeff_dx*( power_cache_bis[2*ix_bis] - power_cache_bis[2*pow_x] );
                    
                    if ( pow_y == 0 )
                        auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                    else if ( pow_y == basis_degree )
                        auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                    else
                        auto dy = coeff_dy*( power_cache_bis[2*iy_bis+1] - power_cache_bis[2*pow_y+1] );
#else
                    
                    size_t j_x = basis_degree - pow_x;
                    size_t j_y = basis_degree - pow_y;
                    
                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                    
                    auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                    auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                    
                    // DERIVATIVES
                    //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                    auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 ;
                    VT dx = 0.0 , dy = 0.0 ;
                    
                    if ( pow_x == 0 ){
                        //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                        size_t j_x_bis0 = N - pow_x;
                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                        dx = -coeff_dx * px_bis0 ;
                        //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                        
                    }
                    else if ( pow_x == basis_degree ){
                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                        size_t j_x_bis1 = N - ix_bis;
                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                        dx = coeff_dx * px_bis1 ;
                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                    }
                    else{
                        size_t j_x_bis0 = N - pow_x;
                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                        
                        size_t j_x_bis1 = N - ix_bis;
                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                        dx = coeff_dx * (px_bis1 -px_bis0);
                        //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                    }
                    
                    if ( pow_y == 0 ){
                        //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                        size_t j_y_bis0 = N - pow_y;
                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                       // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                        dy = -coeff_dy *  py_bis0 ;
                        //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                        
                    }
                    else if ( pow_y == basis_degree ){
                        //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                        size_t j_y_bis1 = N - iy_bis;
                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                        dy = coeff_dy * py_bis1  ;
                        //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                    }
                    else{
                        size_t j_y_bis0 = N - pow_y;
                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                        
                        size_t j_y_bis1 = N - iy_bis;
                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                        dy = coeff_dy * ( py_bis1 - py_bis0 );
                    //auto dy=coeff_dy*(power_cache_bis[2*iy_bis+1]-power_cache_bis[2*pow_y+1]);
                    }
                    
                    
#endif
                    
                    ret(pos,0) = dx*py;
                    ret(pos,1) = px*dy;
                    pos++;
                }
                
                N_partial--;
                starter++;
                internal_bases -= 2;
                
            }
            
            // for B_k, with k even, I.E. odd number of bases: there is a central one.
            if( basis_degree % 2 == 0 )
            {
                //  std::cout<<"N is "<<N_partial<< " and starter is "<<starter<<std::endl;
                assert( N_partial == starter );
                size_t pow_x = starter ;
                size_t pow_y = starter ;
                VT ix_bis = pow_x-1;
                VT iy_bis = pow_y-1;
                
#ifdef POWER_CACHE
                auto px = power_cache[2*pow_x];
                auto py = power_cache[2*pow_y+1];
                if ( pow_x == 0 )
                    auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                else if ( pow_x == basis_degree )
                    auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                else
                    auto dx = coeff_dx*( power_cache_bis[2*ix_bis] - power_cache_bis[2*pow_x] );
                
                if ( pow_y == 0 )
                    auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                else if ( pow_y == basis_degree )
                    auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                else
                    auto dy = coeff_dy*( power_cache_bis[2*iy_bis+1] - power_cache_bis[2*pow_y+1] );
#else
                
                size_t j_x = basis_degree - pow_x;
                size_t j_y = basis_degree - pow_y;
                
                auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                
                auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                
                // DERIVATIVES
                //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 ;
                VT dx = 0.0 , dy = 0.0 ;
                
                if ( pow_x == 0 ){
                    //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                    size_t j_x_bis0 = N - pow_x;
                    auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                    auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                    dx = -coeff_dx * px_bis0 ;
                    //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                    
                }
                else if ( pow_x == basis_degree ){
                    //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                    size_t j_x_bis1 = N - ix_bis;
                    auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                    auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                    dx = coeff_dx * px_bis1 ;
                    //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                }
                else{
                    size_t j_x_bis0 = N - pow_x;
                    auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                    auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                    
                    size_t j_x_bis1 = N - ix_bis;
                    auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                    auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                    dx = coeff_dx * (px_bis1 -px_bis0);
                    //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                }
                
                if ( pow_y == 0 ){
                    //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                    size_t j_y_bis0 = N - pow_y;
                    auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                    auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                   // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                    dy = -coeff_dy *  py_bis0 ;
                    //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                    
                }
                else if ( pow_y == basis_degree ){
                    //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                    size_t j_y_bis1 = N - iy_bis;
                    auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                    auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                    dy = coeff_dy * py_bis1  ;
                    //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                }
                else{
                    size_t j_y_bis0 = N - pow_y;
                    auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                    auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                    
                    size_t j_y_bis1 = N - iy_bis;
                    auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                    auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                    dy = coeff_dy * ( py_bis1 - py_bis0 );
                //auto dy=coeff_dy*(power_cache_bis[2*iy_bis+1]-power_cache_bis[2*pow_y+1]);
                }
                
                
#endif
                
                ret(pos,0) = dx*py;
                ret(pos,1) = px*dy;
                pos++;
            }
        }

          //  std::cout<<"GRADIENTE: pos is"<<pos<<" and basis size is "<< basis_size<<std::endl;
            assert(pos == basis_size);

            return ret;
        

    }

    
    // Same approach of eval
    Matrix<VT, Dynamic, 1>
    eval_divergence(const point_type& pt)
    {
                Matrix<VT, Dynamic, 1> ret = Matrix<VT, Dynamic, 1>::Zero(basis_size, 1);

                auto bx0 = pt.x() - min_x ;
                auto bx1 = max_x - pt.x() ;
                auto by0 = pt.y() - min_y ;
                auto by1 = max_y - pt.y() ;
                //The stuff below should be in the constructor, but gradient is used less times than eval.
                auto N = basis_degree - 2 ;
                auto coeff_d2x = basis_degree*(basis_degree-1)/((max_x - min_x)*(max_x - min_x));
                auto coeff_d2y = basis_degree*(basis_degree-1)/((max_y - min_y)*(max_y - min_y));
                auto scaling_x_bis = 1.0/(pow(max_x - min_x,N) ); // 1/(b-a)^(n-1)
                auto scaling_y_bis = 1.0/(pow(max_y - min_y,N) );
                //std::cout<<"scaling_x_bis "<<scaling_x_bis<<std::endl;
                
        #ifdef POWER_CACHE
                if ( power_cache.size() != (basis_degree+1)*2 )
                  power_cache.resize( (basis_degree+1)*2) ;
                
                if ( power_cache_bis.size() != (N+1)*2 )
                    power_cache_bis.resize( (N+1)*2) ;
                
                if ( binomial_coeff.size() != (basis_degree+1) )
                    binomial_coeff.resize( basis_degree+1 ) ;
                    
                if ( binomial_coeff_bis.size() != (N+1) )
                    binomial_coeff_bis.resize( N+1 ) ;
                
               // Construction of the exponenatial term for bernstein basis B^N and B^(N-2) (useful for derivative)
                for (size_t i = 0; i <= basis_degree ; i++)
                {
                    size_t j = basis_degree - i;
                    binomial_coeff[i] = binomial_coeff_fx( basis_degree , i );
                    power_cache[2*i]    = scaling_x *binomial_coeff[i]*iexp_pow(bx0, i)*iexp_pow(bx1, j);
                    power_cache[2*i+1]  = scaling_y *binomial_coeff[i]*iexp_pow(by0, i)*iexp_pow(by1, j);
                     
                    if( i < basis_degree - 1 )
                    {
                        size_t j_bis = N - i;
                        binomial_coeff_bis[i] = binomial_coeff_fx( N , i );
                        power_cache_bis[2*i] = scaling_x_bis * binomial_coeff_bis[i]*iexp_pow(bx0, i)*iexp_pow(bx1, j_bis);
                        power_cache_bis[2*i+1] = scaling_y_bis * binomial_coeff_bis[i]*iexp_pow(by0, i)*iexp_pow(by1, j_bis);
                       
                    }
                    
                }
        #endif

                size_t pos = 0;
                // Case FEM degree = 0
                if (basis_degree == 0 || basis_degree == 1){
                    VT zero = 0.0;
                    while(pos < basis_size){
                        ret(pos) =  zero;
                        pos++;
                    }
                }
                else{
                    // if degree FEM >= 2
                    //std::cout<<"Bernstein Gradient, basis >=2 "<<std::endl;
                    int N_partial = basis_degree ;
                    int starter = 0;
                    int internal_bases = basis_degree + 1 ;
                    
                    while(internal_bases > 1) // for each layout of internal node
                    {
                        for (int pow_y = starter; pow_y <= N_partial ; pow_y +=(N_partial-starter)){
                            
                            if(pow_y == starter ) // Bottom side
                            {
                                
                                for (int pow_x = starter; pow_x <= N_partial; pow_x+=(N_partial-starter))
                                {
                                    
                                    VT ix_bis = pow_x-1; // element i-1 for derivative
                                    VT iy_bis = pow_y-1; // element i-1 for derivative
                                    VT ix2_bis = pow_x-2; // element i-2 for derivative
                                    VT iy2_bis = pow_y-2; // element i-2 for derivative
                                    
        #ifdef POWER_CACHE
                                    std::cout<<"Bernstein Gradient, basis >=2: POWER_CACHE"<<std::endl;
                                    auto px = power_cache[2*pow_x];
                                    auto py = power_cache[2*pow_y+1];
                                    //std::cout<<"px POWER_CACHE "<<px<<std::endl;
                                    if ( pow_x == 0 )
                                        auto dx = coeff_d2x * power_cache_bis[2*pow_x] ;
                                    else if ( pow_x == basis_degree )
                                        auto dx = coeff_d2x * power_cache_bis[2*ix2_bis] ;
                                    else if ( pow_x == 1 )
                                        auto dx = coeff_d2x * ( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] );
                                    else if ( pow_x == basis_degree - 1 )
                                        auto dx = coeff_d2x * ( power_cache_bis[2*ix2_bis] -2*power_cache_bis[2*ix_bis] );
                                    else
                                        auto dx = coeff_d2x*( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] + power_cache_bis[2*ix2_bis] );
                                    
                                    if ( pow_y == 0 )
                                        auto dy = coeff_d2y*power_cache_bis[2*pow_y+1] ;
                                    else if ( pow_y == basis_degree )
                                        auto dy = coeff_d2y*power_cache_bis[2*iy2_bis+1] ;
                                    else if ( pow_y == 1 )
                                        auto dy = coeff_d2y * ( power_cache_bis[2*pow_y+1] -2*power_cache_bis[2*iy_bis+1] );
                                    else if ( pow_y == basis_degree - 1 )
                                        auto dy = coeff_d2y * ( power_cache_bis[2*iy2_bis+1] -2*power_cache_bis[2*iy_bis+1] );
                                    else
                                        auto dy = coeff_d2y*( power_cache_bis[2*pow_y+1] - 2*power_cache_bis[2*iy_bis+1] + power_cache_bis[2*iy2_bis+1] );
                                    //std::cout<<"dx POWER_CACHE "<<dx<<std::endl;
                                    
        #else
                                    
                                    //std::cout<<"Bernstein Gradient, basis >=2: NOT CACHE"<<std::endl;
                                    size_t j_x = basis_degree - pow_x;
                                    size_t j_y = basis_degree - pow_y;
                                    //std::cout<<"j_x PRE "<<j_x<<std::endl;
                                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                    //std::cout<<"coeff_n_x PRE "<<coeff_n_x<<std::endl;
                                    //std::cout<<"scaling_x PRE "<<scaling_x<<std::endl;
                                    //std::cout<<"pow_x PRE "<<pow_x<<std::endl;
                                    //std::cout<<"bx0 PRE "<<bx0<<std::endl;
                                    //std::cout<<"bx1 PRE "<<bx1<<std::endl;
                                    
                                    auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                                    auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                                    //std::cout<<"px PRE "<<px<<std::endl;
                                    
                                    // DERIVATIVES
                                    //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                                    auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 , px_bis2 = 0.0 , py_bis2 = 0.0 ;
                                    VT dx = 0.0 , dy = 0.0 ;
                                    
                                    if ( pow_x == 0 ){
                                        //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                                        size_t j_x_bis0 = N - pow_x;
                                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                        dx = coeff_d2x * px_bis0 ;
                                        //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                                        
                                    }
                                    else if ( pow_x == basis_degree ){
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_x_bis1 = N - ix2_bis;
                                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                        dx = coeff_d2x * px_bis1 ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                    else if ( pow_x == 1 ){
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_x_bis0 = N - pow_x;
                                        size_t j_x_bis1 = N - ix_bis;
                                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                        auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                        dx = coeff_d2x * (px_bis0 - 2*px_bis1) ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                    
                                    else if ( pow_x == basis_degree - 1 ){
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_x_bis0 = N - ix_bis;
                                        size_t j_x_bis1 = N - ix2_bis;
                                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,ix_bis);
                                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                        auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis0);
                                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                        dx = coeff_d2x * (px_bis1 - 2*px_bis0) ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                    else{
                                        size_t j_x_bis0 = N - pow_x;
                                        size_t j_x_bis1 = N - ix_bis;
                                        size_t j_x_bis2 = N - ix2_bis;
                                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                        auto coeff_n_x_bis2 = binomial_coeff_fx(N,ix2_bis);
                                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                        auto px_bis2 = scaling_x_bis * coeff_n_x_bis2 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis2);
                                        
                                        
                                        dx = coeff_d2x * (px_bis0 -2*px_bis1 +px_bis2);
                                        //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                                    }
                                    
                                    if ( pow_y == 0 ){
                                        //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                                        size_t j_y_bis0 = N - pow_y;
                                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                       // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                                        dy = coeff_d2y *  py_bis0 ;
                                        //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                                        
                                    }
                                    else if ( pow_y == basis_degree ){
                                        //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                                        size_t j_y_bis1 = N - iy2_bis;
                                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                        dy = coeff_d2y * py_bis1  ;
                                        //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                                    }
                                    else if ( pow_y == 1 ){
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_y_bis0 = N - pow_y;
                                        size_t j_y_bis1 = N - iy_bis;
                                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                        auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                        dy = coeff_d2y * (py_bis0 - 2*py_bis1) ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                    else if ( pow_y == basis_degree - 1 ){
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_y_bis0 = N - iy_bis;
                                        size_t j_y_bis1 = N - iy2_bis;
                                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,iy_bis);
                                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                        auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis0);
                                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                        dy = coeff_d2y * (py_bis1 - 2*py_bis0) ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                    else{
                                        size_t j_y_bis0 = N - pow_y;
                                        size_t j_y_bis1 = N - iy_bis;
                                        size_t j_y_bis2 = N - iy2_bis;
                                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                        auto coeff_n_y_bis2 = binomial_coeff_fx(N,iy2_bis);
                                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                        auto py_bis2 = scaling_y_bis * coeff_n_y_bis2 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis2);
                                        
                                    
                                        dy = coeff_d2y * ( py_bis0 - 2*py_bis1 + py_bis2 );
                                    //auto dy=coeff_dy*(power_cache_bis[2*iy_bis+1]-power_cache_bis[2*pow_y+1]);
                                    }
                                    
                                    
        #endif
                                    //std::cout<<"pos "<<pos<<std::endl;
                                    //std::cout<<"dx "<<dx<<std::endl;
                                    //std::cout<<"dy "<<dy<<std::endl;
                                    //std::cout<<"px "<<px<<std::endl;
                                    //std::cout<<"py "<<py<<std::endl;
                                    ret(pos) = dx*py + px*dy;
                                    pos++;
                                    
                                    
                                }
                            }
                            else{ // Top side
                                for (int pow_x = N_partial; pow_x >= starter; pow_x-=(N_partial-starter))
                                {
                                    
                                    VT ix_bis = pow_x-1; // element i-1 for derivative
                                    VT iy_bis = pow_y-1; // element i-1 for derivative
                                    VT ix2_bis = pow_x-2; // element i-2 for derivative
                                    VT iy2_bis = pow_y-2; // element i-2 for derivative
                                                                
            #ifdef POWER_CACHE
                                                            
                                    std::cout<<"Bernstein Gradient, basis >=2: POWER_CACHE"<<std::endl;
                                    auto px = power_cache[2*pow_x];
                                    auto py = power_cache[2*pow_y+1];
                                    //std::cout<<"px POWER_CACHE "<<px<<std::endl;
                                    if ( pow_x == 0 )
                                        auto dx = coeff_d2x * power_cache_bis[2*pow_x] ;
                                    else if ( pow_x == basis_degree )
                                        auto dx = coeff_d2x * power_cache_bis[2*ix2_bis] ;
                                    else if ( pow_x == 1 )
                                        auto dx = coeff_d2x * ( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] );
                                    else if ( pow_x == basis_degree - 1 )
                                        auto dx = coeff_d2x * ( power_cache_bis[2*ix2_bis] -2*power_cache_bis[2*ix_bis] );
                                    else
                                        auto dx = coeff_d2x*( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] + power_cache_bis[2*ix2_bis] );
                                                                
                                        
                                    if ( pow_y == 0 )
                                        auto dy = coeff_d2y*power_cache_bis[2*pow_y+1] ;
                                    else if ( pow_y == basis_degree )
                                        auto dy = coeff_d2y*power_cache_bis[2*iy2_bis+1] ;
                                    else if ( pow_y == 1 )
                                        auto dy = coeff_d2y * ( power_cache_bis[2*pow_y+1] -2*power_cache_bis[2*iy_bis+1] );
                                    else if ( pow_y == basis_degree - 1 )
                                        auto dy = coeff_d2y * ( power_cache_bis[2*iy2_bis+1] -2*power_cache_bis[2*iy_bis+1] );
                                    else
                                        auto dy = coeff_d2y*( power_cache_bis[2*pow_y+1] - 2*power_cache_bis[2*iy_bis+1] + power_cache_bis[2*iy2_bis+1] );
                                    //std::cout<<"dx POWER_CACHE "<<dx<<std::endl;
                                                                
            #else
                                                                
                                    //std::cout<<"Bernstein Gradient, basis >=2: NOT CACHE"<<std::endl;
                                    size_t j_x = basis_degree - pow_x;
                                    size_t j_y = basis_degree - pow_y;
                                    //std::cout<<"j_x PRE "<<j_x<<std::endl;
                                    auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                                    auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                                    //std::cout<<"coeff_n_x PRE "<<coeff_n_x<<std::endl;
                                    //std::cout<<"scaling_x PRE "<<scaling_x<<std::endl;
                                    //std::cout<<"pow_x PRE "<<pow_x<<std::endl;
                                    //std::cout<<"bx0 PRE "<<bx0<<std::endl;
                                    //std::cout<<"bx1 PRE "<<bx1<<std::endl;
                                                                
                                    auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                                    auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                                    //std::cout<<"px PRE "<<px<<std::endl;
                                                                
                                    // DERIVATIVES
                                    //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                                    auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 , px_bis2 = 0.0 , py_bis2 = 0.0 ;
                                    VT dx = 0.0 , dy = 0.0 ;
                                                                
                                    if ( pow_x == 0 ){
                                                                
                                        //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                                        size_t j_x_bis0 = N - pow_x;
                                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                        dx = coeff_d2x * px_bis0 ;
                                        //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                                                                    
                                    }
                                    else if ( pow_x == basis_degree ){
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_x_bis1 = N - ix2_bis;
                                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                        dx = coeff_d2x * px_bis1 ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                    else if ( pow_x == 1 ){
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_x_bis0 = N - pow_x;
                                        size_t j_x_bis1 = N - ix_bis;
                                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                        auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                        dx = coeff_d2x * (px_bis0 - 2*px_bis1) ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                                                
                                    
                                    else if ( pow_x == basis_degree - 1 ){
                                        
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_x_bis0 = N - ix_bis;
                                        size_t j_x_bis1 = N - ix2_bis;
                                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,ix_bis);
                                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                        auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis0);
                                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                        dx = coeff_d2x * (px_bis1 - 2*px_bis0) ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                    else{
                                        size_t j_x_bis0 = N - pow_x;
                                        size_t j_x_bis1 = N - ix_bis;
                                        size_t j_x_bis2 = N - ix2_bis;
                                        auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                        auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                        auto coeff_n_x_bis2 = binomial_coeff_fx(N,ix2_bis);
                                        auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                        auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                        auto px_bis2 = scaling_x_bis * coeff_n_x_bis2 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis2);
                                                                    
                                                                    
                                        dx = coeff_d2x * (px_bis0 -2*px_bis1 +px_bis2);
                                        //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                                    }
                                                            
                                                            
                                    if ( pow_y == 0 ){
                                        //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                                        size_t j_y_bis0 = N - pow_y;
                                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                        // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                                        dy = coeff_d2y *  py_bis0 ;
                                        //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                                                                    
                                    }
                                    else if ( pow_y == basis_degree ){
                                        //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                                        size_t j_y_bis1 = N - iy2_bis;
                                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                        dy = coeff_d2y * py_bis1  ;
                                        //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                                    }
                                    else if ( pow_y == 1 ){
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_y_bis0 = N - pow_y;
                                        size_t j_y_bis1 = N - iy_bis;
                                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                        auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                        dy = coeff_d2y * (py_bis0 - 2*py_bis1) ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                    else if ( pow_y == basis_degree - 1 ){
                                        //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                        size_t j_y_bis0 = N - iy_bis;
                                        size_t j_y_bis1 = N - iy2_bis;
                                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,iy_bis);
                                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                        auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis0);
                                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                        dy = coeff_d2y * (py_bis1 - 2*py_bis0) ;
                                        //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                    }
                                    else{
                                        size_t j_y_bis0 = N - pow_y;
                                        size_t j_y_bis1 = N - iy_bis;
                                        size_t j_y_bis2 = N - iy2_bis;
                                        auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                        auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                        auto coeff_n_y_bis2 = binomial_coeff_fx(N,iy2_bis);
                                        auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                        auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                        auto py_bis2 = scaling_y_bis * coeff_n_y_bis2 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis2);
                                                                    
                                                                
                                        dy = coeff_d2y * ( py_bis0 - 2*py_bis1 + py_bis2 );
                                        
                                    }
                                                                
            #endif
                                    //std::cout<<"pos "<<pos<<std::endl;
                                    //std::cout<<"dx "<<dx<<std::endl;
                                    //std::cout<<"dy "<<dy<<std::endl;
                                    //std::cout<<"px "<<px<<std::endl;
                                    //std::cout<<"py "<<py<<std::endl;
                                    ret(pos) = dx*py + px*dy;
                                    pos++;
                                                                
                                                                
                                    
                                    
                                } // loop for top side vertices
                                
                                
                                
                                
                            } // else bottom-top side
                            
                        } // end loop for vertices
                        
                        
                        // face 0 (BOTTOM)
                        for (size_t pow_x = starter + 1 ; pow_x <= N_partial - 1 ; pow_x++)
                        {
                            size_t pow_y = starter;
                            VT ix_bis = pow_x-1; // element i-1 for derivative
                            VT iy_bis = pow_y-1; // element i-1 for derivative
                            VT ix2_bis = pow_x-2; // element i-2 for derivative
                            VT iy2_bis = pow_y-2; // element i-2 for derivative
                                                                                
    #ifdef POWER_CACHE
                                                                            
                            std::cout<<"Bernstein Gradient, basis >=2: POWER_CACHE"<<std::endl;
                            auto px = power_cache[2*pow_x];
                            auto py = power_cache[2*pow_y+1];
                            //std::cout<<"px POWER_CACHE "<<px<<std::endl;
                            if ( pow_x == 0 )
                                auto dx = coeff_d2x * power_cache_bis[2*pow_x] ;
                            else if ( pow_x == basis_degree )
                                auto dx = coeff_d2x * power_cache_bis[2*ix2_bis] ;
                            else if ( pow_x == 1 )
                                auto dx = coeff_d2x * ( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] );
                            else if ( pow_x == basis_degree - 1 )
                                auto dx = coeff_d2x * ( power_cache_bis[2*ix2_bis] -2*power_cache_bis[2*ix_bis] );
                            else
                                auto dx = coeff_d2x*( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] + power_cache_bis[2*ix2_bis] );
                                                                                
                                                    
                            if ( pow_y == 0 )
                                auto dy = coeff_d2y*power_cache_bis[2*pow_y+1] ;
                            else if ( pow_y == basis_degree )
                                auto dy = coeff_d2y*power_cache_bis[2*iy2_bis+1] ;
                            else if ( pow_y == 1 )
                                auto dy = coeff_d2y * ( power_cache_bis[2*pow_y+1] -2*power_cache_bis[2*iy_bis+1] );
                            else if ( pow_y == basis_degree - 1 )
                                auto dy = coeff_d2y * ( power_cache_bis[2*iy2_bis+1] -2*power_cache_bis[2*iy_bis+1] );
                            else
                                auto dy = coeff_d2y*( power_cache_bis[2*pow_y+1] - 2*power_cache_bis[2*iy_bis+1] + power_cache_bis[2*iy2_bis+1] );
                            //std::cout<<"dx POWER_CACHE "<<dx<<std::endl;
                                                                                
            #else
                                                                            
                            //std::cout<<"Bernstein Gradient, basis >=2: NOT CACHE"<<std::endl;
                            size_t j_x = basis_degree - pow_x;
                            size_t j_y = basis_degree - pow_y;
                            //std::cout<<"j_x PRE "<<j_x<<std::endl;
                            auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                            auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                            //std::cout<<"coeff_n_x PRE "<<coeff_n_x<<std::endl;
                            //std::cout<<"scaling_x PRE "<<scaling_x<<std::endl;
                            //std::cout<<"pow_x PRE "<<pow_x<<std::endl;
                            //std::cout<<"bx0 PRE "<<bx0<<std::endl;
                            //std::cout<<"bx1 PRE "<<bx1<<std::endl;
                                                                                
                            auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                            auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                            //std::cout<<"px PRE "<<px<<std::endl;
                                                                            
                            // DERIVATIVES
                            //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                            auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 , px_bis2 = 0.0 , py_bis2 = 0.0 ;
                            VT dx = 0.0 , dy = 0.0 ;
                                                                                
                            if ( pow_x == 0 ){
                                                                                
                                //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                dx = coeff_d2x * px_bis0 ;
                                //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                                                                                    
                            }
                            else if ( pow_x == basis_degree ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis1 = N - ix2_bis;
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * px_bis1 ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else if ( pow_x == 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                size_t j_x_bis1 = N - ix_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * (px_bis0 - 2*px_bis1) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                                                                                
                            else if ( pow_x == basis_degree - 1 ){
                                                    
                            //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                            size_t j_x_bis0 = N - ix_bis;
                            size_t j_x_bis1 = N - ix2_bis;
                            auto coeff_n_x_bis0 = binomial_coeff_fx(N,ix_bis);
                            auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                            auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis0);
                            auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                            dx = coeff_d2x * (px_bis1 - 2*px_bis0) ;
                            //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            
                            }
                            else{
                                size_t j_x_bis0 = N - pow_x;
                                size_t j_x_bis1 = N - ix_bis;
                                size_t j_x_bis2 = N - ix2_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto coeff_n_x_bis2 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                auto px_bis2 = scaling_x_bis * coeff_n_x_bis2 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis2);
                                                                                
                                                                                    
                                dx = coeff_d2x * (px_bis0 -2*px_bis1 +px_bis2);
                                //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                            }
                                                                            
                                                                            
                            
                            if ( pow_y == 0 ){
                                //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                                dy = coeff_d2y *  py_bis0 ;
                                //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                                                                                    
                            }
                            else if ( pow_y == basis_degree ){
                                //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                                size_t j_y_bis1 = N - iy2_bis;
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * py_bis1  ;
                                //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                            }
                            else if ( pow_y == 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                size_t j_y_bis1 = N - iy_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * (py_bis0 - 2*py_bis1) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else if ( pow_y == basis_degree - 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_y_bis0 = N - iy_bis;
                                size_t j_y_bis1 = N - iy2_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,iy_bis);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * (py_bis1 - 2*py_bis0) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else{
                                size_t j_y_bis0 = N - pow_y;
                                size_t j_y_bis1 = N - iy_bis;
                                size_t j_y_bis2 = N - iy2_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto coeff_n_y_bis2 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                auto py_bis2 = scaling_y_bis * coeff_n_y_bis2 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis2);
                                                                                    
                                                                                
                                dy = coeff_d2y * ( py_bis0 - 2*py_bis1 + py_bis2 );
                                                    
                            }
                                                                                
            #endif
                            //std::cout<<"pos "<<pos<<std::endl;
                            //std::cout<<"dx "<<dx<<std::endl;
                            //std::cout<<"dy "<<dy<<std::endl;
                            //std::cout<<"px "<<px<<std::endl;
                            //std::cout<<"py "<<py<<std::endl;
                            
                            ret(pos) = dx*py + px*dy;
                            pos++;
                                                                                
                                                        
                           
                            
                        }
                        
                        // face 1 (RIGHT) size_t pow_x = N_partial;
                        for (size_t pow_y = starter + 1; pow_y <= N_partial - 1 ; pow_y++)
                        {
                            size_t pow_x = N_partial;
                            VT ix_bis = pow_x-1; // element i-1 for derivative
                            
                            VT iy_bis = pow_y-1; // element i-1 for derivative
                            VT ix2_bis = pow_x-2; // element i-2 for derivative
                            VT iy2_bis = pow_y-2; // element i-2 for derivative
                                                                                                        
        #ifdef POWER_CACHE
                                                                                                    
         
                            std::cout<<"Bernstein Gradient, basis >=2: POWER_CACHE"<<std::endl;
                            auto px = power_cache[2*pow_x];
                            auto py = power_cache[2*pow_y+1];
                            //std::cout<<"px POWER_CACHE "<<px<<std::endl;
                            if ( pow_x == 0 )
                                auto dx = coeff_d2x * power_cache_bis[2*pow_x] ;
                            else if ( pow_x == basis_degree )
                                auto dx = coeff_d2x * power_cache_bis[2*ix2_bis] ;
                            else if ( pow_x == 1 )
                                auto dx = coeff_d2x * ( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] );
                            else if ( pow_x == basis_degree - 1 )
                                auto dx = coeff_d2x * ( power_cache_bis[2*ix2_bis] -2*power_cache_bis[2*ix_bis] );
                            else
                                auto dx = coeff_d2x*( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] + power_cache_bis[2*ix2_bis] );
                                                                                                        
                                                                            
                            if ( pow_y == 0 )
                                auto dy = coeff_d2y*power_cache_bis[2*pow_y+1] ;
                            else if ( pow_y == basis_degree )
                                auto dy = coeff_d2y*power_cache_bis[2*iy2_bis+1] ;
                            else if ( pow_y == 1 )
                                auto dy = coeff_d2y * ( power_cache_bis[2*pow_y+1] -2*power_cache_bis[2*iy_bis+1] );
                            else if ( pow_y == basis_degree - 1 )
                                auto dy = coeff_d2y * ( power_cache_bis[2*iy2_bis+1] -2*power_cache_bis[2*iy_bis+1] );
                            else
                                auto dy = coeff_d2y*( power_cache_bis[2*pow_y+1] - 2*power_cache_bis[2*iy_bis+1] + power_cache_bis[2*iy2_bis+1] );
                            //std::cout<<"dx POWER_CACHE "<<dx<<std::endl;
                                                                                                        
            #else
                            //std::cout<<"Bernstein Gradient, basis >=2: NOT CACHE"<<std::endl;
                            size_t j_x = basis_degree - pow_x;
                            size_t j_y = basis_degree - pow_y;
                            //std::cout<<"j_x PRE "<<j_x<<std::endl;
                            auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                            auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                            //std::cout<<"coeff_n_x PRE "<<coeff_n_x<<std::endl;
                            //std::cout<<"scaling_x PRE "<<scaling_x<<std::endl;
                            //std::cout<<"pow_x PRE "<<pow_x<<std::endl;
                            //std::cout<<"bx0 PRE "<<bx0<<std::endl;
                            //std::cout<<"bx1 PRE "<<bx1<<std::endl;
                                                                                                        
                            auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                            auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                            //std::cout<<"px PRE "<<px<<std::endl;
                                                                                                    
                            // DERIVATIVES
                            //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                            auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 , px_bis2 = 0.0 , py_bis2 = 0.0 ;
                            VT dx = 0.0 , dy = 0.0 ;
                                                                                                        
                            
                            if ( pow_x == 0 ){
                                //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                dx = coeff_d2x * px_bis0 ;
                                //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                                                                                                            
                            }
                            else if ( pow_x == basis_degree ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis1 = N - ix2_bis;
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * px_bis1 ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else if ( pow_x == 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                size_t j_x_bis1 = N - ix_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * (px_bis0 - 2*px_bis1) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                                                                                                        
                            else if ( pow_x == basis_degree - 1 ){
                                                                        
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - ix_bis;
                                size_t j_x_bis1 = N - ix2_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,ix_bis);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * (px_bis1 - 2*px_bis0) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                                    
                            }
                            else{
                                size_t j_x_bis0 = N - pow_x;
                                size_t j_x_bis1 = N - ix_bis;
                                size_t j_x_bis2 = N - ix2_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto coeff_n_x_bis2 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                auto px_bis2 = scaling_x_bis * coeff_n_x_bis2 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis2);
                                                                                                        
                                dx = coeff_d2x * (px_bis0 -2*px_bis1 +px_bis2);
                                //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                            }
                            
                            if ( pow_y == 0 ){
                                //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                                dy = coeff_d2y *  py_bis0 ;
                                //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                            }
                            
                            else if ( pow_y == basis_degree ){
                                //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                                size_t j_y_bis1 = N - iy2_bis;
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * py_bis1  ;
                                //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                            }
                            else if ( pow_y == 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                size_t j_y_bis1 = N - iy_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * (py_bis0 - 2*py_bis1) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else if ( pow_y == basis_degree - 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_y_bis0 = N - iy_bis;
                                size_t j_y_bis1 = N - iy2_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,iy_bis);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * (py_bis1 - 2*py_bis0) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else{
                                size_t j_y_bis0 = N - pow_y;
                                size_t j_y_bis1 = N - iy_bis;
                                size_t j_y_bis2 = N - iy2_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto coeff_n_y_bis2 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                auto py_bis2 = scaling_y_bis * coeff_n_y_bis2 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis2);
                                dy = coeff_d2y * ( py_bis0 - 2*py_bis1 + py_bis2 );
                                                                        
                            }
                                                                                                        
            #endif
                            //std::cout<<"pos "<<pos<<std::endl;
                            //std::cout<<"dx "<<dx<<std::endl;
                            //std::cout<<"dy "<<dy<<std::endl;
                            //std::cout<<"px "<<px<<std::endl;
                            //std::cout<<"py "<<py<<std::endl;
                                            
                            ret(pos) = dx*py + px*dy;
                            pos++;
                        }
                        
                        // face 2 (TOP)
                        for (size_t pow_x = N_partial - 1; pow_x >= starter + 1 ; pow_x--)
                        {
                            size_t pow_y = N_partial;
                            VT ix_bis = pow_x-1; // element i-1 for derivative
                            VT iy_bis = pow_y-1; // element i-1 for derivative
                            VT ix2_bis = pow_x-2; // element i-2 for derivative
                            VT iy2_bis = pow_y-2; // element i-2 for derivative
                                                                                                                            
        #ifdef POWER_CACHE
                                                                                                                        
                            std::cout<<"Bernstein Gradient, basis >=2: POWER_CACHE"<<std::endl;
                            auto px = power_cache[2*pow_x];
                            auto py = power_cache[2*pow_y+1];
                            //std::cout<<"px POWER_CACHE "<<px<<std::endl;
                            if ( pow_x == 0 )
                                auto dx = coeff_d2x * power_cache_bis[2*pow_x] ;
                            else if ( pow_x == basis_degree )
                                auto dx = coeff_d2x * power_cache_bis[2*ix2_bis] ;
                            else if ( pow_x == 1 )
                                auto dx = coeff_d2x * ( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] );
                            else if ( pow_x == basis_degree - 1 )
                                auto dx = coeff_d2x * ( power_cache_bis[2*ix2_bis] -2*power_cache_bis[2*ix_bis] );
                            else
                                auto dx = coeff_d2x*( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] + power_cache_bis[2*ix2_bis] );
                                                                                                                            
                                                                                                
                            if ( pow_y == 0 )
                                auto dy = coeff_d2y*power_cache_bis[2*pow_y+1] ;
                            else if ( pow_y == basis_degree )
                                auto dy = coeff_d2y*power_cache_bis[2*iy2_bis+1] ;
                            else if ( pow_y == 1 )
                                auto dy = coeff_d2y * ( power_cache_bis[2*pow_y+1] -2*power_cache_bis[2*iy_bis+1] );
                            else if ( pow_y == basis_degree - 1 )
                                auto dy = coeff_d2y * ( power_cache_bis[2*iy2_bis+1] -2*power_cache_bis[2*iy_bis+1] );
                            else
                                auto dy = coeff_d2y*( power_cache_bis[2*pow_y+1] - 2*power_cache_bis[2*iy_bis+1] + power_cache_bis[2*iy2_bis+1] );
                            //std::cout<<"dx POWER_CACHE "<<dx<<std::endl;
                                                                                                                        
            #else
             
                            //std::cout<<"Bernstein Gradient, basis >=2: NOT CACHE"<<std::endl;
                            size_t j_x = basis_degree - pow_x;
                            size_t j_y = basis_degree - pow_y;
                            //std::cout<<"j_x PRE "<<j_x<<std::endl;
                            auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                            auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                            //std::cout<<"coeff_n_x PRE "<<coeff_n_x<<std::endl;
                            //std::cout<<"scaling_x PRE "<<scaling_x<<std::endl;
                            //std::cout<<"pow_x PRE "<<pow_x<<std::endl;
                            //std::cout<<"bx0 PRE "<<bx0<<std::endl;
                            //std::cout<<"bx1 PRE "<<bx1<<std::endl;
                                                                                                                            
                            auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                            auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                            //std::cout<<"px PRE "<<px<<std::endl;
                                                                                                                    
                            //DERIVATIVES
                            //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                            auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 , px_bis2 = 0.0 , py_bis2 = 0.0 ;
                            VT dx = 0.0 , dy = 0.0 ;
                                                                                                                            
                                                
                            if ( pow_x == 0 ){
                                //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                dx = coeff_d2x * px_bis0 ;
                                //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                            }
                             
                            else if ( pow_x == basis_degree ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis1 = N - ix2_bis;
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * px_bis1 ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                             
                            else if ( pow_x == 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                size_t j_x_bis1 = N - ix_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * (px_bis0 - 2*px_bis1) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                                
                            else if ( pow_x == basis_degree - 1 ){
                                                                                            
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - ix_bis;
                                size_t j_x_bis1 = N - ix2_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,ix_bis);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * (px_bis1 - 2*px_bis0) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                                                    
                            }
                             
                            else{
                                size_t j_x_bis0 = N - pow_x;
                                size_t j_x_bis1 = N - ix_bis;
                                size_t j_x_bis2 = N - ix2_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto coeff_n_x_bis2 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                auto px_bis2 = scaling_x_bis * coeff_n_x_bis2 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis2);
                                                                                                                        
                                dx = coeff_d2x * (px_bis0 -2*px_bis1 +px_bis2);
                                //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                            }
                                                
                            if ( pow_y == 0 ){
                                //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                                dy = coeff_d2y *  py_bis0 ;
                                //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                            }
                                                
                            else if ( pow_y == basis_degree ){
                                //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                                size_t j_y_bis1 = N - iy2_bis;
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * py_bis1  ;
                                //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                            }
                            else if ( pow_y == 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                size_t j_y_bis1 = N - iy_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * (py_bis0 - 2*py_bis1) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            
                            else if ( pow_y == basis_degree - 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_y_bis0 = N - iy_bis;
                                size_t j_y_bis1 = N - iy2_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,iy_bis);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * (py_bis1 - 2*py_bis0) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else{
                                size_t j_y_bis0 = N - pow_y;
                                size_t j_y_bis1 = N - iy_bis;
                                size_t j_y_bis2 = N - iy2_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto coeff_n_y_bis2 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                auto py_bis2 = scaling_y_bis * coeff_n_y_bis2 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis2);
                                dy = coeff_d2y * ( py_bis0 - 2*py_bis1 + py_bis2 );
                                                                                            
                            }
                                                                                                                            
            #endif
            
                            //std::cout<<"pos "<<pos<<std::endl;
                            //std::cout<<"dx "<<dx<<std::endl;
                            //std::cout<<"dy "<<dy<<std::endl;
                            //std::cout<<"px "<<px<<std::endl;
                            //std::cout<<"py "<<py<<std::endl;
                                                                
                            ret(pos) = dx*py + px*dy;
                            pos++;
                        }
                        
                        
                        // face 3 (LEFT)
                        
                        for (size_t pow_y = N_partial - 1; pow_y >= starter + 1 ; pow_y--)
                        {
                            size_t pow_x = starter;
                            VT ix_bis = pow_x-1; // element i-1 for derivative
                            VT iy_bis = pow_y-1; // element i-1 for derivative
                            VT ix2_bis = pow_x-2; // element i-2 for derivative
                            VT iy2_bis = pow_y-2; // element i-2 for derivative
                                                                                                                            
        #ifdef POWER_CACHE
                                                                                                                        
                             
                            std::cout<<"Bernstein Gradient, basis >=2: POWER_CACHE"<<std::endl;
                            auto px = power_cache[2*pow_x];
                            auto py = power_cache[2*pow_y+1];
                            //std::cout<<"px POWER_CACHE "<<px<<std::endl;
                            if ( pow_x == 0 )
                                auto dx = coeff_d2x * power_cache_bis[2*pow_x] ;
                            else if ( pow_x == basis_degree )
                                auto dx = coeff_d2x * power_cache_bis[2*ix2_bis] ;
                            else if ( pow_x == 1 )
                                auto dx = coeff_d2x * ( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] );
                            else if ( pow_x == basis_degree - 1 )
                                auto dx = coeff_d2x * ( power_cache_bis[2*ix2_bis] -2*power_cache_bis[2*ix_bis] );
                            else
                                auto dx = coeff_d2x*( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] + power_cache_bis[2*ix2_bis] );
                                                                                                                            
                                                                                                
                            if ( pow_y == 0 )
                                auto dy = coeff_d2y*power_cache_bis[2*pow_y+1] ;
                            else if ( pow_y == basis_degree )
                                auto dy = coeff_d2y*power_cache_bis[2*iy2_bis+1] ;
                            else if ( pow_y == 1 )
                                auto dy = coeff_d2y * ( power_cache_bis[2*pow_y+1] -2*power_cache_bis[2*iy_bis+1] );
                            else if ( pow_y == basis_degree - 1 )
                                auto dy = coeff_d2y * ( power_cache_bis[2*iy2_bis+1] -2*power_cache_bis[2*iy_bis+1] );
                            else
                                auto dy = coeff_d2y*( power_cache_bis[2*pow_y+1] - 2*power_cache_bis[2*iy_bis+1] + power_cache_bis[2*iy2_bis+1] );
                            //std::cout<<"dx POWER_CACHE "<<dx<<std::endl;
                                                                                                                            
            #else
                            
                            //std::cout<<"Bernstein Gradient, basis >=2: NOT CACHE"<<std::endl;
                            size_t j_x = basis_degree - pow_x;
                            size_t j_y = basis_degree - pow_y;
                            //std::cout<<"j_x PRE "<<j_x<<std::endl;
                            auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                            auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                            //std::cout<<"coeff_n_x PRE "<<coeff_n_x<<std::endl;
                            //std::cout<<"scaling_x PRE "<<scaling_x<<std::endl;
                            //std::cout<<"pow_x PRE "<<pow_x<<std::endl;
                            //std::cout<<"bx0 PRE "<<bx0<<std::endl;
                            //std::cout<<"bx1 PRE "<<bx1<<std::endl;
                                                                                                                            
                            auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                            auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                            //std::cout<<"px PRE "<<px<<std::endl;
                                                                                                                        
                            // DERIVATIVES
                            //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                            auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 , px_bis2 = 0.0 , py_bis2 = 0.0 ;
                            VT dx = 0.0 , dy = 0.0 ;
                                                                                                                            
                                                
                            if ( pow_x == 0 ){
                                //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                dx = coeff_d2x * px_bis0 ;
                                //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                            }
                            
                            else if ( pow_x == basis_degree ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis1 = N - ix2_bis;
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * px_bis1 ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else if ( pow_x == 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - pow_x;
                                size_t j_x_bis1 = N - ix_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * (px_bis0 - 2*px_bis1) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                                 
                            else if ( pow_x == basis_degree - 1 ){
                                                                                            
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_x_bis0 = N - ix_bis;
                                size_t j_x_bis1 = N - ix2_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,ix_bis);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                                dx = coeff_d2x * (px_bis1 - 2*px_bis0) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                                                        
                            }
                            else{
                                size_t j_x_bis0 = N - pow_x;
                                size_t j_x_bis1 = N - ix_bis;
                                size_t j_x_bis2 = N - ix2_bis;
                                auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                                auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                                auto coeff_n_x_bis2 = binomial_coeff_fx(N,ix2_bis);
                                auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                                auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                                auto px_bis2 = scaling_x_bis * coeff_n_x_bis2 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis2);
                                                                                                                            
                                dx = coeff_d2x * (px_bis0 -2*px_bis1 +px_bis2);
                                //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                            }
                                                
                            if ( pow_y == 0 ){
                                //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                                dy = coeff_d2y *  py_bis0 ;
                                //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                            }
                                                
                            else if ( pow_y == basis_degree ){
                                //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                                size_t j_y_bis1 = N - iy2_bis;
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * py_bis1  ;
                                //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                            }
                            
                            else if ( pow_y == 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_y_bis0 = N - pow_y;
                                size_t j_y_bis1 = N - iy_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * (py_bis0 - 2*py_bis1) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                                    
                            else if ( pow_y == basis_degree - 1 ){
                                //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                                size_t j_y_bis0 = N - iy_bis;
                                size_t j_y_bis1 = N - iy2_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,iy_bis);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                                dy = coeff_d2y * (py_bis1 - 2*py_bis0) ;
                                //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                            }
                            else{
                                size_t j_y_bis0 = N - pow_y;
                                size_t j_y_bis1 = N - iy_bis;
                                size_t j_y_bis2 = N - iy2_bis;
                                auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                                auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                                auto coeff_n_y_bis2 = binomial_coeff_fx(N,iy2_bis);
                                auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                                auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                                auto py_bis2 = scaling_y_bis * coeff_n_y_bis2 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis2);
                                dy = coeff_d2y * ( py_bis0 - 2*py_bis1 + py_bis2 );
                                                                                            
                            }
                                                                                                                            
                #endif
                            //std::cout<<"pos "<<pos<<std::endl;
                            //std::cout<<"dx "<<dx<<std::endl;
                            //std::cout<<"dy "<<dy<<std::endl;
                            //std::cout<<"px "<<px<<std::endl;
                            //std::cout<<"py "<<py<<std::endl;
                                                                
                            ret(pos) = dx*py + px*dy;
                            pos++;
                        }
                        
                        N_partial--;
                        starter++;
                        internal_bases -= 2;
                        
                    }
                    
                    // for B_k, with k even, I.E. odd number of bases: there is a central one.
                    if( basis_degree % 2 == 0 )
                    {
                        //  std::cout<<"N is "<<N_partial<< " and starter is "<<starter<<std::endl;
                        assert( N_partial == starter );
                        size_t pow_x = starter ;
                        size_t pow_y = starter ;
                                         
                        VT ix_bis = pow_x-1; // element i-1 for derivative
                        VT iy_bis = pow_y-1; // element i-1 for derivative
                        VT ix2_bis = pow_x-2; // element i-2 for derivative
                        VT iy2_bis = pow_y-2; // element i-2 for derivative
                                                                                                                                            
        #ifdef POWER_CACHE
                                                                                                                                        
                                             
                        std::cout<<"Bernstein Gradient, basis >=2: POWER_CACHE"<<std::endl;
                        auto px = power_cache[2*pow_x];
                        auto py = power_cache[2*pow_y+1];
                        //std::cout<<"px POWER_CACHE "<<px<<std::endl;
                        
                        if ( pow_x == 0 )
                            auto dx = coeff_d2x * power_cache_bis[2*pow_x] ;
                        else if ( pow_x == basis_degree )
                            auto dx = coeff_d2x * power_cache_bis[2*ix2_bis] ;
                        else if ( pow_x == 1 )
                            auto dx = coeff_d2x * ( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] );
                        else if ( pow_x == basis_degree - 1 )
                            auto dx = coeff_d2x * ( power_cache_bis[2*ix2_bis] -2*power_cache_bis[2*ix_bis] );
                        else
                            auto dx = coeff_d2x*( power_cache_bis[2*pow_x] -2*power_cache_bis[2*ix_bis] + power_cache_bis[2*ix2_bis] );
                       
                        if ( pow_y == 0 )
                            auto dy = coeff_d2y*power_cache_bis[2*pow_y+1] ;
                        else if ( pow_y == basis_degree )
                            auto dy = coeff_d2y*power_cache_bis[2*iy2_bis+1] ;
                        else if ( pow_y == 1 )
                            auto dy = coeff_d2y * ( power_cache_bis[2*pow_y+1] -2*power_cache_bis[2*iy_bis+1] );
                        else if ( pow_y == basis_degree - 1 )
                            auto dy = coeff_d2y * ( power_cache_bis[2*iy2_bis+1] -2*power_cache_bis[2*iy_bis+1] );
                        else
                            auto dy = coeff_d2y*( power_cache_bis[2*pow_y+1] - 2*power_cache_bis[2*iy_bis+1] + power_cache_bis[2*iy2_bis+1] );
                        //std::cout<<"dx POWER_CACHE "<<dx<<std::endl;
                    
    #else
                        //std::cout<<"Bernstein Gradient, basis >=2: NOT CACHE"<<std::endl;
                        size_t j_x = basis_degree - pow_x;
                        size_t j_y = basis_degree - pow_y;
                        //std::cout<<"j_x PRE "<<j_x<<std::endl;
                        auto coeff_n_x = binomial_coeff_fx(basis_degree,pow_x);
                        auto coeff_n_y = binomial_coeff_fx(basis_degree,pow_y);
                        //std::cout<<"coeff_n_x PRE "<<coeff_n_x<<std::endl;
                        //std::cout<<"scaling_x PRE "<<scaling_x<<std::endl;
                        //std::cout<<"pow_x PRE "<<pow_x<<std::endl;
                        //std::cout<<"bx0 PRE "<<bx0<<std::endl;
                        //std::cout<<"bx1 PRE "<<bx1<<std::endl;
                                                                                                                                            
                        auto px = scaling_x *coeff_n_x* iexp_pow(bx0, pow_x)*iexp_pow(bx1, j_x);
                        auto py = scaling_y *coeff_n_y* iexp_pow(by0, pow_y)*iexp_pow(by1, j_y);
                        //std::cout<<"px PRE "<<px<<std::endl;
                                                                                                                                        
                        // DERIVATIVES
                        //std::cout<<"pow_x= "<<pow_x<<"pow_y= "<<pow_y<<std::endl;
                        auto px_bis0 = 0.0 , py_bis0 = 0.0 , px_bis1 = 0.0 , py_bis1 = 0.0 , px_bis2 = 0.0 , py_bis2 = 0.0 ;
                        VT dx = 0.0 , dy = 0.0 ;
                                                                                                
                        if ( pow_x == 0 ){
                            //std::cout<<"into bottom pow_x == 0: HO b_i^K-1(x)"<<std::endl;
                            size_t j_x_bis0 = N - pow_x;
                            auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                            auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                            dx = coeff_d2x * px_bis0 ;
                            //auto dx = -coeff_dx*power_cache_bis[2*pow_x] ;
                        }
                                            
                        else if ( pow_x == basis_degree ){
                            //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                            size_t j_x_bis1 = N - ix2_bis;
                            auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                            auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                            dx = coeff_d2x * px_bis1 ;
                            //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                        }
                        else if ( pow_x == 1 ){
                            //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                            size_t j_x_bis0 = N - pow_x;
                            size_t j_x_bis1 = N - ix_bis;
                            auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                            auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                            auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                            auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                            dx = coeff_d2x * (px_bis0 - 2*px_bis1) ;
                            //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                        }
                                                 
                        else if ( pow_x == basis_degree - 1 ){
                            //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                            size_t j_x_bis0 = N - ix_bis;
                            size_t j_x_bis1 = N - ix2_bis;
                            auto coeff_n_x_bis0 = binomial_coeff_fx(N,ix_bis);
                            auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix2_bis);
                            auto px_bis0 = scaling_x_bis *  coeff_n_x_bis0 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis0);
                            auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis1);
                            dx = coeff_d2x * (px_bis1 - 2*px_bis0) ;
                            //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                                                                                    
                        }
                        else{
                            size_t j_x_bis0 = N - pow_x;
                            size_t j_x_bis1 = N - ix_bis;
                            size_t j_x_bis2 = N - ix2_bis;
                            auto coeff_n_x_bis0 = binomial_coeff_fx(N,pow_x);
                            auto coeff_n_x_bis1 = binomial_coeff_fx(N,ix_bis);
                            auto coeff_n_x_bis2 = binomial_coeff_fx(N,ix2_bis);
                            auto px_bis0 = scaling_x_bis * coeff_n_x_bis0 * iexp_pow(bx0, pow_x) * iexp_pow(bx1, j_x_bis0);
                            auto px_bis1 = scaling_x_bis * coeff_n_x_bis1 * iexp_pow(bx0, ix_bis) * iexp_pow(bx1, j_x_bis1);
                            auto px_bis2 = scaling_x_bis * coeff_n_x_bis2 * iexp_pow(bx0, ix2_bis) * iexp_pow(bx1, j_x_bis2);
                          
                            dx = coeff_d2x * (px_bis0 -2*px_bis1 +px_bis2);
                            //auto dx = coeff_dx*( power_cache_bis[2*ix_bis]-power_cache_bis[2*pow_x]);
                        }
                                                                
                        
                        if ( pow_y == 0 ){
                            //std::cout<<"into bottom pow_y = 0: HO b_i^K-1(y)"<<std::endl;
                            size_t j_y_bis0 = N - pow_y;
                            auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                            auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                            // auto dy = coeff_dy * ( py_bis1 - py_bis0 );
                            dy = coeff_d2y *  py_bis0 ;
                            //auto dy = -coeff_dy*power_cache_bis[2*pow_y+1] ;
                        }
                                                                
                        else if ( pow_y == basis_degree ){
                            //std::cout<<"into bottom pow_y != basis_degree: HO b_i-1^K-1(y)"<<std::endl;
                            size_t j_y_bis1 = N - iy2_bis;
                            auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                            auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                            dy = coeff_d2y * py_bis1  ;
                            //auto dy = coeff_dy*power_cache_bis[2*iy_bis+1] ;
                        }
                                        
                        else if ( pow_y == 1 ){
                            //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                            size_t j_y_bis0 = N - pow_y;
                            size_t j_y_bis1 = N - iy_bis;
                            auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                            auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                            auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                            auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                            dy = coeff_d2y * (py_bis0 - 2*py_bis1) ;
                            //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                        }
                                                    
                        else if ( pow_y == basis_degree - 1 ){
                            //std::cout<<"into bottom pow_x != basis_degree: HO b_i-1^K-1(x)"<<std::endl;
                            size_t j_y_bis0 = N - iy_bis;
                            size_t j_y_bis1 = N - iy2_bis;
                            auto coeff_n_y_bis0 = binomial_coeff_fx(N,iy_bis);
                            auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy2_bis);
                            auto py_bis0 = scaling_y_bis *  coeff_n_y_bis0 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis0);
                            auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis1);
                            dy = coeff_d2y * (py_bis1 - 2*py_bis0) ;
                            //auto dx = coeff_dx*power_cache_bis[2*ix_bis] ;
                        }
                        
                        else{
                            size_t j_y_bis0 = N - pow_y;
                            size_t j_y_bis1 = N - iy_bis;
                            size_t j_y_bis2 = N - iy2_bis;
                            auto coeff_n_y_bis0 = binomial_coeff_fx(N,pow_y);
                            auto coeff_n_y_bis1 = binomial_coeff_fx(N,iy_bis);
                            auto coeff_n_y_bis2 = binomial_coeff_fx(N,iy2_bis);
                            auto py_bis0 = scaling_y_bis * coeff_n_y_bis0 * iexp_pow(by0, pow_y) * iexp_pow(by1, j_y_bis0);
                            auto py_bis1 = scaling_y_bis * coeff_n_y_bis1 * iexp_pow(by0, iy_bis) * iexp_pow(by1, j_y_bis1);
                            auto py_bis2 = scaling_y_bis * coeff_n_y_bis2 * iexp_pow(by0, iy2_bis) * iexp_pow(by1, j_y_bis2);
                            dy = coeff_d2y * ( py_bis0 - 2*py_bis1 + py_bis2 );
                            
                        }
                                                                                                
    #endif
                        //std::cout<<"pos "<<pos<<std::endl;
                        //std::cout<<"dx "<<dx<<std::endl;
                        //std::cout<<"dy "<<dy<<std::endl;
                        //std::cout<<"px "<<px<<std::endl;
                        //std::cout<<"py "<<py<<std::endl;
                                                                                
                        ret(pos) = dx*py + px*dy;
                        pos++;

                       
                        
                    }
                }

                  //  std::cout<<"GRADIENTE: pos is"<<pos<<" and basis size is "<< basis_size<<std::endl;
                    assert(pos == basis_size);

                    return ret;
                

    }
    
    
    
    size_t size() const
    {
        return basis_size;
    }

    size_t degree() const
    {
        return basis_degree;
    }
    
    static size_t size(size_t degree)
    {
        return (degree+1)*(degree+1);
    }
};



template<typename Mesh, typename T >
class cell_basis_Lagrangian
{
    typedef typename Mesh::coordinate_type  coordinate_type;
    typedef typename Mesh::point_type       point_type;

    point_type          cell_bar;
    coordinate_type     cell_h;
    size_t              basis_degree, basis_size;
    std::vector<point<T, 2> >          nodes;
    //std::vector<size_t>         indeces;
    
public:
    cell_basis_Lagrangian(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
    {
        nodes           = equidistriduted_nodes<T,Mesh>(msh, cl, degree);
       
        basis_degree    = degree;
        basis_size      = (basis_degree+1)*(basis_degree+1);
    }
    
    /*
    cell_basis_Lagrangian(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, const std::vector<size_t>& indeces)
    {
        
        nodes = equidistriduted_nodes_subcell<T,Mesh>(msh, cl, degree, indeces);
        basis_degree    = degree;
        basis_size      = (basis_degree+1)*(basis_degree+1);
        
    }
    */

    Matrix<T, Dynamic, 1>
    eval_basis(const point_type& pt)
    {
        Matrix<T, Dynamic, 1> rety = Matrix<T, Dynamic, 1>::Zero(basis_size);
        Matrix<T, Dynamic, 1> retx = Matrix<T, Dynamic, 1>::Zero(basis_size);
        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(basis_size);

        // Per la y, trovo la colonna facendo col = l%(degree+1)
        // scorro su tutta la colonna tmpy = [col:(degree+1): col+(degree+1)*degree]
        // faccio base bl moltiplicando tutti tranne quando tmpy = l
        for ( size_t l = 0; l < basis_size ; l++ )
        {
            size_t col = l%(basis_degree+1);
            T bl = 1.0;
            for (size_t tmpy = col; tmpy <= col+(basis_degree+1)*basis_degree; tmpy+=(basis_degree+1))
            {
                if (tmpy!=l)
                {
                    bl *= ( ( pt.y() - (nodes.at(tmpy)).y() )/ ( (nodes.at(l)).y() - (nodes.at(tmpy)).y() ) );
                }
            }
            rety(l) = bl;
        }
        
        // Per la x, trovo la riga facendo riga = floor(k/(degree+1))
        //scorro su tutta la riga tmpx = [(degree+1)*riga: (degree+1)*(riga+1)-1]
        // faccio base bk moltiplicando tutti tranne quando tmpx = k
         
        for (size_t k = 0 ; k < basis_size ; k++ )
        {
            T bk = 1.0;
            size_t row=floor( k/(basis_degree+1) );
            for (size_t tmpx = (basis_degree+1)*row; tmpx <= (basis_degree+1)*(row+1)-1; tmpx++)
            {
                if (tmpx!=k) {
                    bk *= ( ( pt.x() - (nodes.at(tmpx)).x() )/ ( (nodes.at(k)).x() - (nodes.at(tmpx)).x() ) );
                }
            }
            retx(k) = bk;
        }
        
        for (size_t i = 0; i<basis_size; i++)
        {
            ret(i) = rety(i)*retx(i);
        }
        return ret;
        
    }

    Matrix<T, Dynamic, 2>
    eval_gradients(const point_type& pt)
    {
        // Modified Yves Daoust Algorithm (https://math.stackexchange.com/questions/809927/first-derivative-of-lagrange-polynomial)
        
        Matrix<T, Dynamic, 2> ret = Matrix<T, Dynamic, 2>::Zero(basis_size, 2);
       
        Matrix<T, Dynamic, 1> rety = Matrix<T, Dynamic, 1>::Zero(basis_size);
        Matrix<T, Dynamic, 1> retx = Matrix<T, Dynamic, 1>::Zero(basis_size);
        Matrix<T, Dynamic, 1> sy = Matrix<T, Dynamic, 1>::Zero(basis_size);
        Matrix<T, Dynamic, 1> sx = Matrix<T, Dynamic, 1>::Zero(basis_size);


        // for each l, b_l(y)' = {sum(tmpy!=l)[prod(jy!=l,jy!=tmpy)[x-x_jy]]}/prod(tmpy!=l)[x_l-x_tmpy]
        
        for ( size_t l = 0; l < basis_size ; l++ )
        {
            size_t col = l%(basis_degree+1);
            T bl = 1.0 , bl_der = 1.0 ;
            T sumy = 0.0;
            for (size_t tmpy = col; tmpy <= col+(basis_degree+1)*basis_degree; tmpy+=(basis_degree+1))
            {
                 T sumyy = 1.0 ;
                if (tmpy!=l)
                {

                    bl *= ( ( pt.y() - (nodes.at(tmpy)).y() )/ ( (nodes.at(l)).y() - (nodes.at(tmpy)).y() ) );

                    bl_der *= ( 1.0/ ( (nodes.at(l)).y() - (nodes.at(tmpy)).y() ) );
                    for ( size_t jy = col; jy <= col+(basis_degree+1)*basis_degree; jy+=(basis_degree+1) )
                   {
                        if (jy!=tmpy && jy!=l)
                        {
                            sumyy *= ( pt.y()-(nodes.at(jy)).y() );
                        }
                    }
                    sumy +=sumyy;
                }
            }
            rety(l) = bl;
            sy(l) = bl_der*sumy;
        }
        
        // For the x-derivative of b_k(x), same procedure of b_l(y)'
         
        for (size_t k = 0 ; k < basis_size ; k++ )
        {
            size_t row=floor( k/(basis_degree+1) );
            T bk = 1.0 , bk_der = 1.0 ;
            T sumx = 0.0;
            for (size_t tmpx = (basis_degree+1)*row; tmpx <= (basis_degree+1)*(row+1)-1; tmpx++)
            {
                T sumxx = 1.0 ;
                
                if (tmpx!=k) {
                    
                    bk *= ( ( pt.x() - (nodes.at(tmpx)).x() )/ ( (nodes.at(k)).x() - (nodes.at(tmpx)).x() ) );
                    bk_der *= ( 1.0/ ( (nodes.at(k)).x() - (nodes.at(tmpx)).x() ) );
                    for (size_t jx = (basis_degree+1)*row; jx <= (basis_degree+1)*(row+1)-1; jx++)
                    {
                        if (jx!=tmpx && jx!=k)
                        {
                            sumxx *= ( pt.x()-(nodes.at(jx)).x() );
                        }
                    }
                    sumx += sumxx;
                    
                }
            }
            retx(k) = bk;
            sx(k) = bk_der*sumx;
        }
        
        for (size_t i = 0; i<basis_size; i++)
        {
            ret(i,0) = rety(i)*sx(i);
            ret(i,1) = retx(i)*sy(i);
            
        }
        return ret;
        
    }

    size_t size() const
    {
        return basis_size;
    }

    size_t degree() const
    {
        return basis_degree;
    }

    static size_t size(size_t degree)
    {
        return (degree+1)*(degree+1);
    }
};


template<typename T>
std::vector<point<T,1> >
reference_nodes(size_t degree)
{
    auto comp_degree = degree + 1;

    size_t reqd_nodes = comp_degree;

    std::vector<point<T,1> > ret;
    ret.reserve(reqd_nodes);

    point<T,1>  qp;
    T           a1, a2;
    T           delta_x;
    switch(reqd_nodes)
    {
        case 1:
            qp = point<T,1>({0.0});
            ret.push_back(qp);
            return ret;

        case 2:
            qp = point<T,1>({ 1.0 });
            ret.push_back( -qp );
            ret.push_back( qp );
            return ret;

        case 3:
            qp = point<T,1>({ 1.0 });
            ret.push_back( -qp );
            ret.push_back( qp );
            qp = point<T,1>({0.0});
            ret.push_back( qp );
            return ret;

        case 4:
            a1 = 1.0/3.0;
            qp = point<T,1>({ 1.0 });
            ret.push_back( -qp );
            ret.push_back( qp );
            qp = point<T,1>({ a1 });
            ret.push_back( -qp );
            ret.push_back( qp );
            return ret;

        case 5:
            // Be carefull in what order data is inserted in ret!
            // In Gauss Legendre the first one was 0.0, now is the last one
            a2 = 0.5;
            a1 = 1.0;
            qp = point<T,1>({ a1 });
            ret.push_back( -qp );
            ret.push_back( qp );

            qp = point<T,1>({ a2 });
            ret.push_back( -qp );
            ret.push_back( qp );

            qp = point<T,1>({ 0.0 });
            ret.push_back( qp );

            return ret;

        default:

            delta_x = 1.0/degree;
            a1 = 1.0;
            while (a1>0) {
                qp = point<T,1>({ a1 });
                ret.push_back( -qp );
                ret.push_back( qp );
                a1-=delta_x;

            }
            if(a1==0)
            {
                qp = point<T,1>({0.0});
                ret.push_back( qp );
            }
            return ret;
    }

    return ret;
}

template<typename T,typename Mesh>
std::vector< point<T,2> >
equidistriduted_nodes(const Mesh& msh,
          const typename Mesh::cell_type& cl,
          size_t degree)
{
    typedef typename Mesh::point_type    point_type;

    auto qps = reference_nodes<T>(degree);
   

    auto pts = points(msh, cl);
    
    auto v0 = pts[1] - pts[0];
    auto v1 = pts[2] - pts[1];
    auto v2 = pts[3] - pts[2];
    auto v3 = pts[3] - pts[0];

    std::vector< point<T,2> > ret;

    auto P = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].x() * (1-xi)*(1-eta) +
               0.25 * pts[1].x() * (1+xi)*(1-eta) +
               0.25 * pts[2].x() * (1+xi)*(1+eta) +
               0.25 * pts[3].x() * (1-xi)*(1+eta);
    };

    auto Q = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].y() * (1-xi)*(1-eta) +
               0.25 * pts[1].y() * (1+xi)*(1-eta) +
               0.25 * pts[2].y() * (1+xi)*(1+eta) +
               0.25 * pts[3].y() * (1-xi)*(1+eta);
    };

    for (auto jtor = qps.begin(); jtor != qps.end(); jtor++)
    {
        for (auto itor = qps.begin(); itor != qps.end(); itor++)
        {
            auto qp_x = *itor;
            auto qp_y = *jtor;

            auto xi = qp_x.x();
            auto eta = qp_y.x();

            auto px = P(xi, eta);
            auto py = Q(xi, eta);

            ret.push_back( point_type(px, py) );
        }
    }

    return ret;
}


/*
template<typename T,typename Mesh>
std::vector< point<T,2> >
equidistriduted_nodes_subcell(const Mesh& msh,
          const typename Mesh::cell_type& cl,
          size_t degree, const std::vector<size_t>& indeces)
{
    typedef typename Mesh::point_type    point_type;

    auto qps = reference_nodes<T>(degree);
   
    std::vector<typename Mesh::point_type> pts;
    auto all_points = points(msh, cl);
    for (size_t i = 0 ; i < 4 ; i++) {
        pts.push_back( all_points[indeces[i]] );
    }
        
    auto v0 = pts[1] - pts[0];
    auto v1 = pts[2] - pts[1];
    auto v2 = pts[3] - pts[2];
    auto v3 = pts[3] - pts[0];

    std::vector< point<T,2> > ret;

    auto P = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].x() * (1-xi)*(1-eta) +
               0.25 * pts[1].x() * (1+xi)*(1-eta) +
               0.25 * pts[2].x() * (1+xi)*(1+eta) +
               0.25 * pts[3].x() * (1-xi)*(1+eta);
    };

    auto Q = [&](T xi, T eta) -> T {
        return 0.25 * pts[0].y() * (1-xi)*(1-eta) +
               0.25 * pts[1].y() * (1+xi)*(1-eta) +
               0.25 * pts[2].y() * (1+xi)*(1+eta) +
               0.25 * pts[3].y() * (1-xi)*(1+eta);
    };

    for (auto jtor = qps.begin(); jtor != qps.end(); jtor++)
    {
        for (auto itor = qps.begin(); itor != qps.end(); itor++)
        {
            auto qp_x = *itor;
            auto qp_y = *jtor;

            auto xi = qp_x.x();
            auto eta = qp_y.x();

            auto px = P(xi, eta);
            auto py = Q(xi, eta);

            ret.push_back( point_type(px, py) );
        }
    }

    return ret;
}

*/

/*

template<typename T , typename Mesh>
size_t
cells_counter(const Mesh& msh , const typename Mesh::cell_type& cl)
{
    
}
*/


/*
// OLD VERSION OF pt_in_cell
template<typename T, typename Mesh>
bool
pt_in_cell(const Mesh& msh, const point<T,2>& point_to_find, const typename Mesh::cell_type& cl)
{
    bool ret = 0;
    auto pts =points(msh,cl);
    for (size_t i = 0; i < pts.size(); i++)
    {
        if( pts[i].x()>=point_to_find.x() && pts[i].y()>=point_to_find.y() )
            ret = 1;
    }
    return ret;

}
*/
template<typename T, typename Mesh>
bool
pt_in_cell(const Mesh& msh, const point<T,2>& point_to_find, const typename Mesh::cell_type& cl)
{
    auto pts =points(msh,cl);
 
    //std::cout<<"Point to find "<<std::setprecision(15)<<point_to_find.x()<<", "<<point_to_find.y()<<std::endl;

   // std::cout<<"Min x "<<std::setprecision(15)<<pts[0].x()<<", max x "<<pts[1].x()<<std::endl;

    //std::cout<<"Min y "<<std::setprecision(15)<<pts[1].y()<<", max y "<<pts[2].y()<<std::endl;

    T epsilon = 1e-10;
     if( (pts[0].x()-epsilon)<=point_to_find.x() && (pts[1].x()+epsilon)>=point_to_find.x() && (pts[1].y()-epsilon)<=point_to_find.y() && (pts[2].y()+epsilon)>=point_to_find.y() )
         return TRUE;
    else
        return FALSE;
  
}

template<typename T, typename Mesh>
size_t
pt_in_subcell(const Mesh& msh, const point<T,2>& point_to_find, const typename Mesh::cell_type& agglocl)
{
    
    for (auto& offset_subcells : agglocl.user_data.offset_subcells)
    {
        //std::cout<<"OFFSET ORIGINAL CELL "<<offset_subcells<<std::endl;
        auto cl = msh.cells[offset_subcells];
        if( pt_in_cell(msh,point_to_find,cl) )
            return offset_subcells;
    }
    // IF IT ARRIVES HERE, IT DIDN'T FIND THE POINT IN THE CELL.
    std::cout<<"the point doesn't find is "<<point_to_find<<std::endl;
    std::cout<<"IT DIDN'T FIND THE POINT IN SUBCELL "<<offset(msh,agglocl)<<std::endl;
    throw std::invalid_argument("Invalid point-> NOT IN AGGLO_CELL");

}

// SOSPETTO SIA UNITILE
/*
template<typename T , typename Mesh>
std::vector<size_t>
subcell_finder(const Mesh& msh, const point<T,2>& point_to_find, const typename Mesh::cell_type& cl, const mesh_init_params<T>& params)
{
    std::vector<size_t> ret(4);
    auto pts = points(msh,cl);
    auto hx = params.hx();
    auto hy = params.hy();
    bool check_already_found_x = FALSE;
    for (size_t i = 0; i < pts.size()-1; i++)
    {
        if( !check_already_found_x && pts[i].x()!= pts[i+1].x() && pts[i].x()<=point_to_find.x() && pts[i+1].x()>=point_to_find.x() )
        {
            if ( (pts[i].y()+hy)>=point_to_find.y() )
            {
                // point in the first row of subcells
                ret[0]=i;
                ret[1]=i+1;
                check_already_found_x = TRUE;
                for (size_t ii = 0; ii < pts.size(); ii++)
                {
                    if( (pts[i].y()+hy) == pts[ii].y() && pts[i].x() == pts[ii].x() ){
                        ret[2] = ii;
                    }
                    if( (pts[i].y()+hy) == pts[ii].y() && pts[i+1].x() == pts[ii].x() ){
                        ret[3] = ii;
                    }
                }
            }
            else
            {
                // point in the second row of subcells
                check_already_found_x = TRUE;
                for (size_t ii = 0; ii < pts.size(); ii++) // loop to find lower pts
                {
                    if( (pts[i].y()+hy) == pts[ii].y() && pts[i].x() == pts[ii].x() )
                    {
                        ret[0] = ii; // now I look for 3rd point
                        for (size_t j = 0; j < pts.size(); j++) {
                            if( (pts[ii].y()+hy) == pts[j].y() && pts[ii].x() == pts[j].x() ){
                                ret[2] = j;
                            }
                        }
                        
                    } // end of the research for 1st and 3rd points
                    if( (pts[i].y()+hy) == pts[ii].y() && pts[i+1].x() == pts[ii].x() )
                    {
                        ret[1] = ii; // now I look for 4th point
                        for (size_t j = 0; j < pts.size(); j++) {
                            if( (pts[ii].y()+hy) == pts[j].y() && pts[ii+1].x() == pts[j].x() ){
                                ret[3] = j;
                            }
                        }
                        
                    }// end of the research for 2nd and 4th points
                    
                } // END of loop to find lower pts
            }  // END of else "search in the second row of subcells"
            
        } // END of the initial if-> I dont enter if ALREADY CHECKED or consecutive points have: same x or doesn't include pt to find
        
    } // end of the initial loop over the point
    return ret;
}

*/


template<typename Mesh, typename T = typename Mesh::coordinate_type>
Matrix<T, Dynamic, 1>
make_bernstein_local_mass_matrix_lumped(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, size_t di = 2)
{
    
    //auto f_neigh = cl.user_data.f_neighbors;
    //auto d_neigh = cl.user_data.d_neighbors;
    
    cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs, 1);
    //Matrix<T, Dynamic, 1> ret2 = Matrix<T, Dynamic, 1>::Zero(cbs, 1);

    auto qps = integrate(msh, cl, (degree)+2); // integration of order 2k

    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        ret += qp.second * phi ;
        // phi * phi.transpose is (degree+1)^2 x (degree+1)^2 ; qp.second is a scalar
    }
    
    /*
    Matrix<T, Dynamic, 1> ret2 = Matrix<T, Dynamic, 1>::Zero(cbs, 1);
    //Matrix<T, Dynamic, 1> ret2 = Matrix<T, Dynamic, 1>::Zero(cbs, 1);

    auto qps2 = integrate(msh, cl, (degree)+2); // integration of order 2k

    for (auto& qp : qps2)
    {
        auto phi = cb.eval_basis(qp.first);
        ret2 += qp.second * phi ;
        // phi * phi.transpose is (degree+1)^2 x (degree+1)^2 ; qp.second is a scalar
    }
    
    std::cout<<"CHECKING AMOUNT QUADRATURE POINTS IN LOCAL MASS LUMPED"<<'\n'<<ret-ret2<<'\n'<<std::endl;
    std::cout<<"FINE CHECKING AMOUNT QUADRATURE POINTS IN LOCAL MASS LUMPED"<<std::endl;
    
    
    //ret2  =  ret.rowwise().sum(); // sum row; i.e. mi0 + mi1 + mi2 + mi3 , with i = 0 : 3
    */
    return ret;
}


template<typename Mesh, typename T = typename Mesh::coordinate_type>
Matrix<T, Dynamic, Dynamic>
make_bernstein_local_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, size_t di = 1)
{
    
    //auto f_neigh = cl.user_data.f_neighbors;
    //auto d_neigh = cl.user_data.d_neighbors;
    
    cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);
    //Matrix<T, Dynamic, 1> ret2 = Matrix<T, Dynamic, 1>::Zero(cbs, 1);

    auto qps = integrate(msh, cl, 2*(degree+di)); // integration of order 2k

    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        ret += qp.second * phi * phi.transpose();
        // phi * phi.transpose is (degree+1)^2 x (degree+1)^2 ; qp.second is a scalar
    }
    //ret2  =  ret.rowwise().sum(); // sum row; i.e. mi0 + mi1 + mi2 + mi3 , with i = 0 : 3
    
    // CHECK INTEGRATION
    /*
    auto qps2 = integrate(msh, cl, 2*(degree)); // integration of order 2k
    Matrix<T, Dynamic, Dynamic> ret2 = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);
    for (auto& qp : qps2)
    {
        auto phi = cb.eval_basis(qp.first);
        ret2 += qp.second * phi * phi.transpose();
        // phi * phi.transpose is (degree+1)^2 x (degree+1)^2 ; qp.second is a scalar
    }
    std::cout<<"CHECKING AMOUNT QUADRATURE POINTS IN LOCAL MASS MATRIX"<<'\n'<<ret-ret2<<'\n'<<std::endl;
     std::cout<<"FINE CHECKING AMOUNT QUADRATURE POINTS IN LOCAL MASS MATRIX"<<std::endl;
    */
    return ret;
}

template<typename Mesh, typename T = typename Mesh::coordinate_type>
std::pair< Matrix<T, Dynamic, Dynamic> , Matrix<T, Dynamic, Dynamic> >
make_bernstein_local_cij_matrix(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, size_t di = 1)
{
    cell_basis_Bernstein<Mesh,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, Dynamic> ret0 = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);
    Matrix<T, Dynamic, Dynamic> ret1 = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);
    
    // for Q1 , degree = 1-> integration of order 2
    auto qps = integrate(msh, cl, 2*(degree+di) ); // integration of order 2k

    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        
        auto phi_grad = cb.eval_gradients(qp.first);
        
        ret0 += qp.second * phi * ((phi_grad).col(0)).transpose();
    
        ret1 += qp.second * phi * ((phi_grad).col(1)).transpose();
    }
    
    /*
    // CHECKING ORDER INTEGRATION
    Matrix<T, Dynamic, Dynamic> ret1_bis = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);
    
    // for Q1 , degree = 1-> integration of order 2
    auto qps2 = integrate(msh, cl, 2*(degree+di)); // integration of order 2k

    for (auto& qp : qps2)
    {
        auto phi = cb.eval_basis(qp.first);
        
        auto phi_grad = cb.eval_gradients(qp.first);
    
        ret1_bis += qp.second * phi * ((phi_grad).col(1)).transpose();
    }
    
    std::cout<<"THE CHECKING FOR CIJ ORDER IS "<<'\n'<<ret1 - ret1_bis<<std::endl;
    */
    return std::make_pair(ret0,ret1);
}


// Ho due possibilità.. se uso values_bis1 uso una vector notation,
//                      se uso values_bis uso la matrix notation -> USO QUESTA PER ORA, MA I TEST CHE HO FATTO PER N,M=32, 64 RILEVANO CHE SONO PRATICAMENTE LO STESSO
template< typename Mesh , typename Fonction , typename FiniteSpace , typename T = typename Mesh::coordinate_type  >
struct L2projected_level_set_high_order: public level_set<T>
{
     
    T phi_max , phi_min ;
    size_t  last_row_init, last_row_end, number_faces_one_row;
    
    Mesh msh; // Original mesh, NOT agglomerated.
    mesh_init_params<T> params; // mesh parameter
       
    size_t degree_FEM ; // FEM degree
    size_t n_cls ; // #cells
    size_t local_dim; // Local Dimension (degree_FEM+1)*(degree_FEM+1)
    size_t n_vertices ; // #vertices
    size_t      Nx, Ny ; // Number of cells in x and y direciton
    
    // connectivity matrix : for each cell, it stores the global numbering
    std::vector< std::vector<std::pair<size_t,bool>>> connectivity_matrix ;
    std::vector<std::set<size_t>> S_i;
    //std::vector<std::vector<size_t>> connectivity_matrix ;
       
    size_t dim_HHO; // Global dimension Discontinuous framework = Local dimension * #cells
    size_t ndof_FE; // Global dimension FE continuous = #nodes
       
    int mapped = 0 ; // = 0 not mapped ; = 1 mapped ; = 2 inverse mapping
    
    SparseMatrix<T>                 Global_Mass; // Global mass, saved for FEM problem
       
    Eigen::Matrix<T, Dynamic, Dynamic> sol_HHO ; // projection saved in HHO format: cell by cell
    Matrix<T, Dynamic, 1> sol_FEM ; // projection saved in Continuos FE format: global nodes
    Eigen::Matrix<T, Dynamic, 1> vertices; // saving level_set on vertices mesh

    // Local Vandermonde Matrix for interpolation
    Matrix<T, Dynamic, Dynamic> local_vandermonde ;
       
       
    // Assembling into global matrix
    SparseMatrix<T>                 Global_c_term_x; // Global mass, saved for FEM problem
    SparseMatrix<T>                 Global_c_term_y; // Global mass, saved for FEM problem
    std::vector< Triplet<T> >       triplets_c_term_x; // Position elements: Sparse Matrix Notation
    std::vector< Triplet<T> >       triplets_c_term_y; // Position elements: Sparse Matrix Notation
    Matrix<T, Dynamic, 1>      Global_Mass_Lumped; // Global mass, saved for FEM problem

    L2projected_level_set_high_order(const FiniteSpace& fe_data , const Fonction & level_set, const Mesh & msh)
        : degree_FEM(fe_data.order) , local_dim(fe_data.local_ndof), msh(msh), Nx(fe_data.Nx),Ny(fe_data.Ny), params(fe_data.params) , dim_HHO(fe_data.ndof_disc) , n_cls(fe_data.n_cls) ,n_vertices(fe_data.n_vertices) , connectivity_matrix(fe_data.connectivity_matrix) , ndof_FE(fe_data.ndof_FE)
    {
        
        Matrix<T, Dynamic, 1>           RHS;    // Known term
        std::vector< Triplet<T> >       triplets; // Position elements: Sparse Matrix Notation
        
        last_row_init = Ny*(2*Nx+1); // There are 2 faces for each row of cells + Ny
        last_row_end = last_row_init + Nx-1;
        number_faces_one_row = 2*Nx+1; // for each cell I count only the low and sx faces, respectevely 0-1 2-3 4-5 6-7 8 + the last on the right boundary
        
        Matrix<T, Dynamic, 1> RHS_vandermonde = Matrix<T, Dynamic, 1>::Zero( local_dim );
       
        Global_Mass = SparseMatrix<T>( ndof_FE, ndof_FE ); //(b_i,b_j)_ij , b_i Lagrange basis fx
        RHS = Matrix<T, Dynamic, 1>::Zero( ndof_FE ); // Known term (f,b_i)_i , b_i Lagrange basis fx
        
        Global_c_term_x = SparseMatrix<T>( ndof_FE, ndof_FE ); //(b_i,b_j)_ij , b_i Lagrange basis fx
        Global_c_term_y = SparseMatrix<T>( ndof_FE, ndof_FE ); //(b_i,b_j)_ij , b_i Lagrange basis fx
         
        Global_Mass_Lumped = Eigen::Matrix<T, Dynamic, 1>::Zero( ndof_FE , 1 );
         
        local_vandermonde = Matrix<T, Dynamic, Dynamic>::Zero ( local_dim , local_dim ) ;
        
        
        // Saving the projection in HHO discontinuous format (MATRIX NOTATION)
        sol_HHO = Eigen::Matrix<T, Dynamic, Dynamic>::Zero( local_dim, n_cls );

        
        //Eigen::Matrix<T, Dynamic, Dynamic> sol_HHO_vandermonde =  Eigen::Matrix<T, Dynamic, Dynamic>::Zero( local_dim, n_cls ) ;
        
        
        CompleteOrthogonalDecomposition<Matrix<T, Dynamic, Dynamic > > cod;
        
        std::cout<<yellow<<bold<<"-----> IN 'L2projected_level_set_high_order': Vandermonde interpolation of the level set with BERNSTEIN basis. -----"<<reset<<std::endl;
      
        for( const auto& cl : msh.cells )
        {
            size_t cell_offset = offset(msh, cl) ;
            size_t cell_LHS_offset = cell_offset * degree_FEM ;
            std::vector<std::pair<size_t , bool >> asm_map(local_dim) ;
            
            if(cell_offset == 0)
            {
                auto pts = equidistriduted_nodes_ordered_bis<T,Mesh>( msh, cl, degree_FEM);
                cell_basis_Bernstein<Mesh,T> cb(msh, cl, degree_FEM);
                for (size_t i = 0; i < local_dim; i++)
                {
                    // Local vandermonde matrix
                    local_vandermonde.block(i,0,1,local_dim) = (cb.eval_basis(pts[i])).transpose() ;
                }
                
                cod.compute(local_vandermonde);
            }
            
            // Adding assembly_index for cl
            for (size_t i = 0; i < local_dim; i++)
            {
                asm_map[i] =   std::make_pair( connectivity_matrix[cell_offset][i].first, connectivity_matrix[cell_offset][i].second ) ;
            }
              
       
       
            auto local_mass = make_bernstein_local_mass_matrix( msh, cl , degree_FEM );
            
            //auto local_RHS = make_bernstein_local_RHS( msh , cl , degree_FEM , level_set );
            // Local c_ij = b_i nabla(b_j) -> USEFUL FOR TRANSPORT PROBLEM
            auto local_cij = make_bernstein_local_cij_matrix (msh, cl, degree_FEM);
                
            auto local_mass_lumped = make_bernstein_local_mass_matrix_lumped( msh , cl , degree_FEM ) ;
           
            
            // Costruction of the coefficients of the Bernstein basis
            auto qps = equidistriduted_nodes_ordered_bis<T,Mesh>(msh, cl, degree_FEM);
            
            
            // Assembling triplets for global problem
            for (size_t i = 0; i < local_dim; i++)
            {
                for (size_t j = 0; j < local_dim; j++)
                {
                    triplets.push_back( Triplet<T>(asm_map[i].first, asm_map[j].first, local_mass(i,j) ) );
                
                    triplets_c_term_x.push_back( Triplet<T>(asm_map[i].first, asm_map[j].first, local_cij.first(i,j) ) );
                    triplets_c_term_y.push_back( Triplet<T>(asm_map[i].first, asm_map[j].first, local_cij.second(i,j) ) );
                        
                }
                Global_Mass_Lumped(asm_map[i].first) += local_mass_lumped(i);
                //RHS(asm_map[i].first) += local_RHS(i) ;
                RHS_vandermonde(i) = level_set( qps[i]) ;
            }
            
            
            sol_HHO.col(cell_offset) = cod.solve(RHS_vandermonde) ;
                
        } // end of cl loop
            
        // Finalisation global assembling
        Global_Mass.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();
            
            
        // Finalisation global assembling
        Global_c_term_x.setFromTriplets( triplets_c_term_x.begin(), triplets_c_term_x.end() );
        triplets_c_term_x.clear();
        Global_c_term_y.setFromTriplets( triplets_c_term_y.begin(), triplets_c_term_y.end() );
        triplets_c_term_y.clear();
            
        
        //std::cout<<"local_vandermonde"<<'\n'<<local_vandermonde<<std::endl;
            
        
        // CALCULATION OF THE SIZE + PLOTTING
        /*
        size_t size_supp_nodes = 0;
        std::cout<<"Supporting nodes IN L2:"<<std::endl;
        size_t jjjj = 0;
        for (auto& i: S_i) {
            size_supp_nodes+=i.size();
            std::cout <<"Node "<<jjjj<<":";
            for (auto it=i.begin(); it != i.end(); ++it)
                std::cout << ' ' << *it;
                // std::cout<<ii;
            std::cout<<'\n';
            jjjj++;
        }
        std::cout<<std::endl;
        std::cout<<"Supporting nodes size:"<<size_supp_nodes<<std::endl;
        */
        
            
        // Resolution L2-projection over a Bernstein basis
        sol_FEM = Matrix<T, Dynamic, 1>::Zero(RHS.rows());
            
        //Matrix<T, Dynamic, 1> sol_FEM_vandermonde  = Matrix<T, Dynamic, 1>::Zero(RHS.rows()); ;
        
        std::cout<<"sol_FEM size "<<sol_FEM.size()<<std::endl;
        std::cout<<"local_dim size "<<local_dim<<std::endl;
        std::cout<<"n_cls "<<n_cls<<std::endl;
            
        // Global solution saved as discontinuous HHO approach
        // Also saved min & max coefficients + position in HHO notation
       
        /*
        ConjugateGradient<SparseMatrix<T> > solver_global_mass;
        //SparseLU<SparseMatrix<T>, COLAMDOrdering<int> >solver_global_mass;
        //Notice: for this step the numerical values of A are not used
        //solver_global_mass.analyzePattern(Global_Mass);
        //solver_global_mass.factorize(Global_Mass);
         
        solver_global_mass.compute(Global_Mass); // SAVE INVERSE OF GLOBAL MASS
        if(solver_global_mass.info()!=Success) {
            std::cout<<"FAILED SOLVER 0"<<std::endl;
            return;
        }
       
        sol_FEM = solver_global_mass.solve(RHS);
       */
        for(size_t counter_bis = 0 ; counter_bis < n_cls ;counter_bis++)
        {
            for (size_t i = 0; i < local_dim; i++)
            {
                size_t asm_map =  connectivity_matrix[counter_bis][i].first ;
                //sol_HHO(i,counter_bis) = sol_FEM( asm_map ) ;
                sol_FEM( asm_map ) = sol_HHO(i,counter_bis) ;
            }
        }
       
        //std::cout<<"CHECK L2 proj vs Vandermonde interpolation for phi^0"<<'\n'<<sol_FEM-sol_FEM_vandermonde<<std::endl;
        
        
        //sol_FEM = sol_FEM_vandermonde ;
        //sol_HHO = sol_HHO_vandermonde ;
        //std::cout<<"CHECK L2 proj vs Vandermonde interpolation for phi^0"<<'\n'<<sol_FEM-sol_FEM_vandermonde<<std::endl;
        // Set of maximum and minimum
        set_max_min();
        
        vertices = Eigen::Matrix<T, Dynamic, 1>::Zero( n_vertices , 1 );
      
        for( size_t i_global = 0; i_global < n_cls; i_global++)
        {
            size_t i_vertex = i_global+floor(i_global/Nx);
            vertices(i_vertex) = sol_HHO(0,i_global) ;
            vertices(i_vertex+1) = sol_HHO(1,i_global) ;
            vertices(i_vertex+Nx+2) = sol_HHO(2,i_global) ;
            vertices(i_vertex+Nx+1) = sol_HHO(3,i_global) ;
        }
    
       
    }
            
            

   
    L2projected_level_set_high_order()=default;
    
    
    void
    coefficients_mapping( )
    {

        
        Matrix<T, Dynamic, 1> mapped_phi = Matrix<T, Dynamic, 1>::Zero( ndof_FE );
        
        CompleteOrthogonalDecomposition<Matrix<T, Dynamic, Dynamic > > vandermonde_interpolant( local_vandermonde );
        Matrix<T, Dynamic, 1> RHS_vandermonde = Matrix<T, Dynamic, 1>::Zero( local_dim );
        Matrix<T, Dynamic, 1> sol_vandermonde = Matrix<T, Dynamic, 1>::Zero( local_dim );
        
        for( const auto& cl : msh.cells )
        {
            size_t cell_offset = offset(msh, cl) ;
            cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
            auto nodes = equidistriduted_nodes_ordered_bis<T,Mesh>(msh , cl , degree_FEM);
            size_t ct = 0;
            for (auto& pt : nodes )
            {
                RHS_vandermonde(ct) =  ( (*this)(pt , msh , cl ) - phi_min )/( phi_max - phi_min );
                ct++;
            }
            
            sol_vandermonde = vandermonde_interpolant.solve(RHS_vandermonde);
         
            for (size_t i = 0; i < local_dim; i++)
            {
                auto asm_map = connectivity_matrix[cell_offset][i].first ;
                mapped_phi(asm_map) = sol_vandermonde(i) ;
            }
             
        } // end of cl loop
       
        sol_FEM = mapped_phi ;
        converting_into_HHO_formulation( sol_FEM );
        
    }
    
    void
    coefficients_inverse_mapping( )
    {
        
        Matrix<T, Dynamic, 1> mapped_phi = Matrix<T, Dynamic, 1>::Zero( ndof_FE );
        
        CompleteOrthogonalDecomposition<Matrix<T, Dynamic, Dynamic > > vandermonde_interpolant( local_vandermonde );
        Matrix<T, Dynamic, 1> RHS_vandermonde = Matrix<T, Dynamic, 1>::Zero( local_dim );
        Matrix<T, Dynamic, 1> sol_vandermonde = Matrix<T, Dynamic, 1>::Zero( local_dim );
        
        for( const auto& cl : msh.cells )
        {
            size_t cell_offset = offset(msh, cl) ;
            
            cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
            auto nodes = equidistriduted_nodes_ordered_bis<T,Mesh>(msh , cl , degree_FEM);
            size_t ct = 0;
            for (auto& pt : nodes )
            {
                RHS_vandermonde(ct) = phi_min + (*this)(pt , msh , cl )*( phi_max - phi_min );
                ct++;
            }
            sol_vandermonde = vandermonde_interpolant.solve(RHS_vandermonde);
            
            for (size_t i = 0; i < local_dim; i++)
            {
                auto asm_map =  connectivity_matrix[cell_offset][i].first;
                mapped_phi(asm_map) = sol_vandermonde(i) ;
            }
             
            
        } // end of cl loop
        
        sol_FEM = mapped_phi ;
        converting_into_HHO_formulation( sol_FEM );
           
    }
        
        
        
        
        
    void  set_discrete_points( Eigen::Matrix<T, Dynamic, Dynamic>& values_new)
    {
        sol_HHO = values_new ;
        std::cout<<" Using set_discrete_points check that sol_FEM already uploaded!"<<std::endl;
        
    }
        
        
    void converting_into_HHO_formulation( const Eigen::Matrix<T, Dynamic, 1>& values_new )
    {
        // SAVE BOTH SOL_HHO AND VERTICES
        for(size_t counter_bis = 0 ; counter_bis < n_cls ;counter_bis++)
        {
            for (size_t i = 0; i < local_dim; i++){
                size_t asm_map =  connectivity_matrix[counter_bis][i].first ;
                sol_HHO(i,counter_bis) = values_new( asm_map );
                
            }
            size_t i_vertex = counter_bis+floor(counter_bis/Nx);
            vertices(i_vertex) = sol_HHO(0,counter_bis) ;
            vertices(i_vertex+1) = sol_HHO(1,counter_bis) ;
            vertices(i_vertex+Nx+2) = sol_HHO(2,counter_bis) ;
            vertices(i_vertex+Nx+1) = sol_HHO(3,counter_bis) ;
                
        }
        std::cout<<" Using converting_into_HHO_formulation. TO BE CHECKED that sol_FEM already uploaded!"<<std::endl;
        
    }
        
    void converting_into_FE_formulation( const Eigen::Matrix<T, Dynamic, Dynamic>& values_new )
    {
          
        for(size_t counter_bis = 0 ; counter_bis < n_cls ;counter_bis++)
        {
            for (size_t i = 0; i < local_dim; i++){
                size_t asm_map =  connectivity_matrix[counter_bis][i].first ;
                sol_FEM( asm_map ) = values_new(i,counter_bis) ;
            }
            size_t i_vertex = counter_bis+floor(counter_bis/Nx);
            vertices(i_vertex) = values_new(0,counter_bis) ;
            vertices(i_vertex+1) = values_new(1,counter_bis) ;
            vertices(i_vertex+Nx+2) = values_new(2,counter_bis) ;
            vertices(i_vertex+Nx+1) = values_new(3,counter_bis) ;

        }
        std::cout<<" Using converting_into_FE_formulation. TO BE CHECKED that sol_HHO already uploaded!"<<std::endl;
       
    }
    
    
    
           
    void set_max_min()
    {
        
        T ret0 = -10.0;
        T ret1 = 10.0;
        size_t counter_ret0 = 0;
        
        for(auto& cl:msh.cells)
        {
            cell_basis_Bernstein<Mesh,T> cb(msh, cl, degree_FEM);
            auto nodes = equidistriduted_nodes_ordered_bis<T,Mesh>(msh, cl , degree_FEM);
            
            for(auto& nd : nodes ){
                auto phi_tmp = cb.eval_basis(nd);
                auto values_cell = (sol_HHO.block(0,counter_ret0,local_dim,1)).col(0);
                auto new_ret = values_cell.dot(phi_tmp) ;
                ret0 = std::max( new_ret , ret0 ) ;
                ret1 = std::min( new_ret , ret1);
            }
            counter_ret0++;
        }
        
        phi_max = ret0;
        phi_min = ret1;
            
    }
        
       
    // OK FINE, IT WORKS ALSO FOR AGGLOMERATED MESHES
    // It should work also for Bernstein Basis
    T operator()( const typename Mesh::node_type& node ) const
    {
        return vertices(node.ptid);
            
    }
        
     
    // OK FINE, IT WORKS ALSO FOR AGGLOMERATED MESHES -> BUT SLOW
    T operator()(const point<T,2>& pt) const
    {
        //std::cout<<"I AM IN OPERATOR() SLOW !!!!"<<std::endl;
        size_t counter=0;
            
        // It looks for in what cell the point is
        for( const auto& cl:msh.cells)
        {
            if( pt_in_cell<T,Mesh>(msh,pt,cl) )
            {
                //cell_basis_Lagrangian_ordered<Mesh,T> cb(msh, cl, degree_FEM);
                cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
                    
                auto values_cell = (sol_HHO.block(0,counter,local_dim,1)).col(0);
               
                return values_cell.dot( cb.eval_basis(pt) );
                
            }
            counter+=1;
        }
        std::cout<<"IF HERE, THERE IS A PROBLEM IN projected_level_set::operator()!!!"<<std::endl;
        return 1e10; //to check if doesn't enter in the loop
    }

        
    // IT WORKS ONLY FOR NOT-AGGLOMERATED MESHES --> FAST
    T operator()( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl ) const
    {
            
        size_t counter = offset(msh,cl) ;
        cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
        //cell_basis_Lagrangian_ordered<Mesh,T> cb(msh, cl, degree_FEM);
        auto values_cell = (sol_HHO.block(0,counter,local_dim,1)).col(0);
        T tmp = values_cell.dot( cb.eval_basis(pt) );
           
        return tmp;
                
    }
        
    // IT WORKS ONLY FOR NOT-AGGLOMERATED MESHES --> FAST
    T operator()( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::face_type& fc ) const
    {
        auto counter_face = offset(msh,fc);
        size_t counter_cell;
        // ATTENTION, ALL THIS WORKS IN STRUCTURED QUADRANGULAR MESHES
               
        // Check if I am in the last row, upper faces (ordered differently)
        if(counter_face>=last_row_init && counter_face<=last_row_end)
        {
            counter_cell = (Ny-1)*Nx + counter_face%(last_row_init);
        }
        else
        {
            // Find in what row the face is
            auto  num_cell_row = floor(counter_face/(number_faces_one_row));
            if ( counter_face!= ( (2*Nx)*(num_cell_row+1)+num_cell_row ) )
            {
                // Faces not on the right boudary, to know in what cell are, it is sufficient to analyse the low and left face of each quadrangular cell.
                counter_cell = floor( (counter_face-num_cell_row)/2.0 );
            }
            else
            {
                // Face on the right boudary,
                counter_cell = ( num_cell_row+1 )*Nx -1;
            }
            
        }
        //std::cout<<"Face->Cell number "<<counter_cell<<std::endl;
        auto cl = msh.cells.at(counter_cell);
        cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
        //cell_basis_Lagrangian_ordered<Mesh,T> cb(msh, cl, degree_FEM);
        auto values_cell = (sol_HHO.block(0,counter_cell,local_dim,1)).col(0);
        T tmp = values_cell.dot( cb.eval_basis(pt) );
        return tmp;

               
    }
        
        
    // IT WORKS FOR ALL THE MESHES --> SLOW
    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;

        ret = gradient(pt);
        return ret/ret.norm();
    }
        
         
    // IT WORKS FOR ALL THE MESHES --> SLOW
    Eigen::Matrix<T,2,1> gradient(const point<T,2>& pt) const
    {
        size_t counter=0;
        //std::cout<<"I AM IN GRADIENT SLOW !!!!"<<std::endl;
        Eigen::Matrix<T,2,1> ret = Matrix<T, 2, 1>::Zero(2, 1);
        for( const auto& cl:msh.cells)
        {
            if(pt_in_cell<T,Mesh>(msh,pt,cl))
            {
                cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
                //cell_basis_Lagrangian_ordered<Mesh,T> cb(msh, cl, degree_FEM);
                    
                auto values_cell = sol_HHO.col(counter);
                auto grad_eval =  cb.eval_gradients(pt);
                ret(0) = values_cell.dot( grad_eval.col(0) );
                // std::cout<<"Value of derivative new along x"<<ret(0)<<std::endl;
                ret(1) = values_cell.dot( grad_eval.col(1) );

                return ret;
            }
            counter+=1;

        }
        std::cout<<"Se compare questo problema in gradient()"<<std::endl;
        ret(0)+=1e10;
        ret(1)+=1e10;
        return ret; //to check if doesn't enter in the loop

    }

         
        // IT WORKS FOR NOT-AGGLOMERATED MESHES --> FAST
    Eigen::Matrix<T,2,1> gradient( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl ) const
    {
        
        // MATRIX NOTATION
        size_t counter = offset(msh,cl);
        Eigen::Matrix<T,2,1> ret = Matrix<T, 2, 1>::Zero(2, 1);
        cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
        //cell_basis_Lagrangian_ordered<Mesh,T> cb(msh, cl, degree_FEM);
        auto values_cell = sol_HHO.col(counter);
        auto grad_eval =  cb.eval_gradients(pt);
        ret(0) = values_cell.dot( grad_eval.col(0) );
        // std::cout<<"Value of derivative new along x"<<ret(0)<<std::endl;
        ret(1) = values_cell.dot( grad_eval.col(1) );
        //values_cell.dot( grad_eval.col(1) );
        // std::cout<<"Value of derivative new along y"<<ret(1)<<std::endl;
        return ret;
            
    }

     
    // IT WORKS FOR NOT-AGGLOMERATED MESHES --> FAST
    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt, const Mesh& msh, const typename Mesh::cell_type& cl) const
    {
        Eigen::Matrix<T,2,1> ret;
        ret = gradient(pt,msh,cl);
        return ret/ret.norm();
            
    }
    
    T divergence( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl ) const
    {
            
        size_t counter = offset(msh,cl) ;
        cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
        //cell_basis_Lagrangian_ordered<Mesh,T> cb(msh, cl, degree_FEM);
        auto values_cell = (sol_HHO.block(0,counter,local_dim,1)).col(0);
        
        T grad_norm = (this->gradient( pt , msh , cl )).norm() ;
    
        //std::cout<<"grad norm is "<<grad_norm<<std::endl;
        //std::cout<<"values_cell is "<<'\n'<<values_cell<<std::endl;
        //std::cout<<"cb.eval_divergence(pt) is "<<'\n'<<cb.eval_divergence(pt)<<std::endl;
        //std::cout<<"( values_cell.dot(cb.eval_divergence(pt)) ) is "<<( values_cell.dot(cb.eval_divergence(pt)) )<<std::endl;
           
        return -( values_cell.dot(cb.eval_divergence(pt)) ) / (2 * grad_norm) ;
                
    }

    
    void smooth_cut_off( T C , T x_centre , T y_centre , T radius )
    {
        
        Matrix<T, Dynamic,1> local_RHS = Matrix<T, Dynamic, 1>::Zero( local_dim , 1);
        //Matrix<T, Dynamic, 1> RHS = Matrix<T, Dynamic, 1>::Zero( ndof_FE , 1 );
        Matrix<T, Dynamic,1> sol_loc = Matrix<T, Dynamic, 1>::Zero( local_dim , 1);
        CompleteOrthogonalDecomposition<Matrix<T, Dynamic, Dynamic > > vandermonde_interpolant( local_vandermonde );
        
        T hx = params.hx();
        T hy = params.hy();
        T pos_r0 = 0.45; //std::min(x_centre , 1 - x_centre );
        T dist = pos_r0 - radius + 2.0*0.0625;
        T r0 = radius + dist/2.0;
        T delta = r0/8.0;
        std::cout<<"radius = "<<radius<<" , r0 = "<<r0<<" , delta = "<<delta<<" , hx = hy = "<<hx<<std::endl;
        std::cout<<"value in alfa in r_int = "<<(radius-r0)/delta<<std::endl;
        std::cout<<"value in alfa in R = "<<(pos_r0-r0)/delta<<std::endl;
        
        // Lambda function to define smooth function
        auto alfa = [=](const point<T,2>& pt)
        { // sol
            return (1 - tanh( (sqrt( pow((pt.x()-x_centre),2) + pow((pt.y()-y_centre),2) ) - r0 ) / delta ))/2;};
        
        
        
            
        for( const auto& cl : msh.cells )
        {
            size_t cell_offset = offset(msh, cl) ;
             
            auto qps = equidistriduted_nodes_ordered_bis<T,Mesh>(msh, cl, degree_FEM);
            for(size_t i = 0 ; i<local_dim ; i++)
            {
                auto pt = qps[i];
                auto asm_map =  connectivity_matrix[cell_offset][i].first;
                local_RHS(i) =  ( 1.0 - alfa(pt)  ) * C + alfa(pt) * (*this)(pt,msh,cl) ;
            }
            
            sol_loc = vandermonde_interpolant.solve(local_RHS); // SAVE Vandermonde interpolation
            
            for (size_t i = 0; i < local_dim; i++)
            {
                auto asm_map =  connectivity_matrix[cell_offset][i].first;
                sol_FEM(asm_map) = sol_loc(i) ;
            }
            
            
            
        } // end of cl loop
             
            
            
        converting_into_HHO_formulation(sol_FEM);

        set_max_min();
        
        postprocess_output<double> postoutput00;
        typedef typename Mesh::point_type       point_type;
        point<double,2> node;
        auto alfa_values = std::make_shared< gnuplot_output_object<double> >("alfa.dat");
        //auto interface_pos = std::make_shared< gnuplot_output_object<double> >("interface_alfa.dat");
        for(auto& pt:msh.points )
        {
            alfa_values->add_data(pt,alfa(pt));
        }
        postoutput00.add_object(alfa_values);
        postoutput00.write();
    }
    
    
    void cut_off( T d )
    {
        
        Matrix<T, Dynamic,1> local_RHS = Matrix<T, Dynamic, 1>::Zero( local_dim , 1);
        Matrix<T, Dynamic,1> sol_loc = Matrix<T, Dynamic, 1>::Zero( local_dim , 1);
        
        CompleteOrthogonalDecomposition<Matrix<T, Dynamic, Dynamic > > vandermonde_interpolant( local_vandermonde );
        
        
        // Known term (f,b_i)_i , b_i Bernstein basis fx
        for( const auto& cl : msh.cells )
        {
             size_t cell_offset = offset(msh, cl) ;
              
             auto qps = equidistriduted_nodes_ordered_bis<T,Mesh>(msh, cl, degree_FEM);
             for(size_t i = 0 ; i<local_dim ; i++)
             {
                 auto pt = qps[i];
                 auto asm_map =  connectivity_matrix[cell_offset][i].first;
                 
                 if( (*this)(pt,msh,cl)>= d )
                     local_RHS(i) =  d ;
                 else if( (*this)(pt,msh,cl)<= -d )
                     local_RHS(i) =  -d ;
                 else
                     local_RHS(i) =  (*this)(pt,msh,cl) ;
        
             }
            
             sol_loc = vandermonde_interpolant.solve(local_RHS); // SAVE Vandermonde interpolation
            
            for (size_t i = 0; i < local_dim; i++)
             {
                 auto asm_map =  connectivity_matrix[cell_offset][i].first;
                 sol_FEM(asm_map) = sol_loc(i) ;
             }
                  
         } // end of cl loop
              
       
        
        
        converting_into_HHO_formulation(sol_FEM);

        set_max_min();
        
        
     }
    
        
        
};


template< typename T , typename Mesh ,typename Level_Set,typename Fonction,typename FiniteSpace >
struct LS_cell_L2proj_high_order: public L2projected_level_set_high_order< Mesh,Fonction,FiniteSpace , T >
{
    
    typedef typename Mesh::cell_type       cell_type;
    cell_type agglo_LS_cl;
    std::vector<cell_type> subcells;
    Mesh agglo_msh;
    Level_Set level_set;
    //LS_cell(const Level_Set & level_set, const Mesh & msh, const typename Mesh::cell_type& cl)
   // : agglo_cl(cl), agglo_msh(msh), level_set(level_set){}
    // I don't know if I have to define a copyconstructor for level_set.. TO BE CHECKED!
    LS_cell_L2proj_high_order(const Level_Set & level_set_, const Mesh & msh)
    : agglo_msh(msh), level_set(level_set_){}
    //LS_cell(const Level_Set & level_set )
    //: level_set(level_set){}

    //LS_cell()=default;
    
    T operator()(const point<T,2>& pt) const
       {
           if (subcells.size()<1)
           {
               assert(agglo_LS_cl.user_data.offset_subcells.size()==1);
               auto offset_old = agglo_LS_cl.user_data.offset_subcells[0];
               auto cl_old = level_set.msh.cells[offset_old];
               return level_set( pt , level_set.msh , cl_old );
           }
       
           else
           {
               auto offset = pt_in_subcell(level_set.msh,pt,agglo_LS_cl);
               auto subcl = level_set.msh.cells[offset];
               return level_set( pt , level_set.msh , subcl );
           }
               
       }
    
    
    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt) const
    {
        if (subcells.size()<1)
            {
                assert(agglo_LS_cl.user_data.offset_subcells.size()==1);
                auto offset_old = agglo_LS_cl.user_data.offset_subcells[0];
                auto cl_old = level_set.msh.cells[offset_old];
                return level_set.normal( pt , level_set.msh , cl_old );
            }
        
            else
            {
                auto offset = pt_in_subcell(level_set.msh,pt,agglo_LS_cl);
                auto subcl = level_set.msh.cells[offset];
                return level_set.normal( pt , level_set.msh , subcl );
            }
    }
    
    Eigen::Matrix<T,2,1> gradient(const point<T,2>& pt) const
    {
        if (subcells.size()<1)
            {
                std::cout<<"I m here 2.5 and size is "<<agglo_LS_cl.user_data.offset_subcells.size()<<std::endl;
                assert(agglo_LS_cl.user_data.offset_subcells.size()==1);
                auto offset_old = agglo_LS_cl.user_data.offset_subcells[0];
                auto cl_old = level_set.msh.cells[offset_old];
                return level_set.gradient( pt , level_set.msh , cl_old );
            }
        
            else
            {
                auto offset = pt_in_subcell(level_set.msh,pt,agglo_LS_cl);
                auto subcl = level_set.msh.cells[offset];
                return level_set.gradient( pt , level_set.msh , subcl );
            }
        
    }
    
    
    T divergence( const point<T,2>& pt )const
    {
        if (subcells.size()<1)
            {
                //std::cout<<"I m here 2.5 and size is "<<agglo_LS_cl.user_data.offset_subcells.size()<<std::endl;
                assert(agglo_LS_cl.user_data.offset_subcells.size()==1);
                auto offset_old = agglo_LS_cl.user_data.offset_subcells[0];
                auto cl_old = level_set.msh.cells[offset_old];
                return level_set.divergence( pt , level_set.msh , cl_old );
            }
        
            else
            {
                auto offset = pt_in_subcell(level_set.msh,pt,agglo_LS_cl);
                auto subcl = level_set.msh.cells[offset];
                return level_set.divergence( pt , level_set.msh , subcl );
            }
    }
  
    void cell_assignment(const cell_type& cl)
    {
        subcells.clear();
        agglo_LS_cl = cl;
        if (agglo_LS_cl.user_data.offset_subcells.size()>1)
        {

            for (auto& offset_subcells:agglo_LS_cl.user_data.offset_subcells)
            {
                subcells.push_back( level_set.msh.cells[offset_subcells] );
            
            }
        }
            
    }
    
    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt , const cell_type& cl)
    {
        //agglo_LS_cl = cl;
        cell_assignment(cl);
        return normal( pt );
    }
    
    T operator()(const point<T,2>& pt,const cell_type& cl )
    {
        //agglo_LS_cl = cl;
        cell_assignment(cl);
        return operator()( pt );
    }
    
    
  //  void get_cell()
  //  {
  //      auto crr_cl = download_current_cell();
  //      std::cout<<"IN get_cell.. size cell "<<crr_cl.user_data.offset_subcells.size()<<std::endl;
  //      this->cell_assignment(crr_cl);
 //   }
    
};




template<typename T, typename Mesh ,typename Fonction >
struct projected_level_set: public level_set<T>
{
    std::vector< T > values;
    Eigen::Matrix<T, Dynamic, Dynamic> values_bis; // MATRIX NOTATION
    //Eigen::Matrix<T, Dynamic, 1> values_bis1; // VECTOR NOTATION
    Eigen::Matrix<T, Dynamic, 1> vertices; // saving level_set on vertices mesh
    std::vector< size_t > boundary_cells;
    // In the case in which I use FEM with Qk , k>1 I need another vector for all the values in each node
    
    size_t degree_FEM;
    size_t number_elements;
    Mesh msh;
    size_t      Nx, Ny ;
    mesh_init_params<T> params;
    size_t  last_row_init, last_row_end, number_faces_one_row;
  //  size_t counter_cell , counter_face, num_cell_row;
    
    T phi_max , phi_min ;
    T cut_level;
    
    
    projected_level_set(const Fonction & level_set, const Mesh & msh, size_t degree_k , const mesh_init_params<T>& params)
        : number_elements((degree_k+1)*(degree_k+1)), msh(msh),degree_FEM(degree_k),Nx(params.Nx),Ny(params.Ny), params(params)
    {
        
        last_row_init = Ny*(2*Nx+1); // There are 2 faces for each row of cells + Ny
        last_row_end = last_row_init + Nx-1;
        number_faces_one_row = 2*Nx+1; // for each cell I count only the low and sx faces, respectevely 0-1 2-3 4-5 6-7 8 + the last on the right boundary
        vertices = Eigen::Matrix<T, Dynamic, 1>::Zero( ((Nx+1)*(Ny+1)), 1 );
       
    //#ifdef NODES
        // MATRIX NOTATION
        values_bis= Eigen::Matrix<T, Dynamic, Dynamic>::Zero(number_elements, msh.cells.size());
        // VECTOR NOTATION
        //values_bis1= Eigen::Matrix<T, Dynamic, 1>::Zero(number_elements*msh.cells.size(), 1 );
        
        //std::cout<<"Number of cells "<<msh.cells.size()<<std::endl;
        
        // MATRIX NOTATION
        // std::cout<<"Dimension of the basis "<<values_bis.size()<<std::endl;
        // VECTOR NOTATION
        // std::cout<<"Dimension of the basis "<<values_bis1.size()<<std::endl;
        size_t i_global = 0 , i_local=0 , i_vertex=0;
        for(auto& cl:msh.cells)
        {
            /*
            bool boundary_bool = FALSE;
            for (auto& fc:faces(msh,cl)) {
                if (boundary_bool)
                    break;
    
                if(fc.is_boundary && !boundary_bool){
                    boundary_cells.push_back(offset(msh,cl));
                    boundary_bool = TRUE;
                }
            }
            */
            
            auto qps = equidistriduted_nodes<T,Mesh>(msh, cl, degree_FEM);
            i_local = 0;
            for ( const auto & qp : qps)
            {
                
                values.push_back( level_set(qp) ); // I DONT KNOW IF IT IS USEFUL
                
               // if (boundary_bool) {
               //      values_bis(i_local,i_global) = 1 ;  // MATRIX NOTATION
               // }
                values_bis(i_local,i_global) = level_set(qp) ;  // MATRIX NOTATION
                //values_bis1(i_local+i_global) = level_set(qp) ; // VECTOR NOTATION
                i_vertex = i_global+floor(i_global/Nx);
                if( i_local==0 )
                    vertices(i_vertex) = level_set(qp) ;
                
                if( i_local==1 )
                    vertices(i_vertex+1) = level_set(qp) ;
                
                if( i_local==(degree_FEM+2) )
                    vertices(i_vertex+Nx+2) = level_set(qp) ;
                
                if( i_local==(degree_FEM+1) )
                    vertices(i_vertex+Nx+1) = level_set(qp) ;
                i_local++;
            }
             i_global++;  // MATRIX NOTATION
           //  i_global+=number_elements;       // VECTOR NOTATION
        }
    //#endif
        
        
    }
    
    projected_level_set()=default;
    
    
    void  set_discrete_points( Eigen::Matrix<T, Dynamic, Dynamic>& values_new)
    {
        values_bis=values_new;
    }
    
    
    void converting_into_HHO_formulation( const Eigen::Matrix<T, Dynamic, 1>& vertices )
    {
        assert(degree_FEM == 1);
        for (size_t j = 0; j < values_bis.cols() ; j++)
        {

            for (size_t i = 0; i < values_bis.rows() ; i++)
            {
                auto i_vertex = j + floor(j/Nx);
                if ( i == 0 )
                    values_bis(i,j) = vertices(i_vertex);
                
                if( i ==1 )
                    values_bis(i,j) = vertices(i_vertex+1) ;
                
                if( i ==(degree_FEM+2) )
                    values_bis(i,j) = vertices(i_vertex+Nx+2) ;
                
                if( i ==(degree_FEM+1) )
                    values_bis(i,j) = vertices(i_vertex+Nx+1) ;
            }
        }
                 
    }
    
    void converting_into_FE_formulation( const Eigen::Matrix<T, Dynamic, Dynamic>& values )
    {
        assert(degree_FEM == 1);
        for (size_t j = 0; j < values_bis.cols() ; j++)
        {

            for (size_t i = 0; i < values_bis.rows() ; i++)
            {
                auto i_vertex = j + floor(j/Nx);
                if ( i == 0 )
                    vertices(i_vertex) = values(i,j) ;
                
                if( i ==1 )
                    vertices(i_vertex+1) = values(i,j) ;
                
                if( i ==(degree_FEM+2) )
                    vertices(i_vertex+Nx+2) = values(i,j) ;
                
                if( i ==(degree_FEM+1) )
                    vertices(i_vertex+Nx+1) = values(i,j) ;
            }
        }
                 
    }
       
    
    
    // OK FINE, IT WORKS ALSO FOR AGGLOMERATED MESHES
    T operator()( const typename Mesh::node_type& node ) const
    {
        // Optimised to check the value of the level set only in the vertices
        // std::cout<<"Value in vertices "<<vertices(node.ptid)<<", at position "<<node.ptid<<std::endl;
        return vertices(node.ptid);
        
    }
    
 
    // OK FINE, IT WORKS ALSO FOR AGGLOMERATED MESHES
    T operator()(const point<T,2>& pt) const
    {
        //std::cout<<"I AM IN OPERATOR() SLOW !!!!"<<std::endl;
        size_t counter=0;
        
        // It looks for in what cell the point is
        for( const auto& cl:msh.cells)
        {
            //std::cout<<"pt_in_cell operator in slow levelset"<<std::endl;
            if( pt_in_cell<T,Mesh>(msh,pt,cl) )
            {
                //cell_basis_Qk<Mesh,T> cb(msh, cl, degree_FEM);
                cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
               
                auto values_cell = (values_bis.block(0,counter,number_elements,1)).col(0);
                return values_cell.dot( cb.eval_basis(pt) );
                
               // T tmp=0;
               // for(auto i = 0; i<number_elements ; i++)
               // {
               //     tmp += (values.at(i+counter))*(cb.eval_basis(pt))[i];
               // }
               // return tmp;
            }
            //counter+=number_elements; // OLD VERSION OF THE CODE
            counter+=1;
        }
        std::cout<<"IF HERE, THERE IS A PROBLEM IN projected_level_set::operator()!!!"<<std::endl;
        return 1e10; //to check if doesn't enter in the loop
    }

    // IT WORKS ONLY FOR NOT-AGGLOMERATED MESHES
    T operator()( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl ) const
    {
        
        size_t counter = offset(msh,cl) ;
        // Checking if cell is agglomerated or not
        /*
        size_t cell_points = (points(msh,cl)).size();
        size_t cells_number = msh.cells.size();
        bool agglomeration=FALSE;
        if ( cells_number < (Nx*Ny-1) ) {
            agglomeration = TRUE;
        }
        if ( cells_number > (Nx*Ny-1) ) {
                  // throw std::logic_error("shouldn't have arrived here...");
               }
        */
        
        //if(cell_points<5)
        //{
            // MATRIX NOTATION
         //   if (!agglomeration) {
         //       counter = offset(msh,cl);
         //    }
         //   else{
         //       counter = 1;//cells_counter(msh_origin,msh,cl);
         //   }
        
            //std::cout<<"Value of offset "<<counter<<std::endl;
            //cell_basis_Qk<Mesh,T> cb(msh, cl, degree_FEM);
        cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
        auto values_cell = (values_bis.block(0,counter,number_elements,1)).col(0);
        T tmp = values_cell.dot( cb.eval_basis(pt) );
        return tmp;
            
            // VECTOR NOTATION
            
            //size_t counter = offset(msh,cl)*number_elements;
            //cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
            //auto values_cell = (values_bis1.segment(counter,number_elements));
            //T tmp = values_cell.dot( cb.eval_basis(pt) );
            //return tmp;
       // }
        /*
        else
        {
            std::vector<size_t> indices = subcell_finder<T,Mesh>(msh, pt , cl , params);
            cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM,indices);
            size_t counter = 1; // DA TROVAREEE!!!!!!!!!!!!!!!
            auto values_cell = (values_bis.block(0,counter,number_elements,1)).col(0);
            T tmp = values_cell.dot( cb.eval_basis(pt) );
            return tmp;
        }
        */
        
        
    }
    
 // IT WORKS ONLY FOR NOT-AGGLOMERATED MESHES
    T operator()( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::face_type& fc ) const
       {
           // MATRIX NOTATION
           auto counter_face = offset(msh,fc);
           size_t counter_cell;
           // da fc devo trovare la cella in cui sono per la base
           //std::cout<<"Face number "<<counter_face<<std::endl;

           // ATTENTION, ALL THIS WORKS IN STRUCTURED QUADRANGULAR MESHES
           
           // Check if I am in the last row, upper faces (ordered differently)
           if(counter_face>=last_row_init && counter_face<=last_row_end)
           {
               counter_cell = (Ny-1)*Nx + counter_face%(last_row_init);
           }
           else
           {
            // Find in what row the face is
             auto  num_cell_row = floor(counter_face/(number_faces_one_row));
               if ( counter_face!= ( (2*Nx)*(num_cell_row+1)+num_cell_row ) )
               {
                // Faces not on the right boudary, to know in what cell are, it is sufficient to analyse the low and left face of each quadrangular cell.
                   counter_cell = floor( (counter_face-num_cell_row)/2.0 );
               }
               else
               {
                // Face on the right boudary,
                   counter_cell = ( num_cell_row+1 )*Nx -1;
               }
           
           }
           //std::cout<<"Face->Cell number "<<counter_cell<<std::endl;
           auto cl = msh.cells.at(counter_cell);
           //cell_basis_Qk<Mesh,T> cb(msh, cl, degree_FEM);
           cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
           auto values_cell = (values_bis.block(0,counter_cell,number_elements,1)).col(0);
           T tmp = values_cell.dot( cb.eval_basis(pt) );
           return tmp;

           
       }
    
    

    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;

        ret = gradient(pt);
        return ret/ret.norm();
    }
    
    
    Eigen::Matrix<T,2,1> gradient(const point<T,2>& pt) const
    {
        size_t counter=0;
        //std::cout<<"I AM IN GRADIENT SLOW !!!!"<<std::endl;
        Eigen::Matrix<T,2,1> ret = Matrix<T, 2, 1>::Zero(2, 1);
        for( const auto& cl:msh.cells)
        {
            // std::cout<<"pt_in_cell gradient in slow levelset"<<std::endl;
            if(pt_in_cell<T,Mesh>(msh,pt,cl))
            {
                //cell_basis_Qk<Mesh,T> cb(msh, cl, degree_FEM);
                cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
                
                auto values_cell = values_bis.col(counter);
                auto grad_eval =  cb.eval_gradients(pt);
                ret(0) = values_cell.dot( grad_eval.col(0) );
                // std::cout<<"Value of derivative new along x"<<ret(0)<<std::endl;
                ret(1) = values_cell.dot( grad_eval.col(1) );

                return ret;
            }
            //counter+=number_elements; // OLD VERSION OF THE CODE
            counter+=1;

        }
        std::cout<<"Se compare questo problema in gradient()"<<std::endl;
        ret(0)+=1e10;
        ret(1)+=1e10;
        return ret; //to check if doesn't enter in the loop

    }

    
    Eigen::Matrix<T,2,1> gradient( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl ) const
    {
    
        // MATRIX NOTATION
        size_t counter = offset(msh,cl);
        //std::cout<<"the cell in NEW is the number "<<counter<<std::endl;
        //std::cout<<"Value of offset "<<counter<<std::endl;
        //cell_basis_Qk<Mesh,T> cb(msh, cl, degree_FEM);
        Eigen::Matrix<T,2,1> ret = Matrix<T, 2, 1>::Zero(2, 1);
        cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
        auto values_cell = values_bis.col(counter);
        auto grad_eval =  cb.eval_gradients(pt);
        ret(0) = values_cell.dot( grad_eval.col(0) );
       // std::cout<<"Value of derivative new along x"<<ret(0)<<std::endl;
        ret(1) = values_cell.dot( grad_eval.col(1) );
        //values_cell.dot( grad_eval.col(1) );
       // std::cout<<"Value of derivative new along y"<<ret(1)<<std::endl;
        return ret;
        
        // VECTOR NOTATION
        
        //size_t counter = offset(msh,cl)*number_elements;
        //cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
        //auto values_cell = (values_bis1.segment(counter,number_elements));
        //T tmp = values_cell.dot( cb.eval_basis(pt) );
        //return tmp;
        
        
    }


    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt, const Mesh& msh, const typename Mesh::cell_type& cl) const
    {
        Eigen::Matrix<T,2,1> ret;
        ret = gradient(pt,msh,cl);
        return ret/ret.norm();
        
    }
    
    
    
    
    
    // template<typename T , typename MATRIX , typename VECTOR >
    void smooth_cut_off( T C , T x_centre , T y_centre , T radius )
    {
        T hx = params.hx();
        T hy = params.hy();
        cut_level = C;
        // T r0 = radius + 2*sqrt(pow(hx,2) + pow(hy,2)); // check value
        //T r0 = radius + 2.5*std::max(hx , hy ); // check value
        //T r0 = radius + 2.5*0.0625; // check value
        T pos_r0 = 0.45; //std::min(x_centre , 1 - x_centre );
        // IF I wanna define like that I need to pass the value all'inizio, sennò cambia per valure le exact solution (cambiando x_centre).
        
        //T max_r0 = pos_r0-(pos_r0 - 2*0.0625); //2.5*std::max(hx,hy); // if mesh too much coarse
        //r0 = std::max(r0 , max_r0);
        //r0 = 0.25*max_r0 + 0.75*r0 ;
        //if(r0<radius)
          //  r0 += 0.0625 ;
        T dist = pos_r0 - radius + 2*0.0625;
        T r0 = radius + dist/2;
        T delta = r0/8;
        //T delta = 0.0625 ; // r0/20.0; // or better std::max(hx,hy) ??
        //T delta = sqrt(pow(hx,2) + pow(hy,2));
        //if(std::abs(r0-delta-radius)<2*std::max(hx,hy) )
         //   delta = delta/2;
       
        std::cout<<"radius = "<<radius<<" , r0 = "<<r0<<" , delta = "<<delta<<" , hx = hy = "<<hx<<std::endl;
        std::cout<<"value in alfa in r_int = "<<(radius-r0)/delta<<std::endl;
        std::cout<<"value in alfa in R = "<<(pos_r0-r0)/delta<<std::endl;
        auto alfa = [x_centre , y_centre , delta , r0](const typename Mesh::point_type& pt) { // sol
            return (1 - tanh( (sqrt(pow((pt.x()-x_centre),2) + pow((pt.y()-y_centre),2)) - r0 ) / delta ))/2;};
 
        size_t counter = 0;
       // Eigen::Matrix<T, Dynamic, 1> tmp_vert = Matrix<T, Dynamic, 1>::Zero(vertices.rows(), 1); ;
        for(auto& pt:msh.points )
        {
           // std::cout<<"the point is "<< pt<<" alfa is "<<alfa(pt)<<std::endl;
            vertices(counter) = ( 1 - alfa(pt)  )*C + alfa(pt)*vertices(counter);
            counter++;
        }
       
        converting_into_HHO_formulation(vertices);
        phi_max = vertices.maxCoeff();
        phi_min = vertices.minCoeff(); // DEVO CAMBIARLO O VA BENE COSI?
        
        postprocess_output<double> postoutput00;
        typedef typename Mesh::point_type       point_type;
        point<double,2> node;
        auto alfa_values = std::make_shared< gnuplot_output_object<double> >("alfa.dat");
        //auto interface_pos = std::make_shared< gnuplot_output_object<double> >("interface_alfa.dat");
        for(auto& pt:msh.points )
               {
                  alfa_values->add_data(pt,alfa(pt));
               }
        postoutput00.add_object(alfa_values);
        postoutput00.write();
    }
    
    
    
    
   // template<typename T , typename MATRIX , typename VECTOR >
    void cut_off( T d )
    {
        //auto vertices_abs = vertices.cwiseAbs();
        //auto phi_max_abs = d*vertices_abs.maxCoeff();
        //std::cout<<"MAX VALUE OF PHI ABS"<<phi_max_abs<<std::endl;
        
        cut_level = d ;
        //std::cout<<"CUTTING AT d = "<<d<<std::endl;
        T level_set_max = vertices.maxCoeff();
        //std::cout<<"MAX VALUE OF PHI BEFORE CUTTING = "<<level_set_max<<std::endl;
        //assert(degree_FEM == 1)
        
        Eigen::Matrix<T, Dynamic, Dynamic> One_mat = Eigen::Matrix<T, Dynamic, Dynamic>::Ones(values_bis.rows(), values_bis.cols() ); // MATRIX NOTATION
           //Eigen::Matrix<T, Dynamic, 1> values_bis1; // VECTOR NOTATION
        Eigen::Matrix<T, Dynamic, 1> One_vec = Eigen::Matrix<T, Dynamic, 1>::Ones( vertices.rows(), 1 ); // saving level_set on vertices mesh
        
        
       auto cut_off_level_vec = d*One_vec ;
       auto cut_off_level_mat = d*One_mat ;
        
        vertices/=(level_set_max);
        values_bis/=(level_set_max);
       // std::cout<<"MAX VALUE OF NORMALISED PHI = "<<vertices.maxCoeff()<<std::endl;
        auto vertices_prova = vertices;
        auto values_bis_prova = values_bis ;
        
        vertices = vertices.cwiseMin(cut_off_level_vec);
        values_bis = values_bis.cwiseMin(cut_off_level_mat);
        // Cut off also in inner domain
        vertices = vertices.cwiseMax(-cut_off_level_vec);
        values_bis = values_bis.cwiseMax(-cut_off_level_mat);
        
        // NOT NORMALISED CUT
        //vertices = vertices.cwiseMin(cut_off_level_vec*level_set_max);
        //values_bis = values_bis.cwiseMin(cut_off_level_mat*level_set_max);
        
        
        phi_max = vertices.maxCoeff();
        phi_min = vertices.minCoeff(); // DEVO CAMBIARLO O VA BENE COSI?
        //std::cout<<"MAX VALUE OF PHI_CUT AND NORMALISED = "<<phi_max<<std::endl;
        std::cout<<"Cut at "<<phi_max<<" , MIN VALUE OF PHI_CUT AND NORMALISED = "<<phi_min<<" (NEGATIVE cut off ACTIVED)."<<std::endl;
        
        /*
        postprocess_output<double> postoutput2;
        typedef typename Mesh::point_type       point_type;
        point<double,2> node;
        auto test_activated_nodes = std::make_shared< gnuplot_output_object<double> >("test_activated_nodes.dat");
        auto test_before_cut_off = std::make_shared< gnuplot_output_object<double> >("test_before_cut_off.dat");
        
        for (size_t i = 0; i<vertices_prova.rows(); i++) {
            node = msh.points.at(i);
            test_before_cut_off->add_data(node,vertices_prova(i));
    
            if( std::abs(vertices_prova(i) - vertices(i))>1e-10 )
                vertices_prova(i) = 1;
            else
                vertices_prova(i) = 0;
            test_activated_nodes->add_data(node,vertices_prova(i));
            
            
        }
        
      
        for (size_t j = 0; j<values_bis.cols(); j++)
        {
            bool active_cell= FALSE;
            size_t i = 0;
            while (i<values_bis.rows() && !active_cell ) {
                if( values_bis_prova(i,j)!= values_bis(i,j) ){
                    active_cell = TRUE;
                //std::cout<<"values_bis_prova(i,j)= "<<values_bis_prova(i,j)<<" , values_bis(i,j)= "<<values_bis(i,j)<<std::endl;
              //  std::cout<<" the cell num "<<j<<" is active."<<std::endl;
                }
                i++;
            }
        }
        
        
        
        postoutput2.add_object(test_before_cut_off);
        postoutput2.add_object(test_activated_nodes);
        postoutput2.write();
        */
        /*
        for (size_t i = 0 ; i<values_bis.rows() ; i++ ) {
            for (size_t j = 0 ; j<values_bis.cols() ; j++ ) {
                
            }
            
        }
            */
        
    }

    void boundary_con( T d )
    {
        T level_set_max = vertices.maxCoeff();
        std::cout<<"MAX VALUE OF PHI"<<level_set_max<<std::endl;
        //assert(degree_FEM == 1)
        
        Eigen::Matrix<T, Dynamic, Dynamic> One_mat = Eigen::Matrix<T, Dynamic, Dynamic>::Ones(values_bis.rows(), values_bis.cols() ); // MATRIX NOTATION
           //Eigen::Matrix<T, Dynamic, 1> values_bis1; // VECTOR NOTATION
        Eigen::Matrix<T, Dynamic, 1> One_vec = Eigen::Matrix<T, Dynamic, 1>::Ones( vertices.rows(), 1 ); // saving level_set on vertices mesh
        
        
       auto cut_off_level_vec = d*One_vec ;
       auto cut_off_level_mat = d*One_mat ;
        
        vertices/=(level_set_max);
        values_bis/=(level_set_max);
        vertices = vertices.cwiseMin(cut_off_level_vec);
        values_bis = values_bis.cwiseMin(cut_off_level_mat);
        /*
        for (size_t i = 0 ; i<values_bis.rows() ; i++ ) {
            for (size_t j = 0 ; j<values_bis.cols() ; j++ ) {
                
            }
            
        }
            */
        
    }
    
    
    
    
};


template<typename T>
class velocity_field
{
   public:
    virtual Eigen::Matrix<T,2,1> operator()(const point<T,2>& pt) const
    {
    }
    
    /*
     virtual Eigen::Matrix<T,2,1> flux(const point<T,2>& pt) const
    {
    }
    
    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;

        ret = gradient(pt);
        return ret/ret.norm();
    }
    */
   
};



template<typename T>
struct rotational_velocity_field: public velocity_field<T>
{
    T u1, xc , yc;

    rotational_velocity_field(T u1 )
    : u1(u1){}

    rotational_velocity_field(T xc , T yc , T u1 )
    : u1(u1) , xc(xc) , yc(yc)
    {}

    Eigen::Matrix<T,2,1> operator()(const point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;
        ret(0) = (  u1*pt.y() -xc )/xc ;
        ret(1) = (- u1*pt.x() + yc ) / yc ;
       
        return ret;
    }
    
    
};

template<typename T>
struct linear_velocity_field: public velocity_field<T>
{
    T u1, u2 , u3 , u4;

    linear_velocity_field(T u1 , T u2 , T u3 , T u4 )
        : u1(u1), u2(u2) , u3(u3) , u4(u4)
    {}

    Eigen::Matrix<T,2,1> operator()(const point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;
        ret(0) = u1*pt.x() + u2;
        ret(1) = u3*pt.y() + u4;
       
        return ret;
    }
/*
    Eigen::Matrix<T,2,1> gradient(const point<T,2>& pt) const
    {
        Eigen::Matrix<T,2,1> ret;
        ret(0) = 2*pt.x() - 2*alpha;
        ret(1) = 2*pt.y() - 2*beta;
        return ret;
    }
    */
};




template< typename Mesh , typename Fonction , typename FiniteSpace , typename T = typename Mesh::coordinate_type  >
struct projection_velocity_high_order: public velocity_field<T>
{
  
     size_t  last_row_init, last_row_end, number_faces_one_row;
    
    
     // NEW IMPLEMENTATION STUFF
     Mesh msh; // Original mesh, NOT agglomerated.
     mesh_init_params<T> params; // mesh parameter
     
     size_t degree_FEM ; // FEM degree
     size_t n_cls ; // #cells
     size_t local_dim; // Local Dimension (degree_FEM+1)*(degree_FEM+1)
     size_t n_vertices ; // #vertices
      size_t      Nx, Ny ; // Number of cells in x and y direciton
     // connectivity matrix : for each cell, it stores the global numbering
     std::vector< std::vector<std::pair<size_t,bool>>> connectivity_matrix ;
     //std::vector<std::vector<size_t>> connectivity_matrix ;
     
     size_t dim_HHO; // Global dimension Discontinuous framework = Local dimension * #cells
     size_t ndof_FE; // Global dimension FE continuous = #nodes
    
     T u_max0 , u_max1 ;
     
     
     //SparseMatrix<T>                 Global_Mass; // Global mass, saved for FEM problem
     //Matrix<T, Dynamic, 1>           RHS0;    // Known term 1
     //Matrix<T, Dynamic, 1>           RHS1;    // Known term 1
     //std::vector< Triplet<T> >       triplets; // Position elements: Sparse Matrix Notation
     
     std::pair<Eigen::Matrix<T, Dynamic, Dynamic> , Eigen::Matrix<T, Dynamic, Dynamic>> sol_HHO ; // projection saved in HHO format: cell by cell
     std::pair<Matrix<T, Dynamic, 1>,Matrix<T, Dynamic, 1>> sol_FEM ; // projection saved in Continuos FE format: global nodes
     //std::pair<Eigen::Matrix<T, Dynamic, 1>,Matrix<T, Dynamic, 1>> vertices; // saving level_set on vertices mesh

     
    
     projection_velocity_high_order(const FiniteSpace& fe_data , const Mesh & msh )
           : degree_FEM(fe_data.order) , local_dim(fe_data.local_ndof), msh(msh), Nx(fe_data.Nx),Ny(fe_data.Ny), params(fe_data.params) , dim_HHO(fe_data.ndof_disc) , n_cls(fe_data.n_cls) ,n_vertices(fe_data.n_vertices) , connectivity_matrix(fe_data.connectivity_matrix) , ndof_FE(fe_data.ndof_FE)
       {
           std::cout<<"STO USANDO: equidistriduted_nodes_ordered_bis<T,Mesh>(msh, cl, degree)"<<std::endl;
           //last_row_init = Ny*(2*Nx+1); // There are 2 faces for each row of cells + Ny
           //last_row_end = last_row_init + Nx-1;
           //number_faces_one_row = 2*Nx+1; // for each cell I count only the low and sx faces, respectevely 0-1 2-3 4-5 6-7 8 + the last on the right boundary
           
           
           //vertices = std::make_pair( Eigen::Matrix<T, Dynamic, 1>::Zero( ( (Nx+1)*(Ny+1) ), 1 ) , Eigen::Matrix<T, Dynamic, 1>::Zero( ((Nx+1)*(Ny+1)), 1 ) );

           
           //vertices = std::make_pair( Eigen::Matrix<T, Dynamic, 1>::Zero( n_vertices , 1 ),Eigen::Matrix<T, Dynamic, 1>::Zero( n_vertices , 1 ) );
           sol_FEM = std::make_pair( Matrix<T, Dynamic, 1>::Zero(ndof_FE) , Matrix<T, Dynamic, 1>::Zero(ndof_FE) );
          
           sol_HHO = std::make_pair(Eigen::Matrix<T, Dynamic, Dynamic>::Zero( local_dim, n_cls ) , Eigen::Matrix<T, Dynamic, Dynamic>::Zero( local_dim, n_cls ) ) ;

           
       }
    
    projection_velocity_high_order(const FiniteSpace& fe_data , const Fonction & u, const Mesh & msh)
        : degree_FEM(fe_data.order) , local_dim(fe_data.local_ndof), msh(msh), Nx(fe_data.Nx),Ny(fe_data.Ny), params(fe_data.params) , dim_HHO(fe_data.ndof_disc) , n_cls(fe_data.n_cls) ,n_vertices(fe_data.n_vertices) , connectivity_matrix(fe_data.connectivity_matrix) , ndof_FE(fe_data.ndof_FE)
    {
        //last_row_init = Ny*(2*Nx+1); // There are 2 faces for each row of cells + Ny
        //last_row_end = last_row_init + Nx-1;
        //number_faces_one_row = 2*Nx+1; // for each cell I count only the low and sx faces, respectevely 0-1 2-3 4-5 6-7 8 + the last on the right boundary
     
        // Saving the projection in HHO discontinuous format (MATRIX NOTATION)
         //vertices = std::make_pair( Eigen::Matrix<T, Dynamic, 1>::Zero( n_vertices , 1 ),Eigen::Matrix<T, Dynamic, 1>::Zero( n_vertices , 1 ) );
         sol_FEM = std::make_pair( Matrix<T, Dynamic, 1>::Zero(ndof_FE) , Matrix<T, Dynamic, 1>::Zero(ndof_FE) );
         
         sol_HHO = std::make_pair(Eigen::Matrix<T, Dynamic, Dynamic>::Zero( local_dim, n_cls ) , Eigen::Matrix<T, Dynamic, Dynamic>::Zero( local_dim, n_cls ) ) ;

        
        size_t i_global = 0 , i_local=0 , i_vertex=0;
        for(auto& cl:msh.cells)
        {
            auto qps = equidistriduted_nodes_ordered_bis<T,Mesh>(msh, cl, degree_FEM);
            i_local = 0;
            for ( const auto & qp : qps)
            {
                auto asm_map = connectivity_matrix[i_global][i_local].first ;
                auto u0 = (u(qp))(0) ;
                auto u1 = (u(qp))(1) ;
                
                sol_HHO.first(i_local,i_global) = u0 ;
                sol_HHO.second(i_local,i_global) = u1 ;
                sol_FEM.first(asm_map) = u0 ;
                sol_FEM.second(asm_map) = u1 ;
           /*
                i_vertex = i_global+floor(i_global/Nx);
                if( i_local==0 ){
                    vertices.first(i_vertex) = u0 ;
                    vertices.second(i_vertex) = u1 ;
                }
            
                if( i_local==1 ){
                    vertices.first(i_vertex+1) = u0 ;
                    vertices.second(i_vertex+1) = u1 ;
                }
            
                if( i_local==(degree_FEM+2) ){
                    vertices.first(i_vertex+Nx+2) = u0 ;
                    vertices.second(i_vertex+Nx+2) = u1 ;
                }
                
                if( i_local==(degree_FEM+1) ){
                    vertices.first(i_vertex+Nx+1) = u0 ;
                    vertices.second(i_vertex+Nx+1) = u1 ;
                }
            */
                
                i_local++;
                
            
            }
            i_global++;
            
            
        }
        
    }
    

    template< typename MATRIX >
    void  set_discrete_points( MATRIX& values_new)
    {
        sol_HHO = values_new ;
        std::cout<<" Using set_discrete_points check that sol_FEM already uploaded!"<<std::endl;
       
    }
    
    template< typename VECTOR >
    void converting_into_HHO_formulation( const VECTOR& values_new )
    {
        
        for(size_t counter_bis = 0 ; counter_bis < n_cls ;counter_bis++)
        {
            for (size_t i = 0; i < local_dim; i++){
                size_t asm_map =  connectivity_matrix[counter_bis][i].first ;
                sol_HHO.first(i,counter_bis) = values_new.first( asm_map );
                sol_HHO.second(i,counter_bis) = values_new.second( asm_map );
            }
            //  size_t i_vertex = counter_bis+floor(counter_bis/Nx);
            //  vertices(i_vertex) = sol_HHO(0,counter_bis) ;
            //  vertices(i_vertex+1) = sol_HHO(1,counter_bis) ;
            // vertices(i_vertex+Nx+2) = sol_HHO(2,counter_bis) ;
            //  vertices(i_vertex+Nx+1) = sol_HHO(3,counter_bis) ;
        }
        std::cout<<" Using converting_into_HHO_formulation check that sol_FEM already uploaded!"<<std::endl;
        //set_max_min();
        //phi_min = sol_HHO.minCoeff() ;
        //phi_max = sol_HHO.maxCoeff() ;
                 
    }
    
    template< typename MATRIX >
    void converting_into_FE_formulation( const MATRIX& values_new )
    {
      
        for(size_t counter_bis = 0 ; counter_bis < n_cls ;counter_bis++)
        {
            for (size_t i = 0; i < local_dim; i++){
                size_t asm_map =  connectivity_matrix[counter_bis][i].first ;
                sol_FEM.first( asm_map ) = values_new.first(i,counter_bis) ;
                sol_FEM.second( asm_map ) = values_new.second(i,counter_bis) ;
            }
            //  size_t i_vertex = counter_bis+floor(counter_bis/Nx);
            //  vertices(i_vertex) = sol_HHO(0,counter_bis) ;
            //  vertices(i_vertex+1) = sol_HHO(1,counter_bis) ;
            // vertices(i_vertex+Nx+2) = sol_HHO(2,counter_bis) ;
            //  vertices(i_vertex+Nx+1) = sol_HHO(3,counter_bis) ;
        }
        std::cout<<" Using converting_into_FE_formulation check that sol_HHO already uploaded!"<<std::endl;
        //set_max_min();
        //phi_min = sol_FEM.minCoeff() ;
        //phi_max = sol_FEM.maxCoeff() ;
        
    }
       

    
    std::pair<T,T> operator()( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl) const
    {
        size_t counter = offset(msh,cl) ;
        cell_basis_Lagrangian_ordered<Mesh,T> cb(msh, cl, degree_FEM);
        //std::cout<<"cb u proj high order"<<'\n'<<cb.eval_basis(pt)<<'\n'<<std::endl;
        auto values_cell_first = (sol_HHO.first.block(0,counter,local_dim,1)).col(0);
        auto values_cell_second = (sol_HHO.second.block(0,counter,local_dim,1)).col(0);
        //std::cout<<"values_cell_first"<<'\n'<<values_cell_first<<'\n'<<std::endl;
        //std::cout<<"values_cell_second"<<'\n'<<values_cell_second<<'\n'<<std::endl;
        T tmp1 = values_cell_first.dot( cb.eval_basis(pt) );
        T tmp2 = values_cell_second.dot( cb.eval_basis(pt) );
        //std::cout<<"tmp1"<<" "<<tmp1<<'\n'<<"tmp2"<<" "<<tmp2<<std::endl;
        
        return std::make_pair(tmp1,tmp2);
    }
    
    
       
};


template< typename Mesh , typename FiniteSpace , typename T = typename Mesh::coordinate_type  >
struct velocity_high_order
{
  
     size_t  last_row_init, last_row_end, number_faces_one_row;
    
    
     // NEW IMPLEMENTATION STUFF
     Mesh msh; // Original mesh, NOT agglomerated.
     mesh_init_params<T> params; // mesh parameter
     
     size_t degree_FEM ; // FEM degree
     size_t n_cls ; // #cells
     size_t local_dim; // Local Dimension (degree_FEM+1)*(degree_FEM+1)
     size_t n_vertices ; // #vertices
      size_t      Nx, Ny ; // Number of cells in x and y direciton
     // connectivity matrix : for each cell, it stores the global numbering
     std::vector< std::vector<std::pair<size_t,bool>>> connectivity_matrix ;
     //std::vector<std::vector<size_t>> connectivity_matrix ;
     
     size_t dim_HHO; // Global dimension Discontinuous framework = Local dimension * #cells
     size_t ndof_FE; // Global dimension FE continuous = #nodes
    
     T u_max0 , u_max1 ;
     
     
     //SparseMatrix<T>                 Global_Mass; // Global mass, saved for FEM problem
     //Matrix<T, Dynamic, 1>           RHS0;    // Known term 1
     //Matrix<T, Dynamic, 1>           RHS1;    // Known term 1
     //std::vector< Triplet<T> >       triplets; // Position elements: Sparse Matrix Notation
     
     std::pair<Eigen::Matrix<T, Dynamic, Dynamic> , Eigen::Matrix<T, Dynamic, Dynamic>> sol_HHO ; // projection saved in HHO format: cell by cell
     std::pair<Matrix<T, Dynamic, 1>,Matrix<T, Dynamic, 1>> sol_FEM ; // projection saved in Continuos FE format: global nodes
     //std::pair<Eigen::Matrix<T, Dynamic, 1>,Matrix<T, Dynamic, 1>> vertices; // saving level_set on vertices mesh

     
    
     velocity_high_order(const FiniteSpace& fe_data , const Mesh & msh )
           : degree_FEM(fe_data.order) , local_dim(fe_data.local_ndof), msh(msh), Nx(fe_data.Nx),Ny(fe_data.Ny), params(fe_data.params) , dim_HHO(fe_data.ndof_disc) , n_cls(fe_data.n_cls) ,n_vertices(fe_data.n_vertices) , connectivity_matrix(fe_data.connectivity_matrix) , ndof_FE(fe_data.ndof_FE)
       {
           std::cout<<"STO USANDO: equidistriduted_nodes_ordered_bis<T,Mesh>(msh, cl, degree)"<<std::endl;
           //last_row_init = Ny*(2*Nx+1); // There are 2 faces for each row of cells + Ny
           //last_row_end = last_row_init + Nx-1;
           //number_faces_one_row = 2*Nx+1; // for each cell I count only the low and sx faces, respectevely 0-1 2-3 4-5 6-7 8 + the last on the right boundary
           
           
           //vertices = std::make_pair( Eigen::Matrix<T, Dynamic, 1>::Zero( ( (Nx+1)*(Ny+1) ), 1 ) , Eigen::Matrix<T, Dynamic, 1>::Zero( ((Nx+1)*(Ny+1)), 1 ) );

           
           //vertices = std::make_pair( Eigen::Matrix<T, Dynamic, 1>::Zero( n_vertices , 1 ),Eigen::Matrix<T, Dynamic, 1>::Zero( n_vertices , 1 ) );
           sol_FEM = std::make_pair( Matrix<T, Dynamic, 1>::Zero(ndof_FE) , Matrix<T, Dynamic, 1>::Zero(ndof_FE) );
          
           sol_HHO = std::make_pair(Eigen::Matrix<T, Dynamic, Dynamic>::Zero( local_dim, n_cls ) , Eigen::Matrix<T, Dynamic, Dynamic>::Zero( local_dim, n_cls ) ) ;

           
       }
    
    
    

    template< typename MATRIX >
    void  set_discrete_points( MATRIX& values_new)
    {
        sol_HHO = values_new ;
        std::cout<<" Using set_discrete_points check that sol_FEM already uploaded!"<<std::endl;
       
    }
    
    template< typename VECTOR >
    void converting_into_HHO_formulation( const VECTOR& values_new )
    {
        
        for(size_t counter_bis = 0 ; counter_bis < n_cls ;counter_bis++)
        {
            for (size_t i = 0; i < local_dim; i++){
                size_t asm_map =  connectivity_matrix[counter_bis][i].first ;
                sol_HHO.first(i,counter_bis) = values_new.first( asm_map );
                sol_HHO.second(i,counter_bis) = values_new.second( asm_map );
            }
            //size_t i_vertex = counter_bis+floor(counter_bis/Nx);
            //vertices(i_vertex) = sol_HHO(0,counter_bis) ;
            //vertices(i_vertex+1) = sol_HHO(1,counter_bis) ;
            //vertices(i_vertex+Nx+2) = sol_HHO(2,counter_bis) ;
            //vertices(i_vertex+Nx+1) = sol_HHO(3,counter_bis) ;
        }
        std::cout<<" Using converting_into_HHO_formulation check that sol_FEM already uploaded!"<<std::endl;
        //set_max_min();
        //phi_min = sol_HHO.minCoeff() ;
        //phi_max = sol_HHO.maxCoeff() ;
                 
    }
    
    template< typename MATRIX >
    void converting_into_FE_formulation( const MATRIX& values_new )
    {
      
        for(size_t counter_bis = 0 ; counter_bis < n_cls ;counter_bis++)
        {
            for (size_t i = 0; i < local_dim; i++){
                size_t asm_map =  connectivity_matrix[counter_bis][i].first ;
                sol_FEM.first( asm_map ) = values_new.first(i,counter_bis) ;
                sol_FEM.second( asm_map ) = values_new.second(i,counter_bis) ;
            }
           //size_t i_vertex = counter_bis+floor(counter_bis/Nx);
            //vertices(i_vertex) = sol_HHO(0,counter_bis) ;
            //vertices(i_vertex+1) = sol_HHO(1,counter_bis) ;
            //vertices(i_vertex+Nx+2) = sol_HHO(2,counter_bis) ;
            //vertices(i_vertex+Nx+1) = sol_HHO(3,counter_bis) ;
        }
        std::cout<<" Using converting_into_FE_formulation check that sol_HHO already uploaded!"<<std::endl;
        //set_max_min();
        //phi_min = sol_FEM.minCoeff() ;
        //phi_max = sol_FEM.maxCoeff() ;
        
    }
       

    
    std::pair<T,T> operator()( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl) const
    {
        size_t counter = offset(msh,cl) ;
        cell_basis_Lagrangian_ordered<Mesh,T> cb(msh, cl, degree_FEM);
        //std::cout<<"cb u proj high order"<<'\n'<<cb.eval_basis(pt)<<'\n'<<std::endl;
        auto values_cell_first = (sol_HHO.first.block(0,counter,local_dim,1)).col(0);
        auto values_cell_second = (sol_HHO.second.block(0,counter,local_dim,1)).col(0);
        //std::cout<<"values_cell_first"<<'\n'<<values_cell_first<<'\n'<<std::endl;
        //std::cout<<"values_cell_second"<<'\n'<<values_cell_second<<'\n'<<std::endl;
        T tmp1 = values_cell_first.dot( cb.eval_basis(pt) );
        T tmp2 = values_cell_second.dot( cb.eval_basis(pt) );
        //std::cout<<"tmp1"<<" "<<tmp1<<'\n'<<"tmp2"<<" "<<tmp2<<std::endl;
        
        return std::make_pair(tmp1,tmp2);
    }
    
    
       
};




template< typename Mesh>
struct Current_Mesh
{
    Mesh current_mesh ;
    Current_Mesh(const Mesh& msh):current_mesh(msh){}
    
    void set_current_mesh(const Mesh& msh)
    {
        current_mesh = msh;
    }
    
    
};



template< typename T , typename Mesh ,typename Level_Set,typename Fonction >
struct LS_cell: public projected_level_set< T,  Mesh , Fonction>
{
    typedef typename Mesh::cell_type       cell_type;
    cell_type agglo_LS_cl;
    std::vector<cell_type> subcells;
    Mesh agglo_msh;
    Level_Set level_set;
    //LS_cell(const Level_Set & level_set, const Mesh & msh, const typename Mesh::cell_type& cl)
   // : agglo_cl(cl), agglo_msh(msh), level_set(level_set){}
    // I don't know if I have to define a copyconstructor for level_set.. TO BE CHECKED!
    LS_cell(const Level_Set & level_set_, const Mesh & msh)
    : agglo_msh(msh), level_set(level_set_){}
    //LS_cell(const Level_Set & level_set )
    //: level_set(level_set){}

    //LS_cell()=default;
    
    T operator()(const point<T,2>& pt) const
       {
           if (subcells.size()<1)
           {
               assert(agglo_LS_cl.user_data.offset_subcells.size()==1);
               auto offset_old = agglo_LS_cl.user_data.offset_subcells[0];
               auto cl_old = level_set.msh.cells[offset_old];
               return level_set( pt , level_set.msh , cl_old );
           }
       
           else
           {
               //std::cout<<"Ls operator"<<std::endl;
               auto offset = pt_in_subcell(level_set.msh,pt,agglo_LS_cl);
               auto subcl = level_set.msh.cells[offset];
               return level_set( pt , level_set.msh , subcl );
           }
               
       }
    
    
    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt) const
    {
        if (subcells.size()<1)
            {
                assert(agglo_LS_cl.user_data.offset_subcells.size()==1);
                auto offset_old = agglo_LS_cl.user_data.offset_subcells[0];
                auto cl_old = level_set.msh.cells[offset_old];
                return level_set.normal( pt , level_set.msh , cl_old );
            }
        
            else
            {
               
                auto offset = pt_in_subcell(level_set.msh,pt,agglo_LS_cl);
                auto subcl = level_set.msh.cells[offset];
                return level_set.normal( pt , level_set.msh , subcl );
            }
    }
    
    Eigen::Matrix<T,2,1> gradient(const point<T,2>& pt) const
    {
        if (subcells.size()<1)
            {
                std::cout<<"I m here 2.5 and size is "<<agglo_LS_cl.user_data.offset_subcells.size()<<std::endl;
                assert(agglo_LS_cl.user_data.offset_subcells.size()==1);
                auto offset_old = agglo_LS_cl.user_data.offset_subcells[0];
                auto cl_old = level_set.msh.cells[offset_old];
                return level_set.gradient( pt , level_set.msh , cl_old );
            }
        
            else
            {
                std::cout<<"Ls gradient"<<std::endl;
                auto offset = pt_in_subcell(level_set.msh,pt,agglo_LS_cl);
                auto subcl = level_set.msh.cells[offset];
                return level_set.gradient( pt , level_set.msh , subcl );
            }
        
    }
  
    void cell_assignment(const cell_type& cl)
    {
        subcells.clear();
        agglo_LS_cl = cl;
        if (agglo_LS_cl.user_data.offset_subcells.size()>1)
        {

            for (auto& offset_subcells:agglo_LS_cl.user_data.offset_subcells)
            {
                subcells.push_back( level_set.msh.cells[offset_subcells] );
            
            }
        }
            
    }
    
    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt , const cell_type& cl)
    {
        agglo_LS_cl = cl;
        return normal( pt );
    }
    
    T operator()(const point<T,2>& pt,const cell_type& cl )
    {
        agglo_LS_cl = cl;
        return operator()( pt );
    }
    
    
  //  void get_cell()
  //  {
  //      auto crr_cl = download_current_cell();
  //      std::cout<<"IN get_cell.. size cell "<<crr_cl.user_data.offset_subcells.size()<<std::endl;
  //      this->cell_assignment(crr_cl);
 //   }
    
};





/**************MOVING INTERFACE: LEVEL SET METHOD  **************/

template<typename T, typename Mesh >
struct projection
{
    Eigen::Matrix<T, Dynamic, Dynamic> values_bis;
    Eigen::Matrix<T, Dynamic, 1> vertices; // saving level_set on vertices mesh
    
    size_t degree_FEM;
    size_t number_elements;
    Mesh msh;
    size_t      Nx, Ny ;
    mesh_init_params<T> params;
    //size_t  last_row_init, last_row_end, number_faces_one_row;

    projection( const Mesh & msh, size_t degree_k , const mesh_init_params<T>& params)
        : number_elements((degree_k+1)*(degree_k+1)), msh(msh),degree_FEM(degree_k),Nx(params.Nx),Ny(params.Ny), params(params)
    {
        vertices = Eigen::Matrix<T, Dynamic, 1>::Zero( ((Nx+1)*(Ny+1)), 1 );
        values_bis= Eigen::Matrix<T, Dynamic, Dynamic>::Zero(number_elements, msh.cells.size());
    }
    
    projection()=default;
    
    
    void  set_discrete_points( Eigen::Matrix<T, Dynamic, Dynamic>& values_new)
    {
        values_bis=values_new;
    }
    
    
    void converting_into_HHO_formulation( const Eigen::Matrix<T, Dynamic, 1>& vertices )
    {
        assert(degree_FEM == 1);
        for (size_t j = 0; j < values_bis.cols() ; j++)
        {

            for (size_t i = 0; i < values_bis.rows() ; i++)
            {
                auto i_vertex = j + floor(j/Nx);
                if ( i == 0 )
                    values_bis(i,j) = vertices(i_vertex);
                
                if( i ==1 )
                    values_bis(i,j) = vertices(i_vertex+1) ;
                
                if( i ==(degree_FEM+2) )
                    values_bis(i,j) = vertices(i_vertex+Nx+2) ;
                
                if( i ==(degree_FEM+1) )
                    values_bis(i,j) = vertices(i_vertex+Nx+1) ;
            }
        }
                 
    }
    
    void converting_into_FE_formulation( const Eigen::Matrix<T, Dynamic, Dynamic>& values )
    {
        assert(degree_FEM == 1);
        for (size_t j = 0; j < values_bis.cols() ; j++)
        {

            for (size_t i = 0; i < values_bis.rows() ; i++)
            {
                auto i_vertex = j + floor(j/Nx);
                if ( i == 0 )
                    vertices(i_vertex) = values(i,j) ;
                
                if( i ==1 )
                    vertices(i_vertex+1) = values(i,j) ;
                
                if( i ==(degree_FEM+2) )
                    vertices(i_vertex+Nx+2) = values(i,j) ;
                
                if( i ==(degree_FEM+1) )
                    vertices(i_vertex+Nx+1) = values(i,j) ;
            }
        }
                 
    }
       
    
    
   
    T operator()( const typename Mesh::node_type& node ) const
    {
        return vertices(node.ptid);
    }
    
    T operator()( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl ) const
    {
        size_t counter = offset(msh,cl) ;
        cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
        auto values_cell = (values_bis.block(0,counter,number_elements,1)).col(0);
        T tmp = values_cell.dot( cb.eval_basis(pt) );
        return tmp;
    }
  
    Eigen::Matrix<T,2,1> gradient( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl ) const
    {
        size_t counter = offset(msh,cl);
        Eigen::Matrix<T,2,1> ret = Matrix<T, 2, 1>::Zero(2, 1);
        cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
        auto values_cell = values_bis.col(counter);
        auto grad_eval =  cb.eval_gradients(pt);
        ret(0) = values_cell.dot( grad_eval.col(0) );
        ret(1) = values_cell.dot( grad_eval.col(1) );
        return ret;
    }


    Eigen::Matrix<T,2,1> normal(const point<T,2>& pt, const Mesh& msh, const typename Mesh::cell_type& cl) const
    {
        Eigen::Matrix<T,2,1> ret;
        ret = gradient(pt,msh,cl);
        return ret/ret.norm();
        
    }
    
};




template <typename T>
void positive_part( Eigen::Matrix<T, Dynamic, Dynamic>& mat) {
    for (size_t i = 0; i<mat.rows();i++) {
        for (size_t j = 0; j<mat.cols();j++) {
            if( mat(i,j) < 0 )
                mat(i,j) = 0.;
        }
    }

}





/*

template <typename T>
Eigen::Matrix<T, Dynamic, Dynamic> positive_part(const Eigen::Matrix<T, Dynamic, Dynamic>& mat) {
    //Eigen::Matrix<T, Dynamic, Dynamic> ret =mat;
    Eigen::Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(mat.rows(), mat.cols());
    auto positive = mat.cwiseSign();
    std::cout<<"positive: "<<'\n'<<positive<<std::endl;
    for (size_t i = 0; i<mat.rows();i++) {
        for (size_t j = 0; j<mat.cols();j++) {
            if (positive(i,j)==-1) {
                ret(i,j) = 0;
            }
            //if( mat(i,j) < 0 )
            //    ret(i,j) = 0;
        }
    }
    return ret;
}
*/


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template<typename T, typename Mesh>
class entropy
{
public:
    virtual T operator()(const point<T,2>& pt, const typename Mesh::cell_type& cl) const
    {
    }
    
    virtual T derivative(const point<T,2>& pt, const typename Mesh::cell_type& cl) const
    {
    }
   
};

template<typename T , typename Fonction , typename Mesh>
struct non_linear_entropy: public entropy<T,Mesh>
{
    const T eps;
    Fonction phi;
    T phi_max , phi_min ;
    Mesh msh;
    
    non_linear_entropy(const T eps , const Fonction& phi ,const Mesh& msh): eps(eps),phi(phi),msh(msh)
    {
        phi_max = phi.phi_max;
        phi_min = phi.phi_min;
        
    }
    
    /// THIS WAS FOR A GENERAL PHI, RE-WRITTEN FOR PHI BETWEEN 0 AND 1
    /*
    T operator()(const Fonction& f, const point<T,2>& pt, const typename Mesh::cell_type& cl) const
    {
        return -std::log10( std::abs( (phi_max-f(pt,msh,cl))*(f(pt,msh,cl)-phi_min) ) + eps );
    }
    
    T operator()(const point<T,2>& pt , const typename Mesh::cell_type& cl) const
    {
        return -std::log10( std::abs( (phi_max-phi(pt,msh,cl))*(phi(pt,msh,cl)-phi_min) ) + eps );
    }
    
    T derivative(const point<T,2>& pt, const typename Mesh::cell_type& cl) const
    {
       // std::cout<<"il segno di "<<((phi_max-phi(pt,msh,cl))*(phi(pt,msh,cl)-phi_min ))<<" is "<< sgn( (phi_max-phi(pt,msh,cl))*(phi(pt,msh,cl)-phi_min ) )<<std::endl;
        return -1./(( std::abs( (phi_max-phi(pt,msh,cl))*(phi(pt,msh,cl)-phi_min ) ) +eps)*std::log(10)) * sgn( (phi_max-phi(pt,msh,cl))*(phi(pt,msh,cl)-phi_min ) ) * ( phi_max - 2*phi(pt,msh,cl) + phi_min ) ;
    }
    */
    
    T operator()(const Fonction& f, const point<T,2>& pt, const typename Mesh::cell_type& cl) const
    {
        return -std::log10( std::abs( f(pt,msh,cl)*(1-f(pt,msh,cl)) ) + eps );
    }
    
    T operator()(const point<T,2>& pt , const typename Mesh::cell_type& cl) const
    {
        return -std::log10( std::abs( phi(pt,msh,cl)*(1-phi(pt,msh,cl)) ) + eps );
    }
    
    T derivative(const point<T,2>& pt, const typename Mesh::cell_type& cl) const
    {
       // std::cout<<"il segno di "<<((phi_max-phi(pt,msh,cl))*(phi(pt,msh,cl)-phi_min ))<<" is "<< sgn( (phi_max-phi(pt,msh,cl))*(phi(pt,msh,cl)-phi_min ) )<<std::endl;
        return -1./( ( std::abs( (1-phi(pt,msh,cl))*phi(pt,msh,cl)  ) +eps )*std::log(10) ) * sgn( (1-phi(pt,msh,cl))*phi(pt,msh,cl)  ) * ( 1 - 2*phi(pt,msh,cl)  ) ;
    }
   
};




template<typename Mesh, typename T = typename Mesh::coordinate_type>
Matrix<T, Dynamic, Dynamic>
make_lagrange_local_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, size_t di = 0)
{
    
    //auto f_neigh = cl.user_data.f_neighbors;
    //auto d_neigh = cl.user_data.d_neighbors;
    
    cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, Dynamic> ret = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);
    //Matrix<T, Dynamic, 1> ret2 = Matrix<T, Dynamic, 1>::Zero(cbs, 1);

    auto qps = integrate(msh, cl, 2*(degree+di)); // integration of order 2k

    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        ret += qp.second * phi * phi.transpose();
        // phi * phi.transpose is (degree+1)^2 x (degree+1)^2 ; qp.second is a scalar
    }
    //ret2  =  ret.rowwise().sum(); // sum row; i.e. mi0 + mi1 + mi2 + mi3 , with i = 0 : 3
    
    return ret;
}

template<typename Mesh, typename T = typename Mesh::coordinate_type>
Matrix<T, Dynamic, Dynamic>
make_lagrange_lumped_mass_matrix(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, size_t di = 0)
{
    cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs, 1);

    auto qps = integrate(msh, cl, 2*(degree+di));

    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        ret += qp.second * phi;
    }

    return ret;
}

template<typename Mesh, typename T = typename Mesh::coordinate_type>
std::pair< Matrix<T, Dynamic, Dynamic> , Matrix<T, Dynamic, Dynamic> >
make_local_cij_matrix(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree, size_t di = 0)
{
    cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree);
    auto cbs = cb.size();

    Matrix<T, Dynamic, Dynamic> ret0 = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);
    Matrix<T, Dynamic, Dynamic> ret1 = Matrix<T, Dynamic, Dynamic>::Zero(cbs, cbs);
    
    // for Q1 , degree = 1-> integration of order 2
    auto qps = integrate(msh, cl, 2*(degree+di)); // integration of order 2k

    for (auto& qp : qps)
    {
        auto phi = cb.eval_basis(qp.first);
        
        auto phi_grad = cb.eval_gradients(qp.first);
        
        ret0 += qp.second * phi * ((phi_grad).col(0)).transpose();
    
        ret1 += qp.second * phi * ((phi_grad).col(1)).transpose();
    }
    
    return std::make_pair(ret0,ret1);
}





template < typename VECTOR, typename Mesh, typename MATRIX , typename T = typename Mesh::coordinate_type >
Matrix<T, Dynamic, 1>  make_dij_vector(const Mesh & msh, const MATRIX& dij, const VECTOR& phi_old )
{
    Matrix<T, Dynamic, 1> ret = Eigen::Matrix<T, Dynamic, 1>::Zero(dij.rows(), 1);
    
    
    for(size_t i = 0 ; i<dij.rows() ; i++ )
    {
        auto tmp = (dij.row(i)).adjoint();
        ret(i) = tmp.dot( phi_old ) - tmp(i)*( phi_old(i) );
    }
    return ret;
}

template < typename VECTOR, typename T >
Matrix<T, Dynamic, 1> solveFEM( const VECTOR& lumped_mass , const VECTOR& conv , const VECTOR& last_term, const VECTOR& phi_old , T dt , const std::vector<size_t>& bdry_nodes)
{
    /*
    Matrix<T, Dynamic, 1> sol =  Eigen::Matrix<T, Dynamic, 1>::Zero(lumped_mass.rows(), 1);
    size_t j = 0;
    for (size_t i = 0 ; i < phi_old.rows() ; i++ ) {
        if( i == bdry_nodes.at(j) ){
            sol(i) = phi_old(i);
            j++; // I can do since bdry_nodes is ordered
        }
        else
           sol(i) = phi_old(i) - dt * conv(i)/( lumped_mass(i) ) + dt * last_term(i)/( lumped_mass(i) );
    }
    return sol;
    */
    // if I dont consider boundary condition here, but just solve for all
    
    
    return phi_old - dt * conv.cwiseQuotient(lumped_mass) + dt * last_term.cwiseQuotient(lumped_mass);
    
}






template < typename Fonction, typename VECTOR, typename MATRIX , typename T >
Matrix<T, Dynamic, 1> solveFEM_Entropic( const MATRIX& mass,const VECTOR& conv, const VECTOR& last_term, const Fonction& phi_old , T dt )
{
    // if I dont consider boundary condition here, but just solve for all-> look to solveFEM
    auto b = phi_old - dt * conv + dt * last_term;
    
    return mass.completeOrthogonalDecomposition().solve(b);
    
    
}

template < typename Fonction, typename VECTOR, typename MATRIX , typename T >
Matrix<T, Dynamic, 1> solveFEM_Entropic_FAST( const MATRIX& llt,const VECTOR& conv, const VECTOR& last_term, const Fonction& phi_old , T dt )
{
    // if I dont consider boundary condition here, but just solve for all-> look to solveFEM
    auto b = phi_old - dt * conv + dt * last_term;
    
    return llt.solve(b);
    
    
}



/*

std::vector< size_t > find_neighbors(  const typename Mesh::node_type& node )
{
    size_t row = floor( node/(M+1) ) ;
    size_t col = n % (M+1) ;
    
    size_t central_cell = col + row*M ;
    
}
*/

template<typename T , typename Mesh >
class Finite_Element
{
public:
    Mesh msh; // Original mesh, not agglomerated
    // number of degree of freedom (HHO setting) , #vertices , #cells
    size_t ndof_disc , n_vertices , n_cls  ;
    size_t ndof_FE ; // number of degree of freedom (Continuous FEM setting)
    size_t num_faces = 4 ; // #faces
    mesh_init_params<T> params; // parameters mesh
    size_t order; // FE order of polynomial B^k
    size_t local_ndof ; // local degree of freedom
    
    size_t Nx , Ny ; // number of cells in x and y direction
    T hx; // step size along x direction
    T hy; // step size along y direction
    
    std::vector<std::set<size_t>> S_i;
    std::vector< std::vector<std::pair<size_t,bool>>> connectivity_matrix ;
    
    std::vector< bool > Dirichlet_boundary; //(ndof_FE , FALSE );
    
    Finite_Element(const Mesh & msh , size_t order , const mesh_init_params<T>& params ): msh(msh) , n_vertices(msh.nodes.size()) , n_cls(msh.cells.size()) ,  order(order) , local_ndof((order+1)*(order+1)), hx(params.hx() ) , hy(params.hy() ) , connectivity_matrix( n_cls , std::vector<std::pair<size_t,bool>>(local_ndof) ) , Nx(params.Nx),Ny(params.Ny) , params(params)
    {
        ndof_disc =  n_vertices*n_cls ;
        
        //std::cout<<"n_vertices "<<n_vertices<<" , n_cls "<<n_cls<<" , ndof_disc "<<ndof_disc<<std::endl;
       
        size_t counter = 0 , i_global = 0 ;
        for(const auto& cl : msh.cells)
        {
            
            std::vector<size_t> loc_bdry_nodes;
            auto fcs = faces(msh,cl);
            for(size_t face_i = 0; face_i < num_faces; face_i++)
            {
                auto fc = fcs[face_i];
                //std::cout<<"fc is "<<fc<<std::endl;
                // In HHO flow problem I can have different types of boundary;
                // For level set evolution, I just focus on Dirichelet boundary
                // The indication of the type of boundary in fc is relative to HHO problem, so useless now
                if( fc.is_boundary )
                {
                    if(face_i == 0){
                        loc_bdry_nodes.push_back(0);
                        loc_bdry_nodes.push_back(1);
                
                        if(order > 1){
                            for(size_t j = 0 ; j <= order - 2 ; j++)
                                loc_bdry_nodes.push_back(4+j);
                        }
                        
                    }
                    else if(face_i == 1){
                        loc_bdry_nodes.push_back(1);
                        loc_bdry_nodes.push_back(2);
                        if(order > 1){
                        for(size_t j = 0 ; j <= order - 2 ; j++)
                            loc_bdry_nodes.push_back(order+3+j);
                        }
                    }
                    else if(face_i == 2){
                        loc_bdry_nodes.push_back(2);
                        loc_bdry_nodes.push_back(3);
                        if(order > 1){
                        for(size_t j = 0 ; j <= order - 2 ; j++)
                            loc_bdry_nodes.push_back(2*order+2+j);
                        }
                    }
                    else if(face_i == 3){
                        loc_bdry_nodes.push_back(0);
                        loc_bdry_nodes.push_back(3);
                        if(order > 1){
                        for(size_t j = 0 ; j <= order - 2 ; j++)
                            loc_bdry_nodes.push_back(3*order+1+j);
                        }
                                       
                    }
                    else
                        exit(1);
                
                }
            
            }
            
            size_t offset_curr = offset(msh,cl);
            sort(loc_bdry_nodes.begin(), loc_bdry_nodes.end());
            loc_bdry_nodes.erase( unique( loc_bdry_nodes.begin(), loc_bdry_nodes.end() ), loc_bdry_nodes.end() );
            
           
        
            
            
            //std::cout<<"loc_bdry_nodes "<<std::endl;
            //for(auto i : loc_bdry_nodes)
            //    std::cout<<i<<" , ";
            //std::cout<<std::endl;
            
            for( size_t i_local = 0 ; i_local < local_ndof ; i_local++)
            {
                // if boundary = TRUE is node on the boudnary , else not
                bool boundary =  binary_search(loc_bdry_nodes.begin(), loc_bdry_nodes.end(), i_local) ;
                
                // Case: first cell
                if( offset_curr == 0 ){
                connectivity_matrix[counter][i_local].first = i_global;
                connectivity_matrix[counter][i_local].second = boundary;
                i_global++;
                }
                
                // Case: First row of cells, a part of cell[0]
                // Left face enumeration for cell i = Right face enumeration for cell i-1
                else if( offset_curr > 0 && offset_curr < Nx )
                {
                    if( i_local == 0  ){ // vertex 1 of cell[counter-1] = vertex 0 of cell[counter]
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-1][1].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else if(  i_local == 3  ){ // vertex 2 of cell[counter-1] = vertex 3 of cell[counter]
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-1][2].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else if(  i_local >= 3*order + 1 && i_local <= 4*order - 1   ) // nodes in face 3
                    {
                        // how many points i_local is distant from v3 vertices
                        // Idea is f1 of cell[counter-1] = f2 of cell[counter] (i.e. = means same ordering). It means:
                        //    face 1             face 3
                        // 2*ordering + 1 ---- 3*ordering + 1
                        //    |                      |
                        //    |                      |
                        //    |                      |
                        //  ordering + 3  ---- 4*ordering - 1
                        
                        size_t dist_v3 = i_local - (3*order + 1) ;
                        size_t j = 2*order + 1 - dist_v3 ;
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-1][j].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else
                    {
                        connectivity_matrix[counter][i_local].first = i_global;
                        connectivity_matrix[counter][i_local].second = boundary;
                        i_global++;
                    }
               
                }
                
                // Case: Left Boundary cells
                // Bottom face enumeration for cell i = Top face enumeration for cell i-Nx
                else if( offset_curr % Nx == 0 && offset_curr > 0 )
                {
                    if( i_local == 0  ){ // vertex 3 of cell[counter-Nx] = vertex 0 of cell[counter]
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-Nx][3].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else if(  i_local == 1   ){ // vertex 2 of cell[counter-Nx]=vertex 1 of cell[counter]
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-Nx][2].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else if( i_local >= 4 && i_local <= order + 2  ) // nodes in face 0
                    {
                        // how many points i_local is distant from v0 vertices
                        // Idea is f2 of cell[counter-Nx] = f0 of cell[counter] (i.e. = means same ordering). It means:
                        //         face 0 ; cell[counter]
                        //     4        ----     ordering + 2
                        //     |                      |
                        //         face 2 ; cell[counter-Nx]
                        //     |                      |
                        //  3*ordering  ----    2*ordering +2
                        
                        size_t dist_v0 = i_local - 4 ;
                        size_t j = 3*order - dist_v0 ;
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-Nx][j].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    
                    else
                    {
                        connectivity_matrix[counter][i_local].first = i_global;
                        connectivity_matrix[counter][i_local].second = boundary;
                        i_global++;
                    }
                }
                
                // All the other cells ( i.e. internal and right boundary cells ) are both 2 the cases
                else
                {
                    if( i_local == 0  ){ // vertex 1 of cell[counter-1] = vertex 0 of cell[counter]
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-1][1].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else if(  i_local == 3  ){ // vertex 2 of cell[counter-1] = vertex 3 of cell[counter]
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-1][2].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else if(  i_local == 1   ){ // vertex 2 of cell[counter-Nx]=vertex 1 of cell[counter]
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-Nx][2].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else if(  i_local >= 3*order + 1 && i_local <= 4*order - 1   ) // nodes in face 3
                    {
                        size_t dist_v3 = i_local - (3*order + 1) ;
                        size_t j = 2*order + 1 - dist_v3 ;
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-1][j].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else if( i_local >= 4 && i_local <= order + 2  ) // nodes in face 0
                    {
                        size_t dist_v0 = i_local - 4 ;
                        size_t j = 3*order - dist_v0 ;
                        connectivity_matrix[counter][i_local].first = connectivity_matrix[counter-Nx][j].first ;
                        connectivity_matrix[counter][i_local].second = boundary;
                    }
                    else
                    {
                        connectivity_matrix[counter][i_local].first = i_global;
                        connectivity_matrix[counter][i_local].second = boundary;
                        i_global++;
                    }
                }
                
            } // End of local loop for each cell
        counter++;
        } // End of loop over the cells.
        ndof_FE = i_global ;
        
        S_i.resize(ndof_FE);
        Dirichlet_boundary.resize(ndof_FE);
        
        
        
        for(auto&cl : msh.cells)
        {
            auto offset_cl = offset(msh,cl);
            for (size_t i = 0; i < local_ndof; i++)
            {
                auto asm_map_i = connectivity_matrix[offset_cl][i].first;
                for (size_t j = 0; j < local_ndof; j++)
                {
                    auto asm_map_j = connectivity_matrix[offset_cl][j].first;
                    S_i[asm_map_i].insert( asm_map_j );

                }
            }
        }
        
        
        /*
        size_t size_supp_nodes = 0;
        std::cout<<"Supporting nodes IN FE METHOD:"<<std::endl;
        size_t jjjj = 0;
        for (auto& i: S_i) {
            size_supp_nodes+=i.size();
            std::cout <<"Node "<<jjjj<<":";
            for (auto it=i.begin(); it != i.end(); ++it)
                std::cout << ' ' << *it;
               // std::cout<<ii;
            std::cout<<'\n';
            jjjj++;
        }
        std::cout<<std::endl;
        std::cout<<"Supporting nodes size:"<<size_supp_nodes<<std::endl;
        */
        
        //std::cout<<"Beginning of connectivity matrix"<<std::endl;
        
        
        
        // Dirichlet boundary vector
        for(size_t i = 0 ; i< n_cls ; i++)
        {
            for(size_t j = 0 ; j< local_ndof ; j++)
            {
                //std::cout<<"( "<<connectivity_matrix[i][j].first<<" , "<<connectivity_matrix[i][j].second<<" ) , ";
                auto asmap = connectivity_matrix[i][j] ;
                Dirichlet_boundary[asmap.first] = asmap.second ;
                //std::cout<<"Node "<<asmap.first << " is Dirichlet "<<asmap.second<< "  ,  ";
            }
           // std::cout<<std::endl;
        }
        //std::cout<<std::endl;
        //std::cout<<"End of connectivity matrix"<<std::endl;
        
        
        
        
    }
        
        
    //void
    //assembling(SparseMatrix<T>& Global_Mass ,SparseMatrix<T>& Global_c_term_x ,SparseMatrix<T>& Global_c_term_y, DiagonalMatrix<T, Dynamic>& Global_Mass_Lumped , Matrix<T, Dynamic, 1>& RHS){
        
        
   // }
    
    /*
    size_t get_order() const { return order; }
    size_t get_n_nodes() const { return n_nodes; }
    size_t get_n_cells() const { return n_cls; }
    size_t get_hx() const { return hx; }
    size_t get_hy() const { return hy; }
    size_t get_Nx() const { return Nx; }
    size_t get_Ny() const { return Ny; }
    size_t get_local_ndof() const { return local_ndof; }
    size_t get_ndof() const { return ndof; }
    */
};



template< typename VECTOR_LOC , typename VECTOR_GLOB , typename NODES>
void global_update_vector( const VECTOR_LOC& local , VECTOR_GLOB& global , const NODES& nodes_position )
{
    size_t i = 0 ;
   // std::cout<<"numero nodi "<<nodes_position.size()<<std::endl;
    for(auto & k : nodes_position)
    {
        size_t kk = k.ptid;
        //std::cout<<"riga "<<kk<<std::endl;
        // std::cout<<"i "<<i<<std::endl;
        global( kk ) += local ( i ) ;
        i++;
    }

}


template< typename VECTOR_LOC , typename VECTOR_GLOB , typename NODES>
void global_update_vector_NOSUM( const VECTOR_LOC& local , VECTOR_GLOB& global , const NODES& nodes_position )
{
    size_t i = 0 ;
   // std::cout<<"numero nodi "<<nodes_position.size()<<std::endl;
    for(auto & k : nodes_position)
    {
        size_t kk = k.ptid;
        //std::cout<<"riga "<<kk<<std::endl;
        // std::cout<<"i "<<i<<std::endl;
        global( kk ) = local ( i ) ;
        i++;
    }

}


template< typename VECTOR_LOC , typename VECTOR_GLOB , typename NODES>
void global_update_vector_MAX( const VECTOR_LOC& local , VECTOR_GLOB& global , const NODES& nodes_position )
{
    size_t i = 0 ;
   // std::cout<<"numero nodi "<<nodes_position.size()<<std::endl;
    for(auto & k : nodes_position)
    {
        size_t kk = k.ptid;
        //std::cout<<"riga "<<kk<<std::endl;
        // std::cout<<"i "<<i<<std::endl;
        auto max_glob = global( kk );
        global( kk ) = std::max( local ( i ) , max_glob ) ;
        i++;
    }

}

template< typename VECTOR_LOC , typename VECTOR_GLOB , typename NODES>
void global_update_vector_MIN( const VECTOR_LOC& local , VECTOR_GLOB& global , const NODES& nodes_position )
{
    size_t i = 0 ;
   // std::cout<<"numero nodi "<<nodes_position.size()<<std::endl;
    for(auto & k : nodes_position)
    {
        size_t kk = k.ptid;
        //std::cout<<"riga "<<kk<<std::endl;
        // std::cout<<"i "<<i<<std::endl;
        auto min_glob = global( kk );
        global( kk ) = std::min( local ( i ) , min_glob ) ;
        i++;
    }

}




template< typename MATRIX_LOC , typename MATRIX_GLOB , typename NODES>
void global_update_NOSUM( const MATRIX_LOC& local , MATRIX_GLOB& global , const NODES& nodes_position )
{
    size_t i = 0 ;
   // std::cout<<"numero nodi "<<nodes_position.size()<<std::endl;
    for(auto & k : nodes_position)
    {
        size_t j = 0;
        size_t kk = k.ptid;
        //std::cout<<"riga "<<kk<<std::endl;
       // std::cout<<"i "<<i<<std::endl;
        
        for(auto & l : nodes_position)
        {
            size_t ll = l.ptid;
            
            //std::cout<<"colonna "<<ll<<std::endl;
            //std::cout<<"j "<<j<<std::endl;
            global( kk  , ll  ) = local ( i , j ) ;
            j++;
        }
        i++;
        
    }
    
}



template< typename MATRIX_LOC , typename MATRIX_GLOB , typename NODES>
void global_update( const MATRIX_LOC& local , MATRIX_GLOB& global , const NODES& nodes_position )
{
    size_t i = 0 ;
   // std::cout<<"numero nodi "<<nodes_position.size()<<std::endl;
    for(auto & k : nodes_position)
    {
        size_t j = 0;
        size_t kk = k.ptid;
        //std::cout<<"riga "<<kk<<std::endl;
       // std::cout<<"i "<<i<<std::endl;
        
        for(auto & l : nodes_position)
        {
            size_t ll = l.ptid;
            
            //std::cout<<"colonna "<<ll<<std::endl;
            //std::cout<<"j "<<j<<std::endl;
            global( kk  , ll  ) += local ( i , j ) ;
            j++;
        }
        i++;
        
    }
    
}



/*

// OLD IMPLEMENTATION
template< typename MATRIX_LOC , typename MATRIX_GLOB , typename NODES>
void global_update2( const MATRIX_LOC& local , MATRIX_GLOB& global , const NODES& nodes_position )
{
    size_t i = 0 ;
    
   std::cout<<"LOCAL "<<std::endl;
   std::cout<<local<<std::endl;
    
    
    for(auto & k : nodes_position)
    {
        size_t j = 0;
        //size_t kk = k.ptid;
        std::cout<<"riga "<<k<<std::endl;
        std::cout<<"i "<<i<<std::endl;
        
        for(auto & l : nodes_position)
        {
            //size_t ll = l.ptid;
            
            std::cout<<"colonna "<<l<<std::endl;
            std::cout<<"j "<<j<<std::endl;
            
            global( k  , l  ) += local ( i , j ) ;
            j++;
        }
        i++;
        
    }
    std::cout<<"GLOBAL "<<std::endl;
    std::cout<<global<<std::endl;
    
}
*/



std::vector<size_t> boundary_nodes_function( size_t Nx , size_t Ny )
{
    // IT WORKS ONLY FOR Q1!!!
    std::vector<size_t> bdry_nodes ;
    for (size_t j = 1 ; j<=Ny ; j++) {
        bdry_nodes.push_back( j*(Nx+1) -1 ) ;
        bdry_nodes.push_back( j*(Nx+1) ) ;
    }
    for( size_t j = 1 ; j< Nx ; j++){
        bdry_nodes.push_back( j ) ;
        bdry_nodes.push_back( Ny*(Nx+1)+j ) ;
    }
    bdry_nodes.push_back( 0 ) ;
    bdry_nodes.push_back( (Ny+1)*(Nx+1)-1 ) ;
    std::sort(bdry_nodes.begin(), bdry_nodes.end() );
    
    return bdry_nodes;
}


template<typename Node , typename SUPP_NODES>
void supporting_nodes( SUPP_NODES& ret ,  const Node& nodes_position )
{
    
    for( auto & i : nodes_position)
    {
        for( auto & j : nodes_position)
        {
            (ret.at(i.ptid)).insert(j.ptid);
         //    std::cout<<"ret("<<i.ptid<<") is "<<j<<std::endl;
        }
         //std::sort(ret(i.ptid).begin(), ret(i.ptid).end() ); // set already increasing ordered
    }
   
}


template<typename Node , typename MATRIX >
void division_Si( const Node& S_i , const MATRIX& mat1 , const MATRIX& mat2 , MATRIX& sol )
{
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            sol(counter,elem) = mat1(counter,elem)/mat2(counter,elem);
        }
        counter++;
    }
}

template<typename Node , typename MATRIX >
void division_Si_T( const Node& S_i , const MATRIX& mat1 , const MATRIX& mat2 , MATRIX& sol )
{
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            sol(elem,counter) = mat1(elem,counter)/mat2(elem,counter);
        }
        counter++;
    }
}

template<typename Node , typename MATRIX >
void moltiplication_Si_T( const Node& S_i , const MATRIX& mat1 , const MATRIX& mat2 , MATRIX& sol )
{
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            sol(elem,counter) += ( mat1(elem,counter) *mat2(elem,counter) );
        }
        counter++;
    }
}


template<typename Node , typename MATRIX >
void moltiplication_Si( const Node& S_i , const MATRIX& mat1 , const MATRIX& mat2 , MATRIX& sol )
{
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            sol(counter,elem) += ( mat1(counter,elem)*mat2(counter,elem) );
        }
        counter++;
    }
}

template<typename Node , typename MATRIX , typename VECTOR>
void sum_Si( const Node& S_i , const MATRIX& mat , const VECTOR& vec , VECTOR& sol )
{
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            sol(counter) += mat(counter,elem)*vec(elem);
        }
        counter++;
    }
}

template<typename Node , typename MATRIX , typename VECTOR>
void averaged_sum_Si( const Node& S_i , const MATRIX& mat , const VECTOR& vec , VECTOR& sol )
{
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            sol(counter) += mat(counter,elem)*(vec(elem)-vec(counter));
        }
        counter++;
    }
}


/*
template < typename VECTOR, typename MATRIX , typename T , typename POS >
void checking_phi_l( const VECTOR& lumped_mass, const MATRIX& cij_x , const MATRIX& cij_y , const MATRIX& d_ij, const VECTOR& phi_old , const VECTOR& phi_l , T dt , const POS& S_i )
{
 
    
    VECTOR K_term =  Eigen::Matrix<T, Dynamic, 1>::Zero(lumped_mass.rows(), 1);
    VECTOR Dl_term =  Eigen::Matrix<T, Dynamic, 1>::Zero(lumped_mass.rows(), 1);
    // SBAGLIATO NON è PHI_OLD MA IL FLUSSO CHE SERVE QUI DENTRO!
    averaged_sum_Si( S_i , cij_x , phi_old , K_term );
    averaged_sum_Si( S_i , cij_y , phi_old , K_term );
    sum_Si( S_i , d_ij , phi_old , Dl_term );
    //std::cout<<'\n'<<"check Dij "<<'\n'<<Dl_term-d_ij*phi_old<<std::endl;
    VECTOR tmp =  lumped_mass.cwiseProduct( (phi_l - phi_old)/dt ) + K_term - Dl_term;
    std::cout<<'\n'<<"check phi_l vec is "<<'\n'<<tmp<<std::endl;
}
*/
template < typename VECTOR , typename T  >
void checking_phi_lBIS( const VECTOR& lumped_mass, const VECTOR& cij_term , const VECTOR& dij_term, const VECTOR& phi_old , const VECTOR& phi_l , T dt )
{
 
    VECTOR tmp =  lumped_mass.cwiseProduct( (phi_l - phi_old)/dt ) + cij_term - dij_term;
    std::cout<<'\n'<<"CHECK BIS phi_l vec is "<<'\n'<<tmp<<std::endl;
}


/*
template < typename VECTOR , typename MATRIX , typename T , typename POS >
void checking_phi_h( const MATRIX& mass, const MATRIX& cij_x , const MATRIX& cij_y , const MATRIX& dC_ij , const VECTOR& phi_old , const VECTOR& phi_h , T dt , const POS& S_i)
{
 
    VECTOR K_term =  Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    VECTOR Dc_term =  Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    // SBAGLIATO NON è PHI_OLD MA IL FLUSSO CHE SERVE QUI DENTRO!
    averaged_sum_Si( S_i , cij_x , phi_old , K_term );
    averaged_sum_Si( S_i , cij_y , phi_old , K_term );
    sum_Si( S_i , dC_ij , phi_old , Dc_term );
    VECTOR tmp =  mass*( (phi_h - phi_old)/dt ) + K_term  - dC_ij*phi_old;
   // std::cout<<'\n'<<"check DCij "<<'\n'<<Dc_term-dC_ij*phi_old<<std::endl;
    std::cout<<'\n'<<"check phi_h vec is "<<'\n'<<tmp<<std::endl;
}
*/
template < typename VECTOR , typename T , typename MATRIX  >
void checking_phi_hBIS( const MATRIX& mass, const VECTOR& cij_term , const VECTOR& dij_term, const VECTOR& phi_old , const VECTOR& phi_h , T dt )
{
 
    VECTOR tmp =  mass*( (phi_h - phi_old)/dt ) + cij_term - dij_term;
    std::cout<<'\n'<<"CHECK BIS phi_H vec is "<<'\n'<<tmp<<std::endl;
}




template < typename Entropy , typename Fonction, typename Fonction_TILDE , typename Mesh, typename Vel_Field , typename T , typename POSITION , typename VECTOR >
  void r_i_calculator(const Mesh& msh, const typename Mesh::cell_type& cl , const Entropy& E , const Fonction_TILDE& phi_tilde , const Fonction& phi , T dt , const Vel_Field& u , const POSITION& S_i , VECTOR& Emax_global , VECTOR& Emin_global , VECTOR& R_i  )
{
    size_t di = 0;
    cell_basis_Lagrangian<Mesh,T> cb(msh, cl, phi.degree_FEM);
    auto cbs = cb.size();
    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);
    auto qps = integrate(msh, cl, 2*(phi.degree_FEM+di));
    
    
    auto nds = nodes(msh,cl);
    iter_swap(nds.end()-1 , nds.end()-2 );
    
    auto pts = points(msh,cl);
    iter_swap(pts.end()-1 , pts.end()-2 );
   
    
    
    size_t counter = 0;
    /// FARE FUNZIONE max_min_entropy PER RENDERE IL CODICE PULITO
    Eigen::Matrix<T, Dynamic, 1> Emax =  Eigen::Matrix<T, Dynamic, 1>::Ones(pts.size(), 1);
    Eigen::Matrix<T, Dynamic, 1> Emin =  Eigen::Matrix<T, Dynamic, 1>::Ones(pts.size(), 1);
    for(auto& nd : nds )
    {
        auto i = nd.ptid ;
        T max_loc = -1e20;
        T min_loc = 1e20;
        for(auto& pt_j : pts )
        {
            Emax(counter) = std::max ( E(pt_j , cl) , max_loc );
            Emin(counter) = std::min ( E(pt_j , cl) , min_loc );
            max_loc = std::max ( Emax(counter) , max_loc );
            min_loc = std::min ( Emin(counter) , min_loc );
        }
        
        counter++;
    }
    global_update_vector_MAX( Emax , Emax_global , nds );
    global_update_vector_MIN( Emin , Emin_global , nds );
    // std::cout<<"Emax_global: "<<'\n'<<Emax_global<<std::endl;
    //  std::cout<<"Emin_global: "<<'\n'<<Emin_global<<std::endl;
    
    for (auto& qp : qps)
    {
        auto bi = cb.eval_basis(qp.first);
        auto phi_grad0 = phi.gradient(qp.first,msh,cl)(0);
        auto phi_grad1 = phi.gradient(qp.first,msh,cl)(1);
        
        auto f = ( ( ( phi_tilde(qp.first,msh,cl)-phi(qp.first,msh,cl) )/dt + u(qp.first,msh,cl).first*phi_grad0 + u(qp.first,msh,cl).second*phi_grad1 ) *
                  E.derivative(qp.first,cl) ); //.cwiseQuotient( Emax - Emin ) ;
        ret += qp.second * bi * f;
    }
    global_update_vector( ret , R_i , nds );
    
}



template < typename T , typename MATRIX , typename VECTOR ,typename POSITION >
VECTOR alfaf_ij_creator( const VECTOR& lumped_mass , const MATRIX& mass , const VECTOR& delta_phi , const VECTOR& phi_L , T dt , const MATRIX& D_ij , const MATRIX& Dc_ij , const VECTOR& phi_old , const POSITION& S_i)
{
    
    MATRIX new_mass = - mass;
    for(size_t i = 0 ; i < new_mass.rows() ; i++){
        new_mass(i,i) += lumped_mass(i);
    }
    
    MATRIX new_D = (Dc_ij - D_ij); //(D_ij - Dc_ij);
    
    MATRIX f_ij = Eigen::Matrix<T, Dynamic, Dynamic>::Zero( mass.rows(), mass.cols() );
    MATRIX alpha_ij = Eigen::Matrix<T, Dynamic, Dynamic>::Zero( mass.rows(), mass.cols() );
    
    VECTOR ret = Eigen::Matrix<T, Dynamic, 1>::Zero( mass.rows(), 1 );
    
    VECTOR P_plus = Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    VECTOR P_minus = Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    VECTOR Q_plus = Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    VECTOR Q_minus = Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    
    VECTOR R_plus = Eigen::Matrix<T, Dynamic, 1>::Ones(mass.rows(), 1);
    VECTOR R_minus = Eigen::Matrix<T, Dynamic, 1>::Ones(mass.rows(), 1);
    VECTOR phi_max = phi_old;
    VECTOR phi_min = phi_old;
    
               
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            phi_max(counter) = std::max(phi_old(elem),phi_max(counter));
            phi_min(counter) = std::min(phi_old(elem),phi_min(counter));
        }
        counter++;
    }
    
    for(size_t i = 0; i< f_ij.rows() ; i++)
    {
        for(size_t j = 0; j< f_ij.cols() ; j++)
        {
            f_ij(i, j) = new_mass(i,j)*(delta_phi(j)-delta_phi(i)) +dt*new_D(i,j)*(phi_old(j)-phi_old(i)) ;
            T f_tmp = f_ij(i, j);
            P_plus(i) += std::max( 0. , f_tmp );
            P_minus(i) += std::min( 0. , f_tmp );
        }
        Q_plus(i) = lumped_mass(i)*(phi_max(i)-phi_L(i));
        Q_minus(i) = lumped_mass(i)*(phi_min(i)-phi_L(i));
        if( std::abs(P_plus(i)) > 1e-20 ){
            T Q_P_plus = Q_plus(i)/P_plus(i);
            R_plus(i) = std::min( 1. , Q_P_plus );
        }
        if( std::abs(P_minus(i)) > 1e-20 ){
            T Q_P_minus = Q_minus(i)/P_minus(i);
            R_minus(i) = std::min( 1. , Q_P_minus );
        }
    }
    
    for(size_t i = 0; i< f_ij.rows() ; i++)
    {
        for(size_t j = 0; j< f_ij.cols() ; j++)
        {
            if( f_ij(i, j) > 0 )
                alpha_ij(i,j) = std::min( R_plus(i) , R_minus(j) );
            else
                alpha_ij(i,j) = std::min( R_plus(j) , R_minus(i) );
        }
       
    }
    
    size_t counter2 = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            ret(counter2) += ( alpha_ij(counter2,elem)*f_ij(counter2,elem) );
        }
        counter2++;
    }

    
    // CHECKING F and Alpha properties
    /*
    for(size_t i = 0; i< f_ij.rows() ; i++)
    {
        for(size_t j = 0; j< f_ij.cols() ; j++)
        {
            std::cout<<f_ij(i,j)+f_ij(j,i)<<" , ";
        }
       std::cout<<'\n'<<std::endl;
    }
    */
    
    /*
    std::cout<<"alpha_ij checking symmetry"<<'\n'<<std::endl;
    for(size_t i = 0; i< f_ij.rows() ; i++)
    {
        for(size_t j = 0; j< f_ij.cols() ; j++)
        {
            std::cout<<alpha_ij(i,j)-alpha_ij(j,i)<<" , ";
        }
       std::cout<<'\n'<<std::endl;
    }
    */
    
    /*
    size_t counter3 = 0;
    
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            std::cout<<f_ij(counter3,elem)+f_ij(elem,counter3)<<" , ";
        }
        counter3++;
        std::cout<<'\n';
    }
    std::cout<<std::endl;
     
    counter3 = 0;
     std::cout<<"alpha_ij checking symmetry: METODO 2"<<'\n'<<std::endl;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            std::cout<<alpha_ij(counter3,elem)-alpha_ij(elem,counter3)<<" , ";
        }
        counter3++;
        std::cout<<'\n';
    }
    std::cout<<std::endl;
    */
    //  CHECKING F_IJ
    VECTOR ret0 = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    VECTOR ret1 = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    averaged_sum_Si( S_i , new_mass , delta_phi , ret0 );
    averaged_sum_Si( S_i , new_D , phi_old , ret1 );
    VECTOR f_i = (ret0 + dt*ret1);
    
    VECTOR f_i_NEW = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    VECTOR One_Vec = Eigen::Matrix<T, Dynamic, 1>::Ones(new_mass.rows(), 1);
    sum_Si( S_i , f_ij , One_Vec , f_i_NEW );
    
    //std::cout<<'\n'<<"Diff f_i 2 meths is:"<<'\n'<< f_i_NEW - f_i<<std::endl;
    
    return  ret ;
    
}




template < typename T , typename MATRIX , typename VECTOR ,typename POSITION >
VECTOR f_ij_creator( const VECTOR& lumped_mass , const MATRIX& mass , const VECTOR& delta_phi , T dt , const MATRIX& D_ij , const MATRIX& Dc_ij , const VECTOR& phi_old , const POSITION& S_i)
{
    
    MATRIX new_mass = - mass;
    MATRIX new_D = (Dc_ij - D_ij); //(D_ij - Dc_ij);
    for(size_t i = 0 ; i < new_mass.rows() ; i++){
        new_mass(i,i) += lumped_mass(i);
    }
    
    /*
    std::cout<<"new_mass sum check: "<<std::endl;
    for(size_t i = 0 ; i < new_mass.rows() ; i++){
        new_mass(i,i) += lumped_mass(i);
        std::cout<<(new_mass.row(i)).sum()<<std::endl;
    }
    std::cout<<"new_D check: "<<std::endl;
    for(size_t i = 0 ; i < new_D.rows() ; i++){
       std::cout<<(new_D.row(i).sum() )<<std::endl;
    }
    */
    
    VECTOR ret0 = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    VECTOR ret1 = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    averaged_sum_Si( S_i , new_mass , delta_phi , ret0 );
    averaged_sum_Si( S_i , new_D , phi_old , ret1 );
    /*
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            sol(counter) += mat(counter,elem)*(vec(elem)-vec(counter));
        }
        counter++;
    }
    */
    //  std::cout<<'\n'<<"dt: "<<'\n'<<dt<<std::endl;
   
    //    std::cout<<'\n'<<"ret0: "<<'\n'<<ret0<<std::endl;
    //    std::cout<<'\n'<<"ret1: "<<'\n'<<dt*ret1<<std::endl;
  
    return (ret0 + dt*ret1) ;
    
}
template<typename FUNCTION, typename T>
void mapping_phi(FUNCTION& phi , T phi_max , T phi_min)
{
    Matrix<T, Dynamic, 1> Vec_One = Eigen::Matrix<T, Dynamic, 1>::Ones( phi.vertices.rows(), 1 );
    Matrix<T, Dynamic, Dynamic> Mat_One = Eigen::Matrix<T, Dynamic, Dynamic>::Ones( phi.values_bis.rows(), phi.values_bis.cols() );
    // mapping between 0 and 1
    phi.vertices = (phi.vertices-phi_min*Vec_One)/( phi_max - phi_min );
    phi.values_bis = (phi.values_bis-phi_min*Mat_One)/( phi_max - phi_min );
}

template<typename FUNCTION, typename T>
void inverse_mapping_phi(FUNCTION& phi, T phi_max , T phi_min)
{
    Matrix<T, Dynamic, 1> Vec_One = Eigen::Matrix<T, Dynamic, 1>::Ones( phi.vertices.rows() , 1 );
    Matrix<T, Dynamic, Dynamic> Mat_One = Eigen::Matrix<T, Dynamic, Dynamic>::Ones( phi.values_bis.rows(), phi.values_bis.cols() );
    phi.vertices = phi_min*Vec_One + phi.vertices*( phi_max - phi_min );
    phi.values_bis = phi_min*Mat_One + phi.values_bis*( phi_max - phi_min );
}



template<typename T, typename Mesh , typename FiniteSpace>
struct L2_projection
{
    SparseMatrix<T>                 Global_Mass; // Global mass, saved for FEM problem
    Matrix<T, Dynamic, 1>           RHS;    // Known term
    std::vector< Triplet<T> >       triplets; // Position elements: Sparse Matrix Notation
    
    Eigen::Matrix<T, Dynamic, Dynamic> sol_HHO ; // projection saved in HHO format: cell by cell
    Matrix<T, Dynamic, 1> sol_FEM ; // projection saved in Continuos FE format: global nodes
    Eigen::Matrix<T, Dynamic, 1> vertices; // saving level_set on vertices mesh


    size_t number_elements;
    
    size_t n_cls ; // #cells
    size_t local_dim; // Local Dimension (degree_FEM+1)*(degree_FEM+1)
    size_t n_vertices ; // #vertices

    size_t degree_FEM;
    Mesh msh;
    size_t      Nx, Ny ;
    mesh_init_params<T> params;
    
    std::vector< std::vector<std::pair<size_t,bool>>> connectivity_matrix ;
  
    size_t dim_HHO; // Global dimension Discontinuous framework = Local dimension * #cells
    size_t ndof_FE; // Global dimension FE continuous = #nodes
       

    L2_projection(const FiniteSpace& fe_data ,  const Mesh & msh)
           : degree_FEM(fe_data.order) , local_dim(fe_data.local_ndof), msh(msh), Nx(fe_data.Nx),Ny(fe_data.Ny), params(fe_data.params) , dim_HHO(fe_data.ndof_disc) , n_cls(fe_data.n_cls) ,n_vertices(fe_data.n_vertices) , connectivity_matrix(fe_data.connectivity_matrix) , ndof_FE(fe_data.ndof_FE)
       {
           sol_HHO = Eigen::Matrix<T, Dynamic, Dynamic>::Zero( local_dim, n_cls );
           sol_FEM = Matrix<T, Dynamic, 1>::Zero(ndof_FE);
           vertices = Eigen::Matrix<T, Dynamic, 1>::Zero( n_vertices , 1 );
       }
    

    L2_projection()=default;
    
    
    void  set_discrete_points( Eigen::Matrix<T, Dynamic, Dynamic>& values_new)
    {
        std::cout<<"L2 projection for phi_tilde post-resolution."<<std::endl;
        sol_HHO = values_new;
    }
    
    
    void converting_into_HHO_formulation( const Eigen::Matrix<T, Dynamic, 1>& values_new )
    {
        std::cout<<"L2 projection for phi_tilde post-resolution."<<std::endl;
        for(size_t counter_bis = 0 ; counter_bis < n_cls ;counter_bis++)
        {
            for (size_t i = 0; i < local_dim; i++){
                size_t asm_map =  connectivity_matrix[counter_bis][i].first ;
                sol_HHO(i,counter_bis) = values_new( asm_map );
            }
            size_t i_vertex = counter_bis+floor(counter_bis/Nx);
            vertices(i_vertex) = sol_HHO(0,counter_bis) ;
            vertices(i_vertex+1) = sol_HHO(1,counter_bis) ;
            vertices(i_vertex+Nx+2) = sol_HHO(2,counter_bis) ;
            vertices(i_vertex+Nx+1) = sol_HHO(3,counter_bis) ;
        }
                 
    }
    
    void converting_into_FE_formulation( const Eigen::Matrix<T, Dynamic, Dynamic>& values_new )
    {
        std::cout<<"L2 projection for phi_tilde post-resolution."<<std::endl;
        for(size_t counter_bis = 0 ; counter_bis < n_cls ;counter_bis++)
        {
            for (size_t i = 0; i < local_dim; i++){
                size_t asm_map =  connectivity_matrix[counter_bis][i].first ;
                sol_FEM( asm_map ) = values_new(i,counter_bis) ;
            }
        }

                 
    }
       
    
    
    T operator()( const typename Mesh::node_type& node ) const
    {
        return vertices(node.ptid);
    }
    
    T operator()( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl ) const
    {
        
        size_t counter = offset(msh,cl) ;
        cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
        //cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
        auto values_cell = (sol_HHO.block(0,counter,local_dim,1)).col(0);
        T tmp = values_cell.dot( cb.eval_basis(pt) );
        return tmp;
    }
    
    Eigen::Matrix<T,2,1> gradient( const point<T,2>& pt, const Mesh & msh,  const typename Mesh::cell_type& cl ) const
       {
       
           // MATRIX NOTATION
           size_t counter = offset(msh,cl);
           Eigen::Matrix<T,2,1> ret = Matrix<T, 2, 1>::Zero(2, 1);
           cell_basis_Bernstein <Mesh,T> cb(msh, cl, degree_FEM);
           //cell_basis_Lagrangian<Mesh,T> cb(msh, cl, degree_FEM);
           auto values_cell = sol_HHO.col(counter);
           auto grad_eval =  cb.eval_gradients(pt);
           ret(0) = values_cell.dot( grad_eval.col(0) );
          // std::cout<<"Value of derivative new along x"<<ret(0)<<std::endl;
           ret(1) = values_cell.dot( grad_eval.col(1) );
           //values_cell.dot( grad_eval.col(1) );
          // std::cout<<"Value of derivative new along y"<<ret(1)<<std::endl;
           return ret;
           
       }

    
       // IT WORKS FOR NOT-AGGLOMERATED MESHES --> FAST
       Eigen::Matrix<T,2,1> normal(const point<T,2>& pt, const Mesh& msh, const typename Mesh::cell_type& cl) const
       {
           Eigen::Matrix<T,2,1> ret;
           ret = gradient(pt,msh,cl);
           return ret/ret.norm();
           
       }
       
    
};

template< typename VEC , typename DIAG , typename T >
T time_step_CFL_L2_velocity_NEW( const DIAG& dii , const VEC& lumped_mass , const std::vector< bool >& Dirichlet_boundary , T& dt )
{
    T eps = 0.0;
    T tau_min  = dt;
    T CFL_numb = 10.0;
    size_t i = 0 ;
    for (const auto& dir_elem : Dirichlet_boundary )
    {
        if(!dir_elem){
            if( (1 + 2*tau_min*dii.coeff(i)/lumped_mass(i) ) < eps ){
                std::cout<<"tau_min PRE modification is "<<tau_min;
                tau_min = (eps-1.0)/2.0*lumped_mass(i)/dii.coeff(i);
                std::cout<<" and tau_min POST modification is "<<tau_min<<std::endl;
            }
            CFL_numb = std::min( CFL_numb ,(eps-1.0)/2.0*lumped_mass(i)/dii.coeff(i) );
        }
        i++ ;
    }
    std::cout<<"CFL_numb ---> "<<CFL_numb<<std::endl;
    std::cout<<yellow<<bold<<"dt is "<<dt<<" and tau min is "<<tau_min<<reset<<std::endl;
    dt = std::min(dt,tau_min);
    return CFL_numb;
}



template < typename Entropy , typename Fonction, typename Fonction_TILDE , typename Mesh, typename Vel_Field , typename T , typename VECTOR >
  void r_i_calculator_Bernstein(const Mesh& msh, const typename Mesh::cell_type& cl , const Entropy& E , const Fonction_TILDE& phi_tilde , const Fonction& phi , T dt , const Vel_Field& u , VECTOR& Emax_global , VECTOR& Emin_global , VECTOR& R_i  )
{
    size_t di = 1;
    cell_basis_Bernstein<Mesh,T> cb(msh, cl, phi.degree_FEM);
    auto cbs = cb.size();
    Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs);
    auto qps = integrate(msh, cl, 2*(phi.degree_FEM+di));
    
    auto nds = equidistriduted_nodes_ordered_bis<T,Mesh>(msh, cl, phi.degree_FEM );
    
    size_t offset_cell = offset( msh , cl ) ;
    
    T max_loc = -1e20;
    T min_loc = 1e20;
    
    for(auto& ndj : nds )
    {
        max_loc = std::max ( E(ndj , cl) , max_loc );
        min_loc = std::min ( E(ndj , cl) , min_loc );
    }
    
    for (size_t i = 0; i < phi.local_dim; i++)
    {
        size_t asm_map =  phi.connectivity_matrix[offset_cell][i].first ;
        Emax_global(asm_map) = std::max(Emax_global(asm_map) , max_loc );
        Emin_global(asm_map) = std::min(Emin_global(asm_map) , min_loc ) ;
    }
    
    //std::cout<<"Emax_global: "<<'\n'<<Emax_global<<std::endl;
    //std::cout<<"Emin_global: "<<'\n'<<Emin_global<<std::endl;
    
    
    for (auto& qp : qps)
    {
        auto bi = cb.eval_basis(qp.first);
        auto phi_grad0 = phi.gradient(qp.first,msh,cl)(0);
        auto phi_grad1 = phi.gradient(qp.first,msh,cl)(1);
        
        auto f = ( ( ( phi_tilde(qp.first,msh,cl)-phi(qp.first,msh,cl) )/dt + u(qp.first,msh,cl).first * phi_grad0 + u(qp.first,msh,cl).second * phi_grad1 ) *
                  E.derivative(qp.first,cl) ); //.cwiseQuotient( Emax - Emin ) ;
        ret += qp.second * bi * f;
    }
    
    for (size_t i = 0; i < phi.local_dim; i++)
    {
        size_t asm_map =  phi.connectivity_matrix[offset_cell][i].first ;
        R_i(asm_map) += ret(i) ;
    }
   
    
}


template <typename T>
void positive_part_SPARSE( SparseMatrix<T>& mat) {
    for (size_t i = 0; i<mat.rows();i++) {
        for (size_t j = 0; j<mat.cols();j++) {
            if( mat.coeff(i,j) < 0. )
                mat.coeffRef(i,j) = 0.;
        }
    }

}


template < typename T , typename MATRIX , typename VECTOR ,typename POSITION >
VECTOR alfaf_ij_creator_SPARSE( const VECTOR& lumped_mass , const MATRIX& mass , const VECTOR& delta_phi , const VECTOR& phi_L , T dt , const MATRIX& D_ij , const MATRIX& Dc_ij , const VECTOR& phi_old , const POSITION& S_i)
{
    
    MATRIX new_mass = - mass;
    for(size_t i = 0 ; i < new_mass.rows() ; i++){
        new_mass.coeffRef(i,i) += lumped_mass(i);
    }
    
    MATRIX new_D = (Dc_ij - D_ij); //(D_ij - Dc_ij);
    
    SparseMatrix<T> f_ij =  SparseMatrix<T>( mass.rows() , mass.cols() );
    std::vector< Triplet<T> >   triplets0;
    SparseMatrix<T> alpha_ij =  SparseMatrix<T>( mass.rows() , mass.cols() );
    std::vector< Triplet<T> >   triplets1;
    
    
    //MATRIX f_ij = Eigen::Matrix<T, Dynamic, Dynamic>::Zero( mass.rows(), mass.cols() );
    //MATRIX alpha_ij = Eigen::Matrix<T, Dynamic, Dynamic>::Zero( mass.rows(), mass.cols() );
    
    VECTOR ret = Eigen::Matrix<T, Dynamic, 1>::Zero( mass.rows(), 1 );
    
    VECTOR P_plus = Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    VECTOR P_minus = Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    VECTOR Q_plus = Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    VECTOR Q_minus = Eigen::Matrix<T, Dynamic, 1>::Zero(mass.rows(), 1);
    
    VECTOR R_plus = Eigen::Matrix<T, Dynamic, 1>::Ones(mass.rows(), 1);
    VECTOR R_minus = Eigen::Matrix<T, Dynamic, 1>::Ones(mass.rows(), 1);
    VECTOR phi_max = phi_old;
    VECTOR phi_min = phi_old;
    
               
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            phi_max(counter) = std::max(phi_old(elem),phi_max(counter));
            phi_min(counter) = std::min(phi_old(elem),phi_min(counter));
        }
        counter++;
    }
    
    counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            auto value = new_mass.coeff(counter,elem)*(delta_phi(elem)-delta_phi(counter)) + dt*new_D.coeff(counter,elem)*(phi_old(elem)-phi_old(counter)) ;
            
            triplets0.push_back( Triplet<T>(counter, elem , value ));
            P_plus(counter) += std::max( 0. , value );
            P_minus(counter) += std::min( 0. , value );
                                
        }
        
        Q_plus(counter) = lumped_mass(counter)*(phi_max(counter)-phi_L(counter));
        Q_minus(counter) = lumped_mass(counter)*(phi_min(counter)-phi_L(counter));
        if( std::abs(P_plus(counter)) > 1e-20 ){
            T Q_P_plus = Q_plus(counter)/P_plus(counter);
            R_plus(counter) = std::min( 1.0 , Q_P_plus );
        }
        if( std::abs(P_minus(counter)) > 1e-20 ){
            T Q_P_minus = Q_minus(counter)/P_minus(counter);
            R_minus(counter) = std::min( 1.0 , Q_P_minus );
        }
        
        
        counter++;
    }
    

    f_ij.setFromTriplets( triplets0.begin(), triplets0.end() );
    triplets0.clear();
    
    size_t i = 0 ;
    for(auto& row_i:S_i)
    {
        for(auto& j:row_i)
        {
            if( f_ij.coeff(i, j) > 0 ){
                auto value = std::min( R_plus(i) , R_minus(j) );
                triplets1.push_back( Triplet<T>(i, j , value ));
                
            }
            else{
                auto value = std::min( R_plus(j) , R_minus(i) );
                triplets1.push_back( Triplet<T>(i, j , value ));
                
            }
        }
        i++;
    }
    
    alpha_ij.setFromTriplets( triplets1.begin(), triplets1.end() );
    triplets1.clear();
    
    
    size_t counter2 = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            ret(counter2) += ( alpha_ij.coeff(counter2,elem)*f_ij.coeff(counter2,elem) );
        }
        counter2++;
    }

    
    // CHECKING F and Alpha properties
    
    /*
    std::cout<<"f_ij checking symmetry"<<'\n'<<std::endl;
    for(size_t i = 0; i< f_ij.rows() ; i++)
    {
        for(size_t j = 0; j< f_ij.cols() ; j++)
        {
            std::cout<<f_ij.coeff(i,j)+f_ij.coeff(j,i)<<" , ";
        }
       std::cout<<'\n'<<std::endl;
    }
    
    
    
    std::cout<<"alpha_ij checking symmetry"<<'\n'<<std::endl;
    for(size_t i = 0; i< f_ij.rows() ; i++)
    {
        for(size_t j = 0; j< f_ij.cols() ; j++)
        {
            std::cout<<alpha_ij.coeff(i,j)-alpha_ij.coeff(j,i)<<" , ";
        }
       std::cout<<'\n'<<std::endl;
    }
    */
    
    
    
    /*
    std::cout<<"f_ij checking symmetry"<<'\n'<<std::endl;
    size_t counter3 = 0;
    
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            std::cout<<f_ij.coeff(counter3,elem)+f_ij.coeff(elem,counter3)<<" , ";
        }
        counter3++;
        std::cout<<'\n';
    }
    std::cout<<std::endl;
     
    counter3 = 0;
     std::cout<<"alpha_ij checking symmetry: METODO 2"<<'\n'<<std::endl;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            std::cout<<alpha_ij.coeff(counter3,elem)-alpha_ij.coeff(elem,counter3)<<" , ";
        }
        counter3++;
        std::cout<<'\n';
    }
    std::cout<<std::endl;
    
    //  CHECKING F_IJ
    
    VECTOR ret0 = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    VECTOR ret1 = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    size_t counter4 = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            ret0(counter4) += new_mass.coeff(counter4,elem)*(delta_phi(elem)-delta_phi(counter4));
            ret1(counter4) += new_D.coeff(counter4,elem)*(phi_old(elem)-phi_old(counter4));
        }
        counter++;
    }
    
    VECTOR f_i = (ret0 + dt*ret1);
    
    VECTOR f_i_NEW = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    VECTOR One_Vec = Eigen::Matrix<T, Dynamic, 1>::Ones(new_mass.rows(), 1);
    //sum_Si( S_i , f_ij , One_Vec , f_i_NEW );
    size_t counter5 = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            f_i_NEW(counter5) += f_ij.coeff(counter5,elem)*(One_Vec(elem));
        }
        counter5++;
    }
    
    std::cout<<'\n'<<"Diff f_i 2 meths is:"<<'\n'<< f_i_NEW - f_i<<std::endl;
    */
    return  ret ;
    
}



template < typename T , typename MATRIX , typename VECTOR ,typename POSITION >
VECTOR f_ij_creator_SPARSE( const VECTOR& lumped_mass , const MATRIX& mass , const VECTOR& delta_phi , T dt , const MATRIX& D_ij , const MATRIX& Dc_ij , const VECTOR& phi_old , const POSITION& S_i)
{
    
    MATRIX new_mass = - mass;
    MATRIX new_D = (Dc_ij - D_ij); //(D_ij - Dc_ij);
    for(size_t i = 0 ; i < new_mass.rows() ; i++){
        new_mass.coeffRef(i,i) += lumped_mass(i);
    }
    
    /*
    std::cout<<"new_mass sum check: "<<std::endl;
    for(size_t i = 0 ; i < new_mass.rows() ; i++){
        new_mass(i,i) += lumped_mass(i);
        std::cout<<(new_mass.row(i)).sum()<<std::endl;
    }
    std::cout<<"new_D check: "<<std::endl;
    for(size_t i = 0 ; i < new_D.rows() ; i++){
       std::cout<<(new_D.row(i).sum() )<<std::endl;
    }
    */
    
    VECTOR ret0 = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    VECTOR ret1 = Eigen::Matrix<T, Dynamic, 1>::Zero(new_mass.rows(), 1);
    //averaged_sum_Si( S_i , new_mass , delta_phi , ret0 );
    //averaged_sum_Si( S_i , new_D , phi_old , ret1 );
    
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            ret0(counter) += new_mass.coeff(counter,elem)*(delta_phi(elem)-delta_phi(counter));
            ret1(counter) += new_D.coeff(counter,elem)*(phi_old(elem)-phi_old(counter));
            
            
        }
        counter++;
    }
    
    //  std::cout<<'\n'<<"dt: "<<'\n'<<dt<<std::endl;
   
    //    std::cout<<'\n'<<"ret0: "<<'\n'<<ret0<<std::endl;
    //    std::cout<<'\n'<<"ret1: "<<'\n'<<dt*ret1<<std::endl;
  
    return (ret0 + dt*ret1) ;
    
}








template < typename Fonction, typename Mesh, typename Vel_Field , typename FiniteSpace , typename T = typename Mesh::coordinate_type >
void
run_FEM_BERNSTEIN_CORRECT(const Mesh & msh, const FiniteSpace& fe_data, Fonction & phi , Vel_Field& u , T& dt )
{
    // Starting time for FE calculation
    std::cout<<yellow<<bold<<"----------- STARTING TRANSPORT PROBLEM -----------"<<reset<<std::endl;
    //std::cout<<yellow<<bold<<"PROVA--- USO MIXC LAGRANGE- BERNSTEIN"<<reset<<std::endl;
    
    timecounter tc;
    tc.tic();
    
    size_t degree = fe_data.order; // finite element order
    size_t dim = fe_data.ndof_FE ;
    size_t n_cls = fe_data.n_cls ;
    size_t local_ndof = fe_data.local_ndof ; // local degrees of freedom
    auto S_i = fe_data.S_i;
    
    // TOLTO PER ORA IL MAPPING NON SERVE NEL CASO BASE
    
    phi.coefficients_mapping(); // mapping of phi to have a phi between 0 and 1
   
    
    /// PLOT ---> PHI MAPPED INTO 0-1
    postprocess_output<double> postoutput4;
    auto test_phi_mapped = std::make_shared< gnuplot_output_object<double> >("phi_mapped.dat");
    for (auto cl:msh.cells) {
        auto pts = points(msh,cl) ;
        for (auto pt:pts) {
            T value = phi(pt,msh,cl);
            test_phi_mapped->add_data(pt,value);
        }
    }
    postoutput4.add_object(test_phi_mapped);
    postoutput4.write();
    
    
    // SAVING PHI AND VELOCITY COEFFS
    auto phi_FEM = phi.sol_FEM ;
    auto u0 = u.sol_FEM.first ;
    auto u1 = u.sol_FEM.second ;
    //u.set_max_vel(); // UPLOADING MAX VELOCITY OF u -> Bernstein Basis
    
    //std::cout<<"CONSTANT VELOCITY: u0 is "<<u0<<" and u1 is "<<u1<<std::endl;
    
    // NON LINEAR ENTROPY INITIALISATION
    const T eps = 1e-14 ; //constant into entropy
    non_linear_entropy<T,Fonction,Mesh> E(eps , phi ,msh );
    //SONOQUA1
    // PHI TILDE INITIALISATION --> (FOR HIGH ORDER METHOD)
    auto phi_tilde = L2_projection< T, Mesh , FiniteSpace> ( fe_data , msh );
    
    //L2_projection< T, Mesh , FiniteSpace> ( fe_data , msh );
    
    /// ENTROPY PLOT
    /*
    postprocess_output<double> postoutput6;
    auto cut_entropy = std::make_shared< gnuplot_output_object<double> >("entropy.dat");
    for (auto& cl:msh.cells)
    {
        auto pts = equidistriduted_nodes_ordered_bis<T,Mesh>( msh, cl, degree );
        for (auto& pt:pts)
        {
            T value = E(pt,cl);
            cut_entropy->add_data(pt,value);
        }
    }
    postoutput6.add_object(cut_entropy);
    postoutput6.write();
    */
    
   
    
    // SAVING OF USEFUL MATRICES
    auto global_mass = phi.Global_Mass ;
    Matrix<T, Dynamic, 1> global_lumped_mass = phi.Global_Mass_Lumped;
    
    auto global_cij_x = phi.Global_c_term_x ;
    auto global_cij_y = phi.Global_c_term_y ;
    auto local_vandermonde = phi.local_vandermonde ;
    
    // THE GLOBAL MATRIX c_ij varies each time: dependes on u -> c_ij = int(u*b_i*grad(b_j))
    /*
    SparseMatrix<T>                 global_cij_x = SparseMatrix<T>( dim, dim );
    SparseMatrix<T>                 global_cij_y = SparseMatrix<T>( dim, dim );
    std::vector< Triplet<T> >       triplets_c_term_x;
    std::vector< Triplet<T> >       triplets_c_term_y;
    for( const auto& cl : msh.cells )
    {
        size_t cell_offset = offset(msh, cl) ;
    
        
        // Local c_ij = u(x)*b_i* grad(b_j)
        auto local_cij = make_bernstein_local_cij_matrix_with_velocity(msh, cl, degree, u );
        
        for (size_t i = 0; i < local_ndof; i++)
        {
            auto asm_map_i = phi.connectivity_matrix[cell_offset][i].first ;
            for (size_t j = 0; j < local_ndof; j++)
            {
                auto asm_map_j = phi.connectivity_matrix[cell_offset][j].first ;
                triplets_c_term_x.push_back( Triplet<T>(asm_map_i, asm_map_j, local_cij.first(i,j) ) );
                triplets_c_term_y.push_back( Triplet<T>(asm_map_i, asm_map_j, local_cij.second(i,j) ) );
           }
            
        }
    }
    // FINALISATION c_ij GLOBAL ASSEMBLING
    global_cij_x.setFromTriplets( triplets_c_term_x.begin(), triplets_c_term_x.end() );
    triplets_c_term_x.clear();
    global_cij_y.setFromTriplets( triplets_c_term_y.begin(), triplets_c_term_y.end() );
    triplets_c_term_y.clear();
    */
    
    
    
    // INITIALISATION OF THE SOLVER (CONJUGATE GRADIENT)
    //ConjugateGradient<SparseMatrix<T> > solver_global_mass;
    
    timecounter tc1;
    tc1.tic();
    
    SimplicialLLT<SparseMatrix<T> >solver_global_mass;
    solver_global_mass.compute(global_mass); // use solver_global_mass to solve M^-1
    
    if(solver_global_mass.info()!=Success){
        std::cout<<"FAILED SOLVER 0->phi_tilde"<<std::endl;
        exit(1);
    }
    
    tc1.toc();
    std::cout << bold << yellow << "Inversion of the mass matrix: " << tc1 << " seconds" << reset << std::endl;
    
    /*
    SparseLU<SparseMatrix<T>, AMDOrdering<int> > solver2 ;
    solver2.compute(global_mass);
    if(solver2.info()!=Success) {
        std::cout<<"FAILED SOLVER2 PROVA ->phi_tilde"<<std::endl;
        return;
    }
    */
    
    /// CHECKING OF MASS MATRICES' PROPERTIES
    /*
    std::cout<<"First checking: Lumped Mass"<<std::endl;
    for(size_t i = 0; i<global_lumped_mass.rows();i++)
        if(global_lumped_mass(i)<0)
            std::cout<<"PROBLEM into lumped mass"<<std::endl;
    
    Matrix<T, Dynamic, Dynamic> mass_check = - global_mass;
    for(size_t i = 0 ; i < mass_check.rows() ; i++){
        mass_check(i,i) += global_lumped_mass(i);
        std::cout<<(mass_check.row(i)).sum()<<std::endl;
    }

    /// CHECK OF THE SUM PROPERTY OF CIJ
     Matrix<T, Dynamic, 1> sum0 = Eigen::Matrix<T, Dynamic, 1>::Zero(global_cij_x.rows(), 1);
     Matrix<T, Dynamic, 1> sum1 = Eigen::Matrix<T, Dynamic, 1>::Zero(global_cij_x.rows(), 1);
    
    for (size_t k = 0; k<global_cij_x.cols(); k++)
    {
        sum0 += (global_cij_x).col(k);
        sum1 += (global_cij_y).col(k);
    }
    
    for (size_t i=0; i<global_cij_y.rows(); i++) {
         std::cout<<"The sum is "<<sum1(i)<<" and "<<sum0(i)<<std::endl;
    }
    
    
    
    for(size_t counter_sum = 0 ; counter_sum < global_cij_x.rows() ; counter_sum++ )
    {
        for(size_t counter_col = 0 ; counter_col < global_cij_x.cols() ; counter_col++)
        {
    
           if(counter_col==counter_sum)
               std::cout<<"In ("<<counter_sum<<" , "<<counter_col<<" ), c^0 = "<<global_cij_x.coeff( counter_sum , counter_col )<<" and c1 = "<<global_cij_y.coeff( counter_sum , counter_col )<<std::endl;
           else
                std::cout<<"In ("<<counter_sum<<" , "<<counter_col<<" ), c^0_ij + c^0_ji = "<<global_cij_x.coeff( counter_sum , counter_col ) + global_cij_x.coeff( counter_col , counter_sum )<<" and c^1_ij + c^1_ji = "<<global_cij_y.coeff( counter_sum , counter_col ) + global_cij_y.coeff( counter_col , counter_sum ) <<std::endl;
        }
    }
    
 */
    
    // FLUX TERM
    /*
    Matrix<T, Dynamic, 1>  RHS_flux0 =  Matrix<T, Dynamic, 1>::Zero(dim) ;
    Matrix<T, Dynamic, 1>  RHS_flux1 =  Matrix<T, Dynamic, 1>::Zero(dim) ;
    size_t counter_flux = 0 ;
    for(auto& cl : msh.cells)
    {
        // FLUX TERM : flux is a pair flux0 and flux1
        auto local_flux = make_bernstein_local_RHS_FLUX( msh , cl , degree , u , phi ) ;
        //std::cout<<"flux_loc0  "<<'\n'<<flux_loc.first<<std::endl;
        //std::cout<<"flux_loc1  "<<'\n'<<flux_loc.second<<std::endl;
        for (size_t i = 0; i < local_ndof ; i++)
        {
            size_t asm_map =  phi.connectivity_matrix[counter_flux][i].first ;
            RHS_flux0(asm_map) += local_flux.first(i) ;
            RHS_flux1(asm_map) += local_flux.second(i) ;
        }
        counter_flux++;
    }
    Matrix<T, Dynamic, 1> flux0 = solver_global_mass.solve(RHS_flux0);
    Matrix<T, Dynamic, 1> flux1 = solver_global_mass.solve(RHS_flux1);
     */
    
    
    // ALTERNATIVE VANDERMONDE MATRIX
    
    //Matrix<T, Dynamic, 1> flux0 = u0.cwiseProduct( phi_FEM ) ;
    //Matrix<T, Dynamic, 1> flux1 = u1.cwiseProduct( phi_FEM );
    
    size_t i_fl = 0 ;
    Matrix<T, Dynamic, 1> flux0_loc = Matrix<T, Dynamic, 1>::Zero(local_ndof) ;
    Matrix<T, Dynamic, 1> flux1_loc = Matrix<T, Dynamic, 1>::Zero(local_ndof) ;
    
    Matrix<T, Dynamic, 1> flux0 = Matrix<T, Dynamic, 1>::Zero(dim) ;
    Matrix<T, Dynamic, 1> flux1 = Matrix<T, Dynamic, 1>::Zero(dim) ;
    
    //Matrix<T, Dynamic, Dynamic> local_vandermonde = Matrix<T, Dynamic, Dynamic>::Zero ( local_ndof , local_ndof) ;

    
    //T error0_L2 = 0.0 , error1_L2 = 0.0 ;
    
    CompleteOrthogonalDecomposition<Matrix<T, Dynamic, Dynamic > > cod( local_vandermonde );
    
    
    for(auto& cl : msh.cells)
    {
        // FLUX TERM : flux is a pair flux0 and flux1
        auto pts = equidistriduted_nodes_ordered_bis<T,Mesh>( msh, cl, degree );
        //cell_basis_Bernstein<Mesh,T> cb(msh, cl, degree);
        for (size_t i = 0; i < local_ndof ; i++)
        {
            flux0_loc(i) = u(pts[i] , msh , cl ).first * phi(pts[i] , msh , cl );
            flux1_loc(i) = u(pts[i] , msh , cl ).second * phi(pts[i] , msh , cl );
            // LOCAL VANDERMONDE MATRIX -> è uguale penso basti farla tra [0,1] e riportarla qua
            //local_vandermonde.block(i,0,1,local_ndof) = (cb.eval_basis(pts[i])).transpose() ;
            
            size_t asm_map =  phi.connectivity_matrix[i_fl][i].first ;
            
        }
        
        //std::cout<<'\n'<<'\n'<<"local_vandermonde "<<local_vandermonde<<'\n'<<std::endl;
        //CompleteOrthogonalDecomposition<Matrix<T, Dynamic, Dynamic > > cod( local_vandermonde );
        Matrix<T, Dynamic, 1> sol0 = cod.solve(flux0_loc);
        Matrix<T, Dynamic, 1> sol1 = cod.solve(flux1_loc);
        //std::cout<<"flux0_loc -sol0"<<'\n'<<flux0_loc-sol0<<'\n'<<std::endl;
        //if (cod.info() != Success)
        //{
        //    std::cout<<"Not positive"<<std::endl;
        //    assert(0);
        //}

        
        
        for (size_t i = 0; i < local_ndof ; i++)
        {
            
            size_t asm_map =  phi.connectivity_matrix[i_fl][i].first ;
            flux0(asm_map) = sol0(i) ;
            flux1(asm_map) = sol1(i) ;
            
        }
        
        
        i_fl++;
    }
    
    
    
    
    //std::cout<<"flux0_bis:"<<'\n'<<flux0_bis <<'\n' <<std::endl;
    //std::cout<<"flux0:"<<'\n'<< flux0<<'\n' <<std::endl;
/*
    T sum0 = 0.0 , sum1 = 0.0 ;
    for (int ii = 0 ; ii < dim ; ii++ )
    {
        sum0 += pow( (flux0_bis(ii) - flux0(ii)) , 2 );
        sum1 += pow( (flux1_bis(ii) - flux1(ii)) , 2 );
    }
    
    error0_L2 = sqrt( 1.0/dim * sum0 );
    error1_L2 = sqrt( 1.0/dim * sum1 );
    std::cout<<"ERRORE L2 (flux0_bis - flux0):"<<'\n'<<error0_L2<<'\n' <<std::endl;
    std::cout<<"ERRORE L2 (flux1_bis - flux1):"<<'\n'<<error1_L2<<'\n' <<std::endl;
*/
    
    
    // RESOLUTION OF phi_tilde (GLOBALLY)    ( with cij(phi_j) )
    Matrix<T, Dynamic, 1> mass_phi_old = global_mass * phi_FEM ; // USEFUL ALSO FOR HO PHI_H
    //std::cout<<"vec1  "<<'\n'<<mass_phi_old<<std::endl;
    
    // CONVOLUTION TERM
    Matrix<T, Dynamic, 1> conv_global = global_cij_x * flux0  + global_cij_y * flux1 ;
    //Matrix<T, Dynamic, 1> conv_global = global_cij_x * phi_FEM  + global_cij_y * phi_FEM ;
    /*
    Matrix<T, Dynamic, 1> flux0 = 10 * phi_FEM ;
    Matrix<T, Dynamic, 1> flux1 = 0 * phi_FEM ;
    Matrix<T, Dynamic, 1> conv_global2 = phi.Global_c_term_x * flux0  + phi.Global_c_term_y * flux1 ;
    std::cout<<"CHECK cij conv_global  "<<'\n'<<conv_global-conv_global2<<std::endl;
    */
    
    
    

       
    
    // NORM of c_ij
    SparseMatrix<T> cij_norm = ( global_cij_x.cwiseProduct(global_cij_x) + global_cij_y.cwiseProduct(global_cij_y) ).cwiseSqrt() ;
    //std::cout<<"cij norm "<<'\n'<<cij_norm<<std::endl;
    
    // MATRIX n_ij
    SparseMatrix<T> nij0 = global_cij_x.cwiseQuotient( cij_norm );
    SparseMatrix<T> nij1 = global_cij_y.cwiseQuotient( cij_norm );
    
    //std::cout<<"nij1  "<<'\n'<<nij1<<std::endl;
    

    // MATRIX c_ji
    SparseMatrix<T> cji_x = global_cij_x.adjoint() ;
    SparseMatrix<T> cji_y = global_cij_y.adjoint() ;
     
    // NORM of c_ji -> i.e. c_ij transposed
    SparseMatrix<T> cji_norm = (cji_x.cwiseProduct(cji_x)+cji_y.cwiseProduct(cji_y)).cwiseSqrt();
    
    // MATRIX n_ij (TRANSPOSED)
    SparseMatrix<T> nji0 = cji_x.cwiseQuotient( cji_norm );
    SparseMatrix<T> nji1 = cji_y.cwiseQuotient( cji_norm );
    
    
    // NORMAL VELOCITY
    SparseMatrix<T> normal_vel = SparseMatrix<T>( dim , dim );
    std::vector< Triplet<T> >   triplets;
    
    // NORMAL VELOCITY TRANSPOSED
    SparseMatrix<T> normal_vel_adj = SparseMatrix<T>( dim , dim );
    std::vector< Triplet<T> >   triplets_adj ;
    /*
    // NORMAL VELOCITY
    SparseMatrix<T> normal_vel_bis = SparseMatrix<T>( dim , dim );
    std::vector< Triplet<T> >   triplets0_bis;
    
    // NORMAL VELOCITY TRANSPOSED
    SparseMatrix<T> normal_vel_bis_adj = SparseMatrix<T>( dim , dim );
    std::vector< Triplet<T> >   triplets1_bis ;
    */
    size_t counter_row = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            auto value0 = std::abs( u0(counter_row) * nij0.coeff(counter_row,elem) + u1(counter_row) * nij1.coeff(counter_row,elem) );
            auto value1 = std::abs( u0(elem) * nij0.coeff(counter_row,elem) + u1(elem) * nij1.coeff(counter_row,elem) );
            auto value = std::max(value0 , value1);
            
            auto value_adj0 = std::abs( u0(counter_row) * nji0.coeff(counter_row,elem) + u1(counter_row) * nji1.coeff(counter_row,elem) );
            auto value_adj1 = std::abs( u0(elem) * nji0.coeff(counter_row,elem) + u1(elem) * nji1.coeff(counter_row,elem) );
            auto value_adj = std::max(value_adj0 , value_adj1);
               
            triplets.push_back( Triplet<T>(counter_row, elem, value ) );
            triplets_adj.push_back( Triplet<T>(counter_row, elem, value_adj ) );
            
            
            // CHECK IT
            /*
            auto value0_bis = std::max( u0(counter_row) , u0(elem) );
            auto value1_bis = std::max( u1(counter_row) , u1(elem) );
            auto val = value0_bis * nij0.coeff(counter_row,elem) + value1_bis * nij1.coeff(counter_row,elem) ;
            auto val_adj = value0_bis * nji0.coeff(counter_row,elem) + value1_bis * nji1.coeff(counter_row,elem) ;
            
            triplets0_bis.push_back( Triplet<T>(counter_row, elem, val ) );
            triplets1_bis.push_back( Triplet<T>(counter_row, elem, val_adj ) );
            */
            
               
        }
        counter_row++;
    }
    
    normal_vel.setFromTriplets( triplets.begin(), triplets.end() );
    triplets.clear();
    normal_vel_adj.setFromTriplets( triplets_adj.begin(), triplets_adj.end() );
    triplets_adj.clear();
    
    /*
    normal_vel_bis.setFromTriplets( triplets0_bis.begin(), triplets0_bis.end() );
    triplets0_bis.clear();
    normal_vel_bis_adj.setFromTriplets( triplets1_bis.begin(), triplets1_bis.end() );
    triplets1_bis.clear();
    SparseMatrix<T> normal_vel_abs_bis = normal_vel_bis.cwiseAbs() ;  // o lambda_max=normal_vel ?
    SparseMatrix<T> normal_vel_abs_adj_bis = normal_vel_bis_adj.cwiseAbs();
    */
    // CHECKING NEW METHOD NORMAL VEL
    /*
    for(size_t i = 0 ; i < normal_vel_abs_bis.rows() ; i++)
    {
        for(size_t j = 0 ; j < normal_vel_abs_bis.cols() ; j++){
            if(std::abs(normal_vel.coeff(i,j)-normal_vel_abs_bis.coeff(i,j))>1e-5)
                std::cout<<"In ("<<i<<" , "<<j<<" ), normal vel OLD is "<<normal_vel_abs_bis.coeff(i,j)<<" , normal vel NEW is "<<normal_vel.coeff(i,j)<<std::endl;
            
            if(std::abs(normal_vel_adj.coeff(i,j)-normal_vel_abs_adj_bis.coeff(i,j))>1e-5)
                std::cout<<"In ("<<i<<" , "<<j<<" ), normal vel ADJ OLD is "<<normal_vel_abs_adj_bis.coeff(i,j)<<" , normal vel  ADJ NEW is "<<normal_vel_adj.coeff(i,j)<<std::endl;
            
        }
    }
    */
    //SparseMatrix<T> normal_vel2 = u0 * nij0 + u1 * nij1 ;
    //SparseMatrix<T> normal_vel_adj2 = u0 * nji0 + u1 * nji1 ;
    
    //std::cout<<"normal_vel - nij0"<<'\n'<<normal_vel - nij0<<std::endl;
    //std::cout<<"normal_vel_adj - nij0"<<'\n'<<normal_vel_adj - nji0<<std::endl;
    //std::cout<<"CHECK normal_vel2  : "<<'\n'<<normal_vel-normal_vel2<<'\n'<<std::endl;
    //std::cout<<"CHECK normal_vel_adj2 : "<<'\n'<<normal_vel_adj-normal_vel_adj2<<'\n'<<std::endl;
   
    // MATRIX Dij CALCULATION:
    // normval_vel and the adjoint one are already absolute value
    
    //SparseMatrix<T> normal_vel_abs = normal_vel.cwiseAbs() ;  // o lambda_max=normal_vel ?
    //SparseMatrix<T> normal_vel_abs_adj = normal_vel_adj.cwiseAbs();// o lambda_max=normal_vel?
 
    
    /// POSSO SCEGLIERE SE METTERE normal_vel_abs_bis O normal_vel !!!!!!!!!!!!!!!!
    //std::cout<<bold<<yellow<<"--->  sto usando normal_vel_abs_bis, i.e. con max su u(pt)!!"<<reset<<std::endl;
     std::cout<<bold<<yellow<<"--->  sto usando normal_vel, i.e. con max su u*n !!" <<reset<<std::endl;
    
    
    //SparseMatrix<T> lambda_max = normal_vel_abs_bis.cwiseProduct( cij_norm );
    //SparseMatrix<T> lambda_max_adj = normal_vel_abs_adj_bis.cwiseProduct( cji_norm );
    SparseMatrix<T> lambda_max = normal_vel.cwiseProduct( cij_norm );
    SparseMatrix<T> lambda_max_adj = normal_vel_adj.cwiseProduct( cji_norm );
    
    
    SparseMatrix<T> dij = lambda_max.cwiseMax(lambda_max_adj);
    //SparseMatrix<T> dij2 = lambda_max.cwiseMax(lambda_max_adj);
    
    for(size_t i = 0 ; i < dij.rows() ; i++)
    {
        dij.coeffRef(i,i) = 0;
        dij.coeffRef(i,i) = -dij.row(i).sum();
    }
    //std::cout<<"dij-dij2"<<'\n'<<dij-dij2<<std::endl;
    
    // CHECK TIME STEP dt
    T dt_old = dt ;
    std::cout<<bold<<yellow<<"---> COND IN TEMPO CFL, ALEXANDRE BOOK"<<reset<<std::endl;
    T CFL_numb = time_step_CFL_L2_velocity_NEW( dij.diagonal() , global_lumped_mass , fe_data.Dirichlet_boundary , dt );
    
    //T nu_max0 = CFL_numb/fe_data.hx;
    //T nu0 = dt_old/fe_data.hx;
    //T nu1 = dt/fe_data.hx;
    
    //std::cout<<"VALID FOR u = (1,0). nu_max VERO = "<<nu_max0<<" , nu max con dt assegnato = "<<nu0<< " and with dt appeared by CFL COND "<<nu1<<std::endl;
    //T nu_max = dt/fe_data.hx;
    //std::cout<<"VALID FOR u = (1,0). nu_max = "<<nu_max<<std::endl;
    if(dt_old != dt )
    {
        std::cout<<"dt is "<<dt_old<<" and dt CFL is "<<dt<<" . STOP!"<<std::endl;
        exit(10);
    }
    
    
    
    
    // CONSTANT TERM (PHI TILDE PROBLEM)
    //Matrix<T,Dynamic,1> vec2 = Matrix<T, Dynamic, 1>::Zero(dim) ;
    //averaged_sum_Si_SPARSE(S_i , global_cij_x , flux0 , vec2);
    //averaged_sum_Si_SPARSE(S_i , global_cij_y , flux1 , vec2);
    
    //std::cout<<"CHECK phi TILDE "<<'\n'<<conv_global-vec2<<'\n'<<"FINE PHI TILDE CHECKING"<<std::endl;
    Matrix<T, Dynamic, 1> b = mass_phi_old - dt*conv_global.cwiseQuotient(global_lumped_mass);
    //std::cout<<"TERMINE NOTO:b  "<<'\n'<<b<<std::endl;
    
    // RESOLUTION OF PHI_TILDE
    phi_tilde.sol_FEM = solver_global_mass.solve(b); // SAVE THE L2 projection
    
    //auto prova_phi_tilde = solver2.solve(b);
    //std::cout<<"phi_tilde.sol_FEM  "<<'\n'<<phi_tilde.sol_FEM <<std::endl;
    
    //std::cout<<"phi_tilde.sol_FEM - prova_phi_tilde "<<'\n'<<phi_tilde.sol_FEM - prova_phi_tilde <<std::endl;
    
    
    if(solver_global_mass.info()!=Success) {
        std::cout<<"FAILED SOLVER 1->phi_tilde"<<std::endl;
        exit(1);
    }
    
    // SAVING BOTH SOL_HHO AND VERTICES OF PHI_TILDE
    std::cout<<"CONVERTING phi_tilde"<<std::endl;
    phi_tilde.converting_into_HHO_formulation( phi_tilde.sol_FEM );
    
    
    
    // TERM R_i^n
    Matrix<T, Dynamic, 1> Emax_global = Eigen::Matrix<T, Dynamic, 1>::Ones(dim, 1);
    Matrix<T, Dynamic, 1> Emin_global = Eigen::Matrix<T, Dynamic, 1>::Ones(dim, 1);
    Emax_global *= -1e20;
    Emin_global *= 1e20;
    Matrix<T, Dynamic, 1> R_i = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    
    for( auto& cl: msh.cells )
        r_i_calculator_Bernstein( msh , cl , E , phi_tilde , phi , dt , u , Emax_global , Emin_global , R_i);
    
    // R_i FINALISATION:
    //std::cout<<"Emax_global: "<<'\n'<<Emax_global<<std::endl;
    //std::cout<<"Emin_global: "<<'\n'<<Emin_global<<std::endl;
    R_i = R_i.cwiseQuotient( Emax_global - Emin_global );
    //std::cout<<"Emax_global - Emin_global = "<<'\n'<<Emax_global-Emin_global<<std::endl;
    //std::cout<<"Ri = "<<'\n'<<R_i<<std::endl;
    
    
    
    // ENTROPIC SOLUTION: MINIMUM BETWEEN d_ij AND R_i --> d^E_ij MATRIX
    T c_e = 1.0;
    SparseMatrix<T> dE_ij  = SparseMatrix<T>( dim , dim );
    std::vector< Triplet<T> >   triplets_dE_ij;
    SparseMatrix<T> phi_ij = SparseMatrix<T>( dim , dim );
    std::vector< Triplet<T> >   triplets_phi_ij;
    
    
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            
            if(elem!=counter)
            {
                auto R_i_j = c_e * std::max(std::abs( R_i(counter) ) , std::abs( R_i(elem) ));
                auto value = std::min( dij.coeff( counter , elem ), R_i_j );
                triplets_dE_ij.push_back( Triplet<T>(counter, elem, value ) );
            //dE_ij( counter_row , elem ) = std::min( dij( counter_row , elem ), R_i_j );
            //std::cout<<"R_i_j "<<R_i_j<<" and dij"<<dij( counter_row , elem )<<std::endl;
            }
            auto value_bis = 0.5*( phi_FEM(counter) + phi_FEM(elem) );
            triplets_phi_ij.push_back( Triplet<T>(counter, elem, value_bis ) );
            
            
        }
        counter++;
    }
    
    //std::cout<<"term_dij_no_entropy : "<<'\n'<<term_dij_no_entropy<<std::endl;
    
    dE_ij.setFromTriplets( triplets_dE_ij.begin(), triplets_dE_ij.end() );
    triplets_dE_ij.clear();
    phi_ij.setFromTriplets( triplets_phi_ij.begin(), triplets_phi_ij.end() );
    triplets_phi_ij.clear();
    
    //std::cout<<"dE_ij: "<<'\n'<<dE_ij<<std::endl;
    
    // MATRIX d^C_ij
    
    T c_comp = 1.0;
    SparseMatrix<T> tmp0 = phi_ij - phi_ij.cwiseProduct(phi_ij);
    //std::cout<<"tmp0 "<<'\n'<<tmp0<<std::endl;
    positive_part_SPARSE( tmp0 );
    //std::cout<<"tmp0 POSITIVE "<<'\n'<<tmp0<<std::endl;
    
    SparseMatrix<T> tmp_dC1 =  SparseMatrix<T>( dim , dim );
    std::vector< Triplet<T> >   triplets_tmp_dC1;
    SparseMatrix<T> Mat_One =  SparseMatrix<T>( dim , dim );
    //Matrix<T, Dynamic, Dynamic> Mat_One = Eigen::Matrix<T, Dynamic, Dynamic>::Ones( dim, dim );
    std::vector< Triplet<T> >   triplets_1;
    
    size_t counter_tmp = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
           // std::cout<<"In ("<<counter<<" , "<<elem<<") tmp_dC0 = "<<tmp_dC0(counter,elem)<<" and phi.ver diff "<<std::abs(phi.vertices(counter) - phi.vertices(elem))<<std::endl;
            if( std::abs(phi_FEM(counter_tmp) - phi_FEM(elem))>1e-15 ){
                auto value = tmp0.coeff(counter_tmp,elem)/( std::abs(phi_FEM(counter_tmp) - phi_FEM(elem)) );
                triplets_tmp_dC1.push_back( Triplet<T>(counter_tmp, elem, value ) );
                triplets_1.push_back( Triplet<T>(counter_tmp, elem, 1.0 ) );
            }
            //tmp_dC1(counter_tmp,elem) = tmp0(counter_tmp,elem)/( std::abs(sol_FEM(counter_tmp) - sol_FEM(elem)) );
        }
        counter_tmp++;
    }
    tmp_dC1.setFromTriplets( triplets_tmp_dC1.begin(), triplets_tmp_dC1.end() );
    triplets_tmp_dC1.clear();
    
    
    Mat_One.setFromTriplets( triplets_1.begin(), triplets_1.end() );
    
    
    Matrix<T , Dynamic , Dynamic > Mat1 = Eigen::Matrix<T, Dynamic, Dynamic>::Ones(dim, dim);
    Matrix<T , Dynamic , Dynamic > tmp1_bis = Mat1 - c_comp*tmp_dC1;
    positive_part( tmp1_bis );
    Matrix<T , Dynamic , Dynamic > dC_ij_bis = dE_ij.cwiseProduct(tmp1_bis) ;


    
    //std::cout<<'\n'<<"tmp_dC1: "<<'\n'<<tmp_dC1<<std::endl;
    SparseMatrix<T> tmp1 = Mat_One - c_comp*tmp_dC1;
    
    //Matrix<T, Dynamic, Dynamic> tmp1 = Mat_One - c_comp*tmp_dC1;
    //positive_part( tmp1  );
    //SparseMatrix<T> dC_ij = SparseMatrix<T>( dim , dim );
    //size_t counter_prova = 0;
    //for(auto& row_i:S_i)
    //{
     //   for(auto& elem:row_i)
     //   {
     //       dC_ij.coeffRef(counter_prova,elem) += ( dE_ij.coeff(counter_prova,elem)*tmp1(counter_prova,elem) );
     //   }
      //  counter_prova++;
    //}
    
   
    
    //std::cout<<'\n'<<"tmp1: "<<'\n'<<tmp1<<std::endl;
    positive_part_SPARSE( tmp1  );
    
    SparseMatrix<T> dC_ij = dE_ij.cwiseProduct(tmp1) ;
    
    /*
    SparseMatrix<T> dC_ij_prova = dC_ij ;
    // NON SONO SICURO SERVA STA COSA
    std::cout<<"CHECK IF SERVE L'IMPOSIZIONE PER L'ELEMENTO DIAGONALE = SOMMA RIGA PER VARI d^I_ij!!!!"<<std::endl;
    for(size_t i = 0 ; i < dC_ij.rows() ; i++)
    {
        dC_ij.coeffRef(i,i) = 0.;
        dij.coeffRef(i,i) = 0.;
        dC_ij.coeffRef(i,i) = -dC_ij.row(i).sum();
        dij.coeffRef(i,i) = -dij.row(i).sum();
    }
    
    std::cout<<'\n'<<"dC_ij_prova - dC_ij: "<<'\n'<<dC_ij_prova - dC_ij<<std::endl;
    */
    
    // ARTIFICIAL VISCOSITY TERM d_ij * (phi_j-phi_i)
    
    Matrix<T, Dynamic, 1> term_dij_no_entropy =  Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    Matrix<T, Dynamic, 1> term_dij_E            =  Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    Matrix<T, Dynamic, 1> term_dij            =  Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    
    //Matrix<T, Dynamic, 1> term_dij_bis =  Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    
    size_t counter_dij = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
             // TERM d_ij * (phi_j-phi_i)
            term_dij_no_entropy(counter_dij) += dij.coeff(counter_dij,elem)*(phi_FEM(elem)-phi_FEM(counter_dij));
             // TERM dC_ij * (phi_j-phi_i)
            term_dij(counter_dij) += dC_ij.coeff(counter_dij,elem)*(phi_FEM(elem)-phi_FEM(counter_dij));
            term_dij_E(counter_dij) += dE_ij.coeff(counter_dij,elem)*(phi_FEM(elem)-phi_FEM(counter_dij));
            
            //term_dij_bis(counter_dij) += dC_ij_bis.coeff(counter_dij,elem)*(phi_FEM(elem)-phi_FEM(counter_dij));
        }
        counter_dij++;
    }
    
    //std::cout<<"term_dij_bis - term_dij"<<'\n'<<term_dij_bis - term_dij<<std::endl;
    
    // FOR CHECKING
    //std::cout<<"CHECK term_dij_no_entropy : "<<'\n'<<term_dij_no_entropy - term_dij<<std::endl;
    
    //std::cout<<"CHECK term_dij_no_entropy : "<<'\n'<<term_dij_no_entropy - term_dij_no_entropy2<<std::endl;
    
    //Matrix<T, Dynamic, 1> term_dij_no_entropy2 = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1 );
    //term_dij_no_entropy2 = dij * phi_FEM ;
    
    //Matrix<T, Dynamic, 1> term_dij2 = dC_ij * phi_FEM ;
    //std::cout<<"CHECK term_dij2 : "<<'\n'<<term_dij2 - term_dij<<std::endl;
    
    //std::cout<<"CHECK term_dij_no_entropy : "<<'\n'<<term_dij_no_entropy - term_dij_no_entropy2<<std::endl;
    
    
    
    /// MATRIX D_IJ CHECKING
    /*
     //The  error is < 1e-15, at max
    std::cout<<"D_ij check: "<<std::endl;
    for(size_t i = 0; i<lambda_max.rows();i++){
        for(size_t j = i+1; j<lambda_max.rows();j++){
            if(std::abs(dij.coeff(i,j)-dij.coeff(j,i))>1e-15)
            std::cout<<dij.coeff(i,j)-dij.coeff(j,i);
        }
        std::cout<<'\n';
    }
    std::cout<<std::endl;
    
    
     std::cout<<'\n'<<"dij: "<<'\n'<<dij<<std::endl;
    
    std::cout<<"D_ij check symmetry "<<std::endl;
    for(size_t counter_row = 0 ; counter_row < global_cij_x.rows() ; counter_row++ )
    {
        
        for(size_t elem = 0 ; elem < global_cij_x.cols() ; elem++)
        {
            std::cout<<dij.coeff(counter_row,elem) - dij.coeff(elem,counter_row)<<std::endl;
        }
        std::cout<<'\n';
    }
    std::cout<<std::endl;
    */
    
    /// TERM D_IJ*(PHI_J-PHI_I)  CHECKING
    /*
    Matrix<T, Dynamic, 1> term_dij_no_entropy = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    for(size_t i = 0 ; i < dij.rows() ; i++)
    {
        for(size_t j = 0 ; j < dij.cols() ; j++)
        {
            term_dij_no_entropy.coeffRef(i) += dij.coeff(i,j) * ( phi_FEM(j)-phi_FEM(i) ) ;
        }
    }
    */
    
    ///********* RESOLUTION OF THE SYSTEM: **********//
   
    
    // RESOLUTION FIRST ORDER
    Matrix<T, Dynamic, 1> phi_L = phi_FEM - dt * conv_global.cwiseQuotient(global_lumped_mass)  + dt * term_dij_no_entropy.cwiseQuotient(global_lumped_mass);
  
    //std::cout<<'\n'<<"phi_L : "<<'\n'<<phi_L<<std::endl;
    
    
    // RESOLUTION HIGH ORDER -> NO MIN MAX PRINCIPLE PRESERVING
    Matrix<T, Dynamic, 1> b_phiH = mass_phi_old - dt * conv_global + dt * term_dij ;
    Matrix<T, Dynamic, 1> phi_H = solver_global_mass.solve(b_phiH);
    
    Matrix<T, Dynamic, 1> b_phiE = mass_phi_old - dt * conv_global + dt * term_dij_E ;
    Matrix<T, Dynamic, 1> phi_E = solver_global_mass.solve(b_phiE);
    
    
    // CHECKING OF phi_h AND phi_l
    //Matrix<T, Dynamic, 1> check0 = global_lumped_mass.cwiseProduct( (phi_L - phi_FEM)/dt ) + conv_global - term_dij_no_entropy;
    //std::cout<<'\n'<<"CHECK phi_L is "<<'\n'<<check0<<std::endl;
    
    //Matrix<T, Dynamic, 1> check1 = global_mass*( (phi_H - phi_FEM)/dt ) + conv_global - term_dij;
    //std::cout<<'\n'<<"CHECK phi_H is "<<'\n'<<check1<<std::endl;
    
    
    // EXTENSION: MAXIMUM PRINCIPLE PRESERVING
    
    Matrix<T, Dynamic, 1> delta_phi = phi_H - phi_FEM;
    //std::cout<<'\n'<<"delta_phi "<<'\n'<<delta_phi<<std::endl;
    
    Matrix<T, Dynamic, 1> f_i = f_ij_creator_SPARSE( global_lumped_mass , global_mass , delta_phi , dt , dij , dC_ij , phi_FEM , S_i );
    
    // CHECKING PHI
    /*
    auto check_phi = phi_H - phi_L - f_i.cwiseQuotient(global_lumped_mass) ;//+dt*Vec_One ;
    
    Matrix<T, Dynamic, 1>  phi_H2 = phi_L + f_i.cwiseQuotient(global_lumped_mass);
    
    std::cout<<'\n'<<"check phi "<<'\n'<<check_phi<<std::endl;
    std::cout<<'\n'<<"check phi_h: "<<'\n'<<phi_H - phi_H2<<std::endl;
    */
    
    
    // CORRECTION TERM
    Matrix<T, Dynamic, 1> correction_fi = alfaf_ij_creator_SPARSE( global_lumped_mass , global_mass , delta_phi , phi_L , dt , dij , dC_ij , phi_FEM , S_i );
    Matrix<T, Dynamic, 1>  phi_new = phi_L + correction_fi.cwiseQuotient(global_lumped_mass);
    
    
    
    
    // IMPOSITION DIRICHLET BOUNDARY CONDITIONS
    size_t counter_dir = 0 ;
    for (const auto& dir_elem : fe_data.Dirichlet_boundary )
    {
        if(dir_elem){
            phi_L(counter_dir) = phi_FEM(counter_dir) ;
            phi_H(counter_dir) = phi_FEM(counter_dir) ;
            phi_new(counter_dir) = phi_FEM(counter_dir) ;
            phi_E(counter_dir) = phi_FEM(counter_dir) ;
        }
        counter_dir++ ;
    }
     
    
    //std::cout<<"phi_L - phi_H"<<'\n'<<phi_L - phi_H<<'\n'<<'\n'<<'\n'<<std::endl;
    //std::cout<<"phi_L - phi_new"<<'\n'<<phi_L - phi_new<<'\n'<<'\n'<<'\n'<<std::endl;
    //std::cout<<"phi_H - phi_new"<<'\n'<<phi_H - phi_new<<'\n'<<'\n'<<'\n'<<std::endl;
    
    
    
    tc.toc();
    std::cout << bold << yellow << "FEM method, time resolution: " << tc << " seconds" << reset << std::endl;
    
    
    
    
    
    /// PLOTTING SOLUTION (GNUPLOT) + SAVING FOR HHO (MISCHIATO PER POTERE PLOTTARE ENTRAMBE).
    postprocess_output<double> postoutput5;
    auto test_phi_h = std::make_shared< gnuplot_output_object<double> >("phi_h.dat");
    auto test_phi_l = std::make_shared< gnuplot_output_object<double> >("phi_l.dat");
    auto test_phi_new = std::make_shared< gnuplot_output_object<double> >("phi_new.dat");
    
    // SAVING AND UPLOAD phi_L  INTO CLASS projected_level_set
    
   /*
   

    phi.sol_FEM = phi_L ;
    
    std::cout<<"USING LOW ORDER SOLUTION PHI_L. Also vertices uploading."<<std::endl;
    //phi.converting_into_HHO_formulation(phi_L);
    phi.converting_into_HHO_formulation(phi_L);
    //std::cout<<"phi max "<<phi.sol_FEM.maxCoeff()<<'\n';
    //std::cout<<"phi min "<<phi.sol_FEM.minCoeff()<<'\n';
    //std::cout<<"SAVED phi max "<<phi.phi_max<<'\n';
    //std::cout<<"SAVED phi min "<<phi.phi_min<<std::endl;
    
    phi.coefficients_inverse_mapping();
    //std::cout<<"phi max "<<phi.sol_FEM.maxCoeff()<<'\n';
    //std::cout<<"phi min "<<phi.sol_FEM.minCoeff()<<'\n';
    //std::cout<<"SAVED phi max "<<phi.phi_max<<'\n';
    //std::cout<<"SAVED phi min "<<phi.phi_min<<std::endl;
    
    //phi.set_max_min(); // SETTING OF THE NEW MAX AND MIN FOR NEXT TRANSPORT PROBLEM
    //std::cout<<"phi max "<<phi.sol_FEM.maxCoeff()<<'\n';
    //std::cout<<"phi min "<<phi.sol_FEM.minCoeff()<<'\n';
    //std::cout<<"SAVED phi max "<<phi.phi_max<<'\n';
    //std::cout<<"SAVED phi min "<<phi.phi_min<<std::endl;
    
    for(auto& cl :msh.cells)
    {
        auto pts = equidistriduted_nodes_ordered_bis<T,Mesh>( msh, cl, degree );
        for (auto pt : pts){
            test_phi_l->add_data( pt ,phi(pt, msh , cl ) );
        }
    }
    
     */
    
    /*
   
 
     
    phi.sol_FEM = phi_E ;
       
    std::cout<<"CONVERTING phi_E in HHO formulation. Also vertices uploading."<<std::endl;
    phi.converting_into_HHO_formulation(phi_E);
       
       
    phi.coefficients_inverse_mapping();
    //phi.set_max_min(); // SETTING OF THE NEW MAX AND MIN FOR NEXT TRANSPORT PROBLEM
       
    
   
    
  
    // SAVING AND UPLOAD phi_L  INTO CLASS projected_level_set
    phi.sol_FEM = phi_H ;
       
    
    std::cout<<"CONVERTING phi_L in HHO formulation. Also vertices uploading."<<std::endl;
    //phi.converting_into_HHO_formulation(phi_L);
    phi.converting_into_HHO_formulation(phi_H);
       
       
    phi.coefficients_inverse_mapping();
    //phi.set_max_min(); // SETTING OF THE NEW MAX AND MIN FOR NEXT TRANSPORT PROBLEM

    for(auto& cl :msh.cells)
    {
        auto pts = equidistriduted_nodes_ordered_bis<T,Mesh>( msh, cl, degree );
        for (auto pt : pts){
            test_phi_h->add_data( pt ,phi(pt, msh , cl ) );
        }
    }
    
    
     */
 
   
    // SAVING AND UPLOAD phi_L  INTO CLASS projected_level_set
    phi.sol_FEM = phi_new ;
    
    std::cout<<"USING HIGH ORDER SOLUTION PHI_HP. Also vertices uploading."<<std::endl;
    //phi.converting_into_HHO_formulation(phi_L);
    phi.converting_into_HHO_formulation(phi_new);
    
    
    phi.coefficients_inverse_mapping();
    //phi.set_max_min(); // SETTING OF THE NEW MAX AND MIN FOR NEXT TRANSPORT PROBLEM
    
    
    
    
    for(auto& cl :msh.cells)
    {
        auto pts = equidistriduted_nodes_ordered_bis<T,Mesh>( msh, cl, degree );
        for (auto pt : pts){
            //std::cout<<pt<<std::endl;
            //test_phi_h->add_data( pt , phi.sol_HHO(iii , counter_cl ) );
            test_phi_new->add_data( pt , phi(pt, msh , cl ) );
        }
    }
    
  
    
    /// PLOTTING SOLUTION (GNUPLOT)
    /*
    postprocess_output<double> postoutput5;
    auto test_phi_h = std::make_shared< gnuplot_output_object<double> >("phi_h.dat");
    auto test_phi_l = std::make_shared< gnuplot_output_object<double> >("phi_l.dat");
    //auto test_phi_new = std::make_shared< gnuplot_output_object<double> >("phi_new.dat");
   
    size_t counter_cl = 0;
    for(auto& cl :msh.cells)
    {
        auto pts = equidistriduted_nodes_ordered_bis<T,Mesh>( msh, cl, degree );
        size_t iii = 0;
        for (auto pt : pts){
            //std::cout<<pt<<std::endl;
            test_phi_h->add_data( pt , phi.sol_HHO(iii , counter_cl ) );
            test_phi_l->add_data( pt , phi.sol_HHO(iii , counter_cl ) );
            iii++;
        }
        counter_cl++;
    }
     */
    postoutput5.add_object(test_phi_h);
    postoutput5.add_object(test_phi_l);
    postoutput5.add_object(test_phi_new);
    postoutput5.write();
    
    std::cout<<yellow<<bold<<"----------- FINE TRANSPORT PROBLEM -----------"<<reset<<std::endl;

    //return phi_tilde;
    
}







template < typename Fonction, typename Mesh, typename Vel_Field , typename T = typename Mesh::coordinate_type >
void run_FEM_levelset(const Mesh & msh, size_t degree, Fonction & phi, const Vel_Field & u , const T& dt , const mesh_init_params<T>& params)
{
    timecounter tc;
    tc.tic();
    size_t ndof = ( degree+1 ) * ( degree+1 ) ; // degrees of freedom
    const T eps = 1e-14; //constant into entropy
 
    // LOCAL MASS MATRIX ==> CHECK DEGREE == 1!!!
    assert(degree==1);
    
    
    mapping_phi(phi , phi.phi_max , phi.phi_min ); // mapping of phi to  have a phi between 0 and 1
    /// PHI MAPPED INTO 0-1 PLOT
    /*
    postprocess_output<double> postoutput4;
    auto test_phi_mapped = std::make_shared< gnuplot_output_object<double> >("phi_mapped.dat");
    for (auto cl:msh.cells) {
        auto pts = points(msh,cl) ;
        for (auto pt:pts) {
            T value = phi(pt,msh,cl);
            test_phi_mapped->add_data(pt,value);
        }
    }
    postoutput4.add_object(test_phi_mapped);
    postoutput4.write();
    */
    
    // ENTROPY -> for phi with codomain in (0,1)
    non_linear_entropy<T,Fonction,Mesh> E(eps , phi ,msh );
    auto phi_tilde = projection<T,Mesh> ( msh, degree , params );
    
    /// ENTROPY PLOT
    /*
    postprocess_output<double> postoutput6;
    auto cut_entropy = std::make_shared< gnuplot_output_object<double> >("entropy.dat");
    for (auto cl:msh.cells) {
        auto pts = points(msh,cl) ;
        for (auto pt:pts) {
            T value = E(pt,cl);
            cut_entropy->add_data(pt,value);
        }
    }
    postoutput6.add_object(cut_entropy);
    postoutput6.write();
    */

    // Global Matrix and Global Vector Definitions:
    size_t dim = msh.nodes.size();
    // Lumped mass vector
    Matrix<T, Dynamic, 1> global_lumped_mass = Eigen::Matrix<T, Dynamic, 1>::Zero( dim, 1 );
    // Mass matrix
    Matrix<T, Dynamic, Dynamic> global_mass = Eigen::Matrix<T, Dynamic, Dynamic>::Zero( dim, dim );
    // C_ij,1 matrix
    Matrix<T, Dynamic, Dynamic> global_cij_x = Eigen::Matrix<T, Dynamic, Dynamic>::Zero( dim, dim );
    // C_ij,2 matrix
    Matrix<T, Dynamic, Dynamic> global_cij_y = Eigen::Matrix<T, Dynamic, Dynamic>::Zero( dim, dim );
    // Convection term (vector)
    Matrix<T, Dynamic, 1> conv_global = Eigen::Matrix<T, Dynamic, 1>::Zero( dim, 1 );
    // D_ij term (vector)
    Matrix<T, Dynamic, 1> term_dij = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1 );
    // N_ij matrix
    std::pair<Matrix<T, Dynamic, Dynamic>,Matrix<T, Dynamic, Dynamic>> nij = std::make_pair(Eigen::Matrix<T, Dynamic, Dynamic>::Zero( dim, dim ) , Eigen::Matrix<T, Dynamic, Dynamic>::Zero( dim, dim ) ) ;
    // N_ji matrix
    std::pair<Matrix<T, Dynamic, Dynamic>,Matrix<T, Dynamic, Dynamic>> nij_transpose = std::make_pair(Eigen::Matrix<T, Dynamic, Dynamic>::Zero( dim, dim ) , Eigen::Matrix<T, Dynamic, Dynamic>::Zero( dim, dim ) ) ;
    // ||c_ij||^2
    Matrix<T, Dynamic, Dynamic> cij_square = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim, dim);
    // ||c_ji||^2
    Matrix<T, Dynamic, Dynamic> cji_square = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim, dim);
    
    // One vector
    Matrix<T, Dynamic, 1> Vec_One = Eigen::Matrix<T, Dynamic, 1>::Ones( dim, 1 );
    // One matrix
    Matrix<T, Dynamic, Dynamic> Mat_One = Eigen::Matrix<T, Dynamic, Dynamic>::Ones( dim, dim );
    
    
    // Boundary Nodes = Vertices (in Q1)
    auto bdry_nodes = boundary_nodes_function( phi.Nx , phi.Ny );

    
    // The Support of Global Shape Function: list of neighborood nodes for each node i
    std::vector<std::set<size_t>> S_i(dim);
    timecounter tc1;
    tc1.tic();
    for (auto cl : msh.cells)
    {
        
        
        auto nodes_position = nodes(msh, cl); // node_pos.ptid to have just the number
        
        // In basis (mesh creation) the order is
        /// 2 - - 3
        /// |       |
        /// 0 - - 1
        
        // In nodes position instead the order is
        /// 3 - - 2
        /// |       |
        /// 0 - - 1
        // I swap last two elements to have the first ordering
        iter_swap(nodes_position.end()-1 , nodes_position.end()-2 );
       
        // Support of global shape function: list of close nodes for each node i
        supporting_nodes( S_i , nodes_position );
        
        /*
        // CALCULATION OF THE SIZE + PLOTTING
        size_t size_supp_nodes = 0;
        std::cout<<"Supporting nodes:"<<std::endl;
        size_t jjjj = 0;
        for (auto& i:S_i) {
            size_supp_nodes+=i.size();
            std::cout <<"Node "<<jjjj<<":";
            for (auto it=i.begin(); it != i.end(); ++it)
                std::cout << ' ' << *it;
               // std::cout<<ii;
            std::cout<<'\n';
            jjjj++;
        }
        std::cout<<std::endl;
        std::cout<<"Supporting nodes size:"<<size_supp_nodes<<std::endl;
        */
        
        
        // Local mass matrix
        Matrix<T, Dynamic, Dynamic> mass = make_lagrange_local_mass_matrix(msh, cl, degree);
        
        // Assembling into global matrix
        //global_update2( mass , global_mass , nodes_position2 );
        global_update( mass , global_mass , nodes_position );
        
        // Local c_ij = b_i nabla(b_j)
        auto cij = make_local_cij_matrix(msh, cl, degree);
        
        // Assembling into global matrix
        global_update( cij.first , global_cij_x , nodes_position );
        global_update( cij.second , global_cij_y , nodes_position );
        
    // USEFUL COMMAND FOR NORMS of a Eigen::matrix A:
    // A.template lpNorm<2>()
    // A.norm()
        
    }
    
    // LUMPED MASS VECTOR (i.e. vector = diagonal of a matrix)
    sum_Si(S_i, global_mass , Vec_One , global_lumped_mass );
    
    
    
    // CHECKING OF MASS MATRICES' PROPERTIES-> OK, perfect precision 1e-19
    /*
    std::cout<<"First checking: Lumped Mass"<<std::endl;
    for(size_t i = 0; i<global_lumped_mass.rows();i++)
        if(global_lumped_mass(i)<0)
            std::cout<<"PROBLEM into lumped mass"<<std::endl;
    
    Matrix<T, Dynamic, Dynamic> mass_check = - global_mass;
    for(size_t i = 0 ; i < mass_check.rows() ; i++){
        mass_check(i,i) += global_lumped_mass(i);
        std::cout<<(mass_check.row(i)).sum()<<std::endl;
    }
    */
    
    
    
    /// CHECK OF THE SUM PROPERTY OF CIJ
    /*
     Matrix<T, Dynamic, 1> sum0 = Eigen::Matrix<T, Dynamic, 1>::Zero(global_cij_x.rows(), 1);
     Matrix<T, Dynamic, 1> sum1 = Eigen::Matrix<T, Dynamic, 1>::Zero(global_cij_x.rows(), 1);
    
    for (size_t k = 0; k<global_cij_x.cols(); k++)
    {
        sum0 += (global_cij_x).col(k);
        sum1 += (global_cij_y).col(k);
    }
    
    for (size_t i=0; i<global_cij_y.rows(); i++) {
         std::cout<<"The sum is "<<sum1(i)<<" and "<<sum0(i)<<std::endl;
    }
    
    size_t counter_sum = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
           sum0( counter_sum ) += global_cij_x( counter_sum , elem );
           sum1( counter_sum ) += global_cij_y( counter_sum , elem );
           if(elem==counter_sum)
               std::cout<<"In ("<<counter_sum<<" , "<<elem<<" ), c^0 = "<<global_cij_x( counter_sum , elem )<<" and c1 = "<<global_cij_y( counter_sum , elem )<<std::endl;
           else
                std::cout<<"In ("<<counter_sum<<" , "<<elem<<" ), c^0_ij + c^0_ji = "<<global_cij_x( counter_sum , elem ) + global_cij_x( elem , counter_sum )<<" and c^1_ij + c^1_ji = "<<global_cij_y( counter_sum , elem ) + global_cij_y( elem , counter_sum ) <<std::endl;
        }
        counter_sum++;
    }
    //std::cout<<"The sum0 is "<<'\n'<<sum0<<" and sum1 "<<'\n'<<sum1<<std::endl;
    */
     
   

    // FLUX TERM
    Matrix<T, Dynamic, 1> flux0 = (u.vertices.first).cwiseProduct(phi.vertices) ;
    Matrix<T, Dynamic, 1> flux1 = (u.vertices.second).cwiseProduct(phi.vertices) ;
    
    
    
    // RESOLUTION OF phi_tilde GLOBALE     ( with cij(phi_j - phi_i) )
    
    Matrix<T, Dynamic, 1> vec1 = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    sum_Si(S_i , global_mass , phi.vertices , vec1);
    
    Matrix<T, Dynamic, 1> vec2 = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    //Matrix<T, Dynamic, 1> vec2_bis = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    averaged_sum_Si(S_i , global_cij_x , flux0 , vec2);
    averaged_sum_Si(S_i , global_cij_y , flux1 , vec2);
    //sum_Si(S_i , global_cij_x , flux0 , vec2_bis);
    //sum_Si(S_i , global_cij_y , flux1 , vec2_bis);
    //std::cout<<"diff vec 2 = "<<'\n'<<vec2-vec2_bis<<std::endl;
    /// THERE IS A DIFFERENCE OF  <1E-18 DEPENDING IF I CONSIDER OR NOT THE AVERAGED SUM
    
    Matrix<T, Dynamic, 1> b = vec1 - dt*vec2.cwiseQuotient(global_lumped_mass);
    
    /// OLD IMPLEMENTATION
    //Matrix<T, Dynamic, 1> prova_FAST = global_mass.completeOrthogonalDecomposition().solve(b);
    //phi_tilde.vertices = global_mass.HouseholderQR().solve(b);
    timecounter tc00;
    tc00.tic();
    LLT<Matrix<T, Dynamic, Dynamic > > llt( global_mass );
    phi_tilde.vertices = llt.solve(b);
    tc00.toc();
    
    
    if (llt.info() != Success)
    {
        std::cout<<"not positive"<<std::endl;
        assert(0);
    }
    //std::cout<<"ERRORE FAST= "<<'\n'<<phi_tilde.vertices - prova_FAST <<std::endl;
   
    phi_tilde.converting_into_HHO_formulation( phi_tilde.vertices );
    std::cout << bold << yellow << "LLT RESOLUTION TIME: " << tc00 << " seconds" << reset << std::endl;
    
    
    //std::cout<<"phi_tilde= "<<'\n'<<phi_tilde.vertices<<std::endl;
    /*
    postprocess_output<double> postoutput3;
    typedef typename Mesh::point_type       point_type;
    point<double,2> node;
    auto test_phi_tilde = std::make_shared< gnuplot_output_object<double> >("phi_tilde.dat");
    for (size_t i = 0; i<phi_tilde.vertices.rows(); i++) {
        node = msh.points.at(i);
        test_phi_tilde->add_data(node,phi_tilde.vertices(i));
    }
    postoutput3.add_object(test_phi_tilde);
    postoutput3.write();
    */
    
    
    // Term R_i^n
    Matrix<T, Dynamic, 1> Emax_global = Eigen::Matrix<T, Dynamic, 1>::Ones(dim, 1);
    Matrix<T, Dynamic, 1> Emin_global = Eigen::Matrix<T, Dynamic, 1>::Ones(dim, 1);
    Emax_global *= -1e20;
    Emin_global *= 1e20;
    Matrix<T, Dynamic, 1> R_i = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    for( auto& cl: msh.cells )
    {
        r_i_calculator( msh , cl , E , phi_tilde , phi , dt , u ,S_i , Emax_global , Emin_global , R_i);
    }
    //std::cout<<"Emax_global - Emin_global = "<<'\n'<<Emax_global-Emin_global<<std::endl;
    R_i = R_i.cwiseQuotient( Emax_global - Emin_global );
    //std::cout<<"Ri = "<<'\n'<<R_i<<std::endl;
    
    
    
    // CONVOLUTION TERM
    //sum_Si(S_i, global_cij_x , flux0 , conv_global );
    //sum_Si(S_i, global_cij_y , flux1 , conv_global );
    averaged_sum_Si(S_i, global_cij_x , flux0 , conv_global );
    averaged_sum_Si(S_i, global_cij_y , flux1 , conv_global );
    // OLD IMPLEMENTATION
    //auto conv_global2 = global_cij_x*( flux0 ) + global_cij_y*( flux1 ) ;
    
    // Norm of c_ij
    moltiplication_Si( S_i , global_cij_x , global_cij_x , cij_square );
    moltiplication_Si( S_i , global_cij_y , global_cij_y , cij_square );
    Matrix<T, Dynamic, Dynamic> cij_norm = cij_square.cwiseSqrt();
    
    // Matrix n_ij
    division_Si( S_i , global_cij_x , cij_norm , nij.first );
    division_Si( S_i , global_cij_y , cij_norm , nij.second );
    // OLD IMPLEMENTATION
    //auto nij2 = std::make_pair( global_cij_x.cwiseQuotient(cij_norm) , (global_cij_y).cwiseQuotient( cij_norm ) ) ;
    
    // NOTICE:
    /// Normal Velocity calculation IS EQUIVALENT TO f', that is EQUIVALENT to lambda_r , lambda_l.  Moreover I don't need the calculation of the normal flux. Since then it has to be derived, it is enough to calculate f' = u dot n
     
    // Normal velocity
    Matrix<T, Dynamic, Dynamic> normal_vel = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim, dim);
    size_t counter_row = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            auto u_max_0 = std::max( u.vertices.first(counter_row),u.vertices.first(elem) );
            auto u_max_1 = std::max( u.vertices.second(counter_row),u.vertices.second(elem) );
            normal_vel(counter_row,elem) = u_max_0 * nij.first(counter_row,elem) + u_max_1*nij.second(counter_row,elem);
        }
        counter_row++;
    }
     
    
    
   
    // C_ji matrix
    std::pair<Matrix<T, Dynamic, Dynamic>,Matrix<T, Dynamic, Dynamic>> cji_global = std::make_pair(global_cij_x.adjoint() , global_cij_y.adjoint() );
    
    // Norm of c_ji -> i.e. c_ij transposed
    moltiplication_Si( S_i , cji_global.first , cji_global.first , cji_square );
    moltiplication_Si( S_i , cji_global.second , cji_global.second , cji_square );
    // CHECK IF WORKS BETTER NOW (  NO TRANSPOSE )
    //moltiplication_Si_T( S_i , cji_global.first , cji_global.first , cji_square );
    //moltiplication_Si_T( S_i , cji_global.second , cji_global.second , cji_square );
    
    //C_ji Norm
    Matrix<T, Dynamic, Dynamic> cji_norm = cji_square.cwiseSqrt();
    //auto cji_square2 = (cji_global.first).cwiseProduct(cji_global.first) + (cji_global.second).cwiseProduct(cji_global.second);
    //Matrix<T, Dynamic, Dynamic> cji_norm2 = cji_square2.cwiseSqrt();
    
   
    // Matrix n_ij (TRANSPOSED)
    division_Si( S_i , cji_global.first , cji_norm , nij_transpose.first );
    division_Si( S_i , cji_global.second , cji_norm , nij_transpose.second );
     // CHECK IF WORKS BETTER NOW (  NO TRANSPOSE )
    //division_Si_T( S_i , cji_global.first , cji_norm , nij_transpose.first );
    //division_Si_T( S_i , cji_global.second , cji_norm , nij_transpose.second );
    //auto  nij_transposex2 = cji_global.first.cwiseQuotient( cji_norm2 );
    //auto  nij_transposey2 = cji_global.second.cwiseQuotient( cji_norm2 );
    
    // Normal velocity (TRANSPOSED)
    Matrix<T, Dynamic, Dynamic> normal_vel_adj = Eigen::Matrix<T,Dynamic,Dynamic>::Zero(dim, dim);
    //Matrix<T, Dynamic, Dynamic> normal_vel_adj2 = Eigen::Matrix<T,Dynamic,Dynamic>::Zero(dim, dim);
    /*
    for(size_t i = 0 ; i< normal_vel_adj2.rows() ; i++)
    {
        for(size_t j = 0 ; j< normal_vel_adj2.cols() ; j++)
        {
            auto u_max_0 = std::max( u.vertices.first(i),u.vertices.first(j) );
            auto u_max_1 = std::max( u.vertices.second(i),u.vertices.second(j) );
            //normal_vel(counter_row,elem) = u_max_0 * nij.first(counter_row,elem) + u_max_1*nij.second(counter_row,elem);
            // CHECK IF WORKS BETTER NOW (  NO TRANSPOSE )
            normal_vel_adj2(i,j) = u_max_0 * nij_transposex2(i,j) + u_max_1*nij_transposey2(i,j);
        }
    }
    
    */

    counter_row = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            auto u_max_0 = std::max( u.vertices.first(elem),u.vertices.first(counter_row) );
            auto u_max_1 = std::max( u.vertices.second(elem),u.vertices.second(counter_row) );
            //normal_vel_adj(elem,counter_row) = u_max_0 * nij.first(elem,counter_row) + u_max_1*nij.second(elem,counter_row);
            // CHECK IF WORKS BETTER NOW (  NO TRANSPOSE )
            normal_vel_adj(counter_row , elem) = u_max_0 * nij_transpose.first(counter_row,elem) + u_max_1*nij_transpose.second(counter_row,elem);
        }
        counter_row++;
    }
    
    
    // Matrix Dij calculation
    Matrix<T, Dynamic, Dynamic> lambda_max = normal_vel.cwiseAbs() ;  // o lambda_max=normal_vel ?
    Matrix<T, Dynamic, Dynamic> tmp = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim, dim);
    moltiplication_Si( S_i , lambda_max , cij_norm , tmp );
    
    Matrix<T, Dynamic, Dynamic> lambda_max_adj=normal_vel_adj.cwiseAbs();// o lambda_max=normal_vel?
    Matrix<T, Dynamic, Dynamic> tmp_adj = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim, dim);
    moltiplication_Si( S_i , lambda_max_adj , cji_norm , tmp_adj );
    // CHECK IF WORKS BETTER NOW (  NO TRANSPOSE )
    //moltiplication_Si_T( S_i , lambda_max_adj , cji_norm , tmp_adj );
    //Matrix<T, Dynamic, Dynamic> lambda_max_adj2=normal_vel_adj2.cwiseAbs();
    //auto tmp_adj2 = lambda_max_adj.cwiseProduct( cji_norm2 ) ;
    
    Matrix<T, Dynamic, Dynamic> dij = tmp.cwiseMax(tmp_adj);
    
    //Matrix<T, Dynamic, Dynamic> dij2 = tmp.cwiseMax(tmp_adj2);
    
    //std::cout<<"Dij: "<<'\n'<<dij<<std::endl;
    
    /*
      /// The  error is < 1e-15, at max
    std::cout<<"D_ij check: "<<std::endl;
    for(size_t i = 0; i<tmp.rows();i++){
        for(size_t j = i+1; j<tmp.rows();j++){
            if(std::abs(dij(i,j)-dij(j,i))>1e-15)
            std::cout<<dij(i,j)-dij(j,i);
        }
        std::cout<<'\n';
    }
    std::cout<<std::endl;
    */
   
    
    tc1.toc();
    std::cout << bold << yellow << "Solving time method original: " << tc1 << " seconds" << reset << std::endl;
    
    timecounter tc2;
    tc2.tic();
       
    

    
    // Minimum between d_ij and R_i
    T c_e = 1.;
    Matrix<T, Dynamic, Dynamic> dE_ij = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim, dim);
    Matrix<T, Dynamic, Dynamic> phi_ij = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim, dim);
    counter_row = 0;
   
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            if(elem!=counter_row)
            {
            auto R_i_j = c_e * std::max( std::abs( R_i(counter_row) ) , std::abs( R_i(elem) ) );
            dE_ij( counter_row , elem ) = std::min( dij( counter_row , elem ), R_i_j );
            //std::cout<<"R_i_j "<<R_i_j<<" and dij"<<dij( counter_row , elem )<<std::endl;
            }
            phi_ij( counter_row , elem ) = 0.5*( phi.vertices(counter_row) + phi.vertices(elem) );
        }
        counter_row++;
    }
    //std::cout<<"dE_ij: "<<'\n'<<dE_ij<<std::endl;
    //std::cout<<"phi_ij"<<'\n'<<phi_ij<<std::endl;
    T c_comp = 1.;
    Matrix<T, Dynamic, Dynamic> dC_ij = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim, dim);
    Matrix<T, Dynamic, Dynamic> tmp0 = phi_ij - phi_ij.cwiseProduct(phi_ij);
    //Matrix<T, Dynamic, Dynamic> tmp0 = (phi_ij-phi.phi_min*Mat_One).cwiseProduct(phi.phi_max*Mat_One-phi_ij);
    //std::cout<<"tmp0: "<<'\n'<<tmp0<<std::endl;
    positive_part( tmp0 );
    
    
   // Matrix<T, Dynamic, Dynamic> tmp_dC0 = positive_part( tmp0 );

   
    Matrix<T, Dynamic, Dynamic> tmp_dC1 =  Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim, dim);
    
    size_t counter = 0;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
           // std::cout<<"In ("<<counter<<" , "<<elem<<") tmp_dC0 = "<<tmp_dC0(counter,elem)<<" and phi.ver diff "<<std::abs(phi.vertices(counter) - phi.vertices(elem))<<std::endl;
            if( std::abs(phi.vertices(counter) - phi.vertices(elem))>1e-15 )
            tmp_dC1(counter,elem) = tmp0(counter,elem)/( std::abs(phi.vertices(counter) - phi.vertices(elem)) );
        }
        counter++;
    }
    //std::cout<<'\n'<<"tmp_dC1: "<<'\n'<<tmp_dC1<<std::endl;
    Matrix<T, Dynamic, Dynamic> tmp1 = Mat_One-c_comp*tmp_dC1;
    
    //std::cout<<'\n'<<"tmp1: "<<'\n'<<tmp1<<std::endl;
    positive_part( tmp1  );
    //std::cout<<'\n'<<"tmp1 POSTIVIVE: "<<'\n'<<tmp1<<std::endl;
    
    moltiplication_Si( S_i , dE_ij , tmp1 , dC_ij );
   
    
    for(size_t i = 0 ; i < dC_ij.rows() ; i++)
    {
        dC_ij(i,i) = 0;
        dij(i,i) = 0;
        dC_ij(i,i) = -dC_ij.row(i).sum();
        dij(i,i) = -dij.row(i).sum();
    }
    //dC_ij = dE_ij.cwiseProduct( positive_part( tmp1  ) );
    //   std::cout<<'\n'<<"dC_ij: "<<'\n'<<dC_ij<<std::endl;
    //   std::cout<<'\n'<<"dij: "<<'\n'<<dij<<std::endl;
    
   
    //std::cout<<"dC_ij: "<<'\n'<<dC_ij<<std::endl;
    
    
    /*
    counter_row = 0;
    std::cout<<"DC_ij check symmetry: "<<std::endl;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            std::cout<<dC_ij(counter_row,elem) - dC_ij(elem,counter_row)<<std::endl;
        }
        counter_row++;
        std::cout<<'\n';
    }
    std::cout<<std::endl;
    
    
    counter_row = 0;
    std::cout<<"D_ij check symmetry "<<std::endl;
    for(auto& row_i:S_i)
    {
        for(auto& elem:row_i)
        {
            std::cout<<dij(counter_row,elem) - dij(elem,counter_row)<<std::endl;
        }
        counter_row++;
        std::cout<<'\n';
    }
    std::cout<<std::endl;
    */
    
    
    // Term D_ij(phi_j - phi_i )
    averaged_sum_Si(S_i, dC_ij , phi.vertices , term_dij );
    //std::cout<<"term_dij : "<<'\n'<<term_dij<<std::endl;
    
    
    //averaged_sum_Si(S_i, dE_ij , phi.vertices , term_dij );
    Matrix<T, Dynamic, 1> term_dij_no_entropy = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1 );
    // VERSION WITHOUT ENTROPY CORRECTION
    averaged_sum_Si(S_i, dij , phi.vertices , term_dij_no_entropy );
    
    //std::cout<<"term_dij_no_entropy : "<<'\n'<<term_dij_no_entropy<<std::endl;
    
    // M*phi^n
    Matrix<T, Dynamic, 1> mass_phi_old = Eigen::Matrix<T, Dynamic, 1>::Zero(dim, 1);
    sum_Si( S_i , global_mass , phi.vertices , mass_phi_old );
    
    
    
  /*
    // **************** OLD IMPLEMENTATION **************** //
    
   for (size_t i = 0; i<global_lumped_mass.rows(); i++) {
          global_lumped_mass(i) = global_mass.row(i).sum() ;
      }
   
   
    Matrix<T, Dynamic, Dynamic> normal_vel2 = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim1, dim2);
    Matrix<T, Dynamic, Dynamic> normal_vel_adj2 = Eigen::Matrix<T, Dynamic, Dynamic>::Zero(dim2, dim1);
    
    //normal velocity
    for (size_t i = 0; i < dim1; i++)
    {
        for (size_t j = 0; j < dim2; j++)
        {
            auto u_max_0 = std::max( u.vertices.first(i),u.vertices.first(j) );
            auto u_max_1 = std::max( u.vertices.second(i),u.vertices.second(j) );
            normal_vel2( i , j ) = u_max_0 * nij.first(i,j) + u_max_1*nij.second(i,j);
            //}
        }
           
    }
    
     //normal velocity (transposed)
    auto nij_transpose2 = std::make_pair( ( global_cij_x.transpose() ).cwiseQuotient(cji_norm) ,                                    ( global_cij_y. transpose() ).cwiseQuotient( cji_norm ) ) ;
    
    for (size_t i = 0; i < dim2; i++)
    {
        for (size_t j = 0; j < dim1; j++)
        {
            auto u_max_0 = std::max( u.vertices.first(i),u.vertices.first(j) );
            auto u_max_1 = std::max( u.vertices.second(i),u.vertices.second(j) );
            normal_vel_adj2( i , j ) = u_max_0 * nij_transpose2.first(i,j) + u_max_1*nij_transpose2.second(i,j);
        }
    }

    Matrix<T, Dynamic, Dynamic> lambda_max2 = normal_vel2.cwiseAbs() ;  // o lambda_max=normal_vel ?
    auto tmp2 = lambda_max2.cwiseProduct( cij_norm ) ;
    
    Matrix<T, Dynamic, Dynamic> lambda_max_adj2 = normal_vel_adj2.cwiseAbs() ;
    auto tmp_adj2 = lambda_max_adj2.cwiseProduct( cji_norm ) ;
  
    auto dij2 = tmp2.cwiseMax(tmp_adj2);
    Matrix<T, Dynamic, 1> term_dij2 = Eigen::Matrix<T, Dynamic, 1>::Zero(dim2, 1);
    /// CHECK OUT THE FUNCTION: make_dij_vector
    // auto phi_old = (phi.values_bis).col(counter);
    // auto last_term = make_dij_vector(msh , dij, phi_old );

    for (size_t i = 0; i<dij.rows(); i++)
    {
        auto diff =  phi.vertices - phi.vertices(i)*Vec_One;
        term_dij2(i) = ( dij2.row(i) ).dot( diff )  ;
    }
    std::cout<<"diff term dij: "<<term_dij-term_dij2<<std::endl;
   
    // ****************END OF OLD IMPLEMENTATION **************** //
    
 */
    
    tc2.toc();
    std::cout << bold << yellow << "Solving time method2: " << tc2 << " seconds" << reset << std::endl;
    
    
    timecounter tc3;
    tc3.tic();
       
    ///********* RESOLUTION OF THE SYSTEM: **********//
    //tc.tic();
    Matrix<T, Dynamic, 1> phi_L = solveFEM(global_lumped_mass , conv_global , term_dij_no_entropy , phi.vertices , dt , bdry_nodes); // VERSION WITHOUT ENTROPY CORRECTION
    tc3.toc();
    std::cout << bold << yellow << "Solving time ORIGINAL METHOD: " << tc3 << " seconds" << reset << std::endl;
    
    timecounter tc4;
    tc4.tic();

     Matrix<T, Dynamic, 1> phi_H = solveFEM_Entropic_FAST(llt , conv_global , term_dij , mass_phi_old , dt );
    
    //Matrix<T, Dynamic, 1> phi_H = solveFEM_Entropic(global_mass , conv_global , term_dij , mass_phi_old , dt );
     //std::cout<<'\n'<<"phi_H -phi_H_fast : "<<'\n'<<phi_H-phi_H_prova<<std::endl;
    tc4.toc();
    
    std::cout << bold << yellow << "Solving time CORRECTED METHOD: " << tc4 << " seconds" << reset << std::endl;
    
    
    
    // std::cout<<'\n'<<"phi_L : "<<'\n'<<phi_L<<std::endl;
    
    
    // CHECKING phi_h and phi_l
    Matrix<T, Dynamic, 1> zero_vec = Eigen::Matrix<T, Dynamic, 1>::Zero(phi.vertices.rows(), 1);
   // checking_phi_l( global_lumped_mass , global_cij_x , global_cij_y ,dij , phi.vertices ,phi_L,dt, S_i); // NO, DA MODIFICARE
    
    //checking_phi_lBIS( global_lumped_mass , conv_global , term_dij_no_entropy , phi.vertices , phi_L , dt );
    
    //checking_phi_h( global_mass , global_cij_x , global_cij_y , dC_ij , phi.vertices  , phi_H ,dt,S_i); // NO, DA MODIFICARE
    
    //checking_phi_hBIS( global_mass , conv_global , term_dij , phi.vertices , phi_H , dt );
    
   
    
    
    

    // EXTENSION: MAXIMUM PRINCIPLE PRESERVING
    
    Matrix<T, Dynamic, 1> delta_phi = phi_H - phi.vertices;
    //std::cout<<'\n'<<"delta_phi "<<'\n'<<delta_phi<<std::endl;
    
    Matrix<T, Dynamic, 1> f_i = f_ij_creator( global_lumped_mass , global_mass , delta_phi , dt , dij , dC_ij , phi.vertices , S_i );
    
    // CHECKING PHI
    auto check_phi = phi_H - phi_L - f_i.cwiseQuotient(global_lumped_mass) ;//+dt*Vec_One ;
    Matrix<T, Dynamic, 1>  phi_H2 = phi_L + f_i.cwiseQuotient(global_lumped_mass);
    //std::cout<<'\n'<<"check phi "<<'\n'<<check_phi<<std::endl;
   // std::cout<<'\n'<<"check phi_h: "<<'\n'<<phi_H - phi_H2<<std::endl;
    
    // CORRECTION TERM
    Matrix<T, Dynamic, 1> correction_fi = alfaf_ij_creator( global_lumped_mass , global_mass , delta_phi , phi_L , dt , dij , dC_ij , phi.vertices , S_i );
    Matrix<T, Dynamic, 1>  phi_new = phi_L + correction_fi.cwiseQuotient(global_lumped_mass);
    
    //T check_phi2 = ((phi_new - phi_L).cwiseProduct(global_lumped_mass)).sum();
   // std::cout<<'\n'<<"check PHI_NEW - PHI_L = "<<check_phi2<<std::endl;
    
    //T check_phi3 = ((phi_new - phi.vertices).cwiseProduct(global_lumped_mass)).sum();
   // std::cout<<'\n'<<"check PHI_NEW - PHI_OLD = "<<check_phi3<<std::endl;
    
   // T check_phi4 = ((phi_L - phi.vertices).cwiseProduct(global_lumped_mass)).sum();
   // std::cout<<'\n'<<"check PHI_L - PHI_OLD = "<<check_phi4<<std::endl;
    
    //T check_phi5 = ((phi_H - phi.vertices).cwiseProduct(global_lumped_mass)).sum();
    //std::cout<<'\n'<<"check PHI_H - PHI_OLD = "<<check_phi5<<std::endl;
    
    //T check_phi6 = ((phi_H - phi_L).cwiseProduct(global_lumped_mass)).sum();
    //std::cout<<'\n'<<"check PHI_H - PHI_L = "<<check_phi6<<std::endl;
    
    

    // BOUNDARY CONDITION IMPOSITION -->   NO ENTROPY
    for (auto& j : bdry_nodes )
        phi_L(j) = phi.vertices(j) ;
    
    // BOUNDARY CONDITION IMPOSITION -> ENTROPIC SOLUTION
    for (auto& j : bdry_nodes )
        phi_H(j) = phi.vertices(j) ;
    
    // BOUNDARY CONDITION IMPOSITION -> ENTROPIC AND MAX PRESERVING SOLUTION
    for (auto& j : bdry_nodes )
        phi_new(j) = phi.vertices(j) ;
    
    
    // SAVING AND UPLOAD phi_new  INTO CLASS projected_level_set
    phi.converting_into_HHO_formulation(phi_new);
    phi.vertices = phi_new;
    
    // INVERSE MAPPING BACK [0,1] --> [PHI_MIN,PHI_MAX]
    inverse_mapping_phi( phi , phi.phi_max , phi.phi_min );
    tc.toc();
    std::cout << bold << yellow << "FEM method, time resolution: " << tc << " seconds" << reset << std::endl;
    
    /// PLOTTING SOLUTION (GNUPLOT)
    /*
    postprocess_output<double> postoutput5;
    // auto test_phi_h = std::make_shared< gnuplot_output_object<double> >("phi_h.dat");
    // auto test_phi_l = std::make_shared< gnuplot_output_object<double> >("phi_l.dat");
    auto test_phi_new = std::make_shared< gnuplot_output_object<double> >("phi_new.dat");
    
    size_t iii = 0;
    for (auto pt : msh.points) {
      //  test_phi_h->add_data(pt,phi_H(iii) );
      //  test_phi_l->add_data(pt,phi_L(iii) );
        test_phi_new->add_data(pt,phi_new(iii) );
        iii++;
    }
    //   postoutput5.add_object(test_phi_h);
    //  postoutput5.add_object(test_phi_l);
    postoutput5.add_object(test_phi_new);
    postoutput5.write();
    */
    
    std::cout << bold << yellow << "!!ATTENCTION: I'M COMPUTING BOTH phi AND phi_entropic!!" << reset <<std::endl;
    
    
    
    /// COUT OF PHI IN CELLWISE NOTATION
/*
    for (size_t i = 0; i<phi.values_bis.rows(); i++) {
        for (size_t j = 0; j<phi.values_bis.cols(); j++) {
            std::cout<<phi.values_bis(i,j);
        }
        std::cout<<'\n';
    }
    std::cout<<std::endl;
*/
  
    

}



template< typename Mesh , typename Fonction , typename T = typename Mesh::coordinate_type ,typename VEC >
void Lp_space_Tfin_error_FEM( const Fonction& level_set_final , const Fonction& level_set_initial , const Mesh& msh , size_t degree ,double p , VEC& error )
{
    T errorLp = 0.;
    for(auto& cl : msh.cells)
    {
        auto qps = integrate(msh, cl, 2*degree + 2 );
            for (auto& qp : qps )
            {
                auto diff_val = std::abs( level_set_final( qp.first , msh , cl ) - level_set_initial( qp.first , msh , cl ) );
                errorLp += qp.second * pow(diff_val , p) ;
            }
        //std::cout<<"The L^2 error squared in cell "<<offset(msh,cl)<<" is "<< errorL2 <<std::endl;
    }
    //std::cout<<"The L^"<<p<<" error is "<< pow (errorLp , 1.0/p )<<std::endl;
    error.push_back(pow (errorLp , 1.0/p ));
}



template< typename Mesh , typename Fonction , typename T = typename Mesh::coordinate_type >
void Lp_space_Tfin_error_FEM( const Fonction& level_set_final , const Fonction& level_set_initial , const Mesh& msh , size_t degree ,double p )
{
    T errorLp = 0.;
    for(auto& cl : msh.cells)
    {
        auto qps = integrate(msh, cl, 2*degree + 2 );
            for (auto& qp : qps )
            {
                auto diff_val = std::abs( level_set_final( qp.first , msh , cl ) - level_set_initial( qp.first , msh , cl ) );
                errorLp += qp.second * pow(diff_val , p) ;
            }
        //std::cout<<"The L^2 error squared in cell "<<offset(msh,cl)<<" is "<< errorL2 <<std::endl;
    }
    std::cout<<"The L^"<<p<<" error is "<< pow (errorLp , 1.0/p )<<std::endl;
}

template< typename Mesh , typename Fonction , typename T = typename Mesh::coordinate_type >
T Lp_space_error_FEM( const Fonction& level_set_final , const Fonction& level_set_initial , const Mesh& msh , size_t degree , double p )
{
    // L^p in space ; l^q in time
    T errorLp = 0.;
    for(auto& cl : msh.cells)
    {
        auto qps = integrate(msh, cl, 2*degree + 2 ); // what orders?
            for (auto& qp : qps )
            {
                auto diff_val = std::abs( level_set_final( qp.first , msh , cl ) - level_set_initial( qp.first , msh , cl ) );
                errorLp += qp.second * pow(diff_val,p) ;
            }
        //std::cout<<"The L^2 error squared in cell "<<offset(msh,cl)<<" is "<< errorL2 <<std::endl;
    }
    errorLp = pow( errorLp , 1.0/ p );
    //std::cout<<"The L^2 error is "<<sqrt( errorL2 )<<std::endl;
    return errorLp ;
}

template< typename Mesh , typename Fonction , typename T = typename Mesh::coordinate_type >
T W1p_error_FEM( const Fonction& level_set_final , const Fonction& level_set_initial , const Mesh& msh , size_t degree , double p )
{
    T errorH1 = 0.;
    for(auto& cl : msh.cells)
    {
        auto qps = integrate(msh, cl, 2*degree + 2 ); // what orders?
            for (auto& qp : qps )
            {
                auto diff_val0 = std::abs( level_set_final.gradient( qp.first , msh , cl )(0) - level_set_initial.gradient( qp.first , msh , cl )(0) );
                auto diff_val1 = std::abs( level_set_final.gradient( qp.first , msh , cl )(1) - level_set_initial.gradient( qp.first , msh , cl )(1) );
                errorH1 += qp.second * (pow(diff_val0,p) + pow(diff_val1,p)) ;
            }
        //std::cout<<"The L^2 error squared in cell "<<offset(msh,cl)<<" is "<< errorL2 <<std::endl;
    }
    errorH1 = pow( errorH1 , 1.0/p );
    //std::cout<<"The L^2 error is "<<sqrt( errorL2 )<<std::endl;
    return errorH1 ;
}



template< typename Mesh , typename Fonction , typename T = typename Mesh::coordinate_type >
T Linf_error_FEM( const Fonction& level_set_final , const Fonction& level_set_initial , const Mesh& msh , size_t degree )
{
    T errorLinf = ((level_set_final.vertices - level_set_initial.vertices).cwiseAbs()).maxCoeff();
    return errorLinf;

}

template< typename Mesh , typename Fonction , typename T = typename Mesh::coordinate_type >
T W1inf_error_FEM( const Fonction& level_set_final , const Fonction& level_set_initial , const Mesh& msh , size_t degree )
{
    T errorLinf0 = 0 , errorLinf1 = 0;
    
    for(auto& cl : msh.cells)
    {
        auto pts = points(msh,cl);
        for (auto& pt : pts )
        {
            auto diff_val0 = std::abs( level_set_final.gradient( pt , msh , cl )(0) - level_set_initial.gradient( pt , msh , cl )(0) );
            errorLinf0 = std::max(errorLinf0 , diff_val0);
            auto diff_val1 = std::abs( level_set_final.gradient( pt , msh , cl )(1) - level_set_initial.gradient( pt , msh , cl )(1) );
            errorLinf1 = std::max(errorLinf1 , diff_val1);
        }
        
    }
    return errorLinf0 + errorLinf1;

}


template< typename T , typename VeloField>
T time_step_CFL( const VeloField& u , const mesh_init_params<T>& mip , T eps ){
    
    auto h_max = std::max(mip.hx() , mip.hy() );
    auto u_max = u.values_bis.first.template lpNorm<Infinity>() + u.values_bis.second.template lpNorm<Infinity>() ;
    if( std::abs(u_max) < 1e-15 )
        return 1e-8;
    else
        return eps*h_max/u_max ;
}

template< typename T , typename VeloField>
T time_step_CFL_new( const VeloField& u , const mesh_init_params<T>& mip ,T eps)
{
    
    auto h_max = std::max(mip.hx() , mip.hy() );
    auto u_max = u.sol_FEM.first.template lpNorm<Infinity>() + u.sol_FEM.second.template lpNorm<Infinity>() ;
    if( std::abs(u_max) < 1e-15 )
        return 1e-8;
    else
        return eps*h_max/(std::abs(u_max)) ;
}
















/*****************************************************************************
*   PREVIOUS CODE STOKES HHO
*****************************************************************************/



///////////////////////   FICTITIOUS DOMAIN METHODS  ///////////////////////////

template<typename T, size_t ET, typename testType>
class stokes_fictdom_method
{
    using Mat  = Matrix<T, Dynamic, Dynamic>;
    using Vect = Matrix<T, Dynamic, 1>;
    using Mesh = cuthho_mesh<T, ET>;

protected:
    bool sym_grad;

    stokes_fictdom_method(bool sym)
        : sym_grad(sym)
        {}

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
        Mat gr2;
        if(sym_grad)
            gr2 = make_hho_gradrec_sym_matrix(msh, cl, hdi).second;
        else
            gr2 = make_hho_gradrec_matrix(msh, cl, hdi).second;
        Mat stab = make_hho_vector_naive_stabilization(msh, cl, hdi);
        Mat lc = gr2 + stab;
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

    gradrec_stokes_fictdom_method(T eta_, bool sym)
        : stokes_fictdom_method<T,ET,testType>(sym), eta(eta_) {}

    std::pair< std::pair<Mat,Mat>, std::pair<Vect,Vect> >
    make_contrib_cut(const Mesh& msh, const typename Mesh::cell_type& cl,
                     const testType test_case, const hho_degree_info hdi,
                     const element_location where = element_location::IN_NEGATIVE_SIDE,
                     const params<T>& parms = params<T>())
    {
        // LHS
        Mat gr1, gr2;
        if( this->sym_grad )
        {
            auto gr = make_hho_gradrec_sym_matrix(msh, cl, test_case.level_set_, hdi, where, 1.0);
            gr1 = gr.first;
            gr2 = gr.second;
        }
        else
        {
            auto gr = make_hho_gradrec_matrix(msh, cl, test_case.level_set_, hdi, where, 1.0);
            gr1 = gr.first;
            gr2 = gr.second;
        }
        Mat stab = make_hho_vector_cut_stabilization(msh, cl, hdi, where)
            + make_hho_cut_interface_vector_penalty(msh, cl, hdi, eta);
        Mat lc = gr2 + stab;
        auto dr = make_hho_divergence_reconstruction(msh, cl, test_case.level_set_,
                                                     hdi, where, 1.0);

        // RHS
        auto celdeg = hdi.cell_degree();
        auto cbs = vector_cell_basis<Mesh,T>::size(celdeg);

        Vect f = Vect::Zero(lc.rows());
        f.block(0, 0, cbs, 1) += make_vector_rhs(msh, cl, celdeg, test_case.rhs_fun, where);
        f.block(0, 0, cbs, 1) += make_vector_rhs_penalty(msh, cl, celdeg, test_case.bcs_vel, eta);
        f += make_vector_GR_rhs(msh, cl, celdeg, test_case.bcs_vel, test_case.level_set_,
                                gr1, this->sym_grad);
        auto p_rhs = make_pressure_rhs(msh, cl, hdi.face_degree(), where,
                                       test_case.level_set_, test_case.bcs_vel);

        return std::make_pair(std::make_pair(lc, dr.second), std::make_pair(f,p_rhs) );
    }
};



template<typename T, size_t ET, typename testType>
auto make_gradrec_stokes_fictdom_method(const cuthho_mesh<T, ET>& msh, const T eta_,
                                        const testType test_case, bool sym)
{
    return gradrec_stokes_fictdom_method<T, ET, testType>(eta_, sym);
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
    auto assembler = make_stokes_fict_assembler(msh, bcs_fun, hdi, where);
    auto assembler_sc = make_stokes_fict_condensed_assembler(msh, bcs_fun, hdi, where);

    // method with gradient reconstruction (penalty-free)
    auto class_meth = make_gradrec_stokes_fictdom_method(msh, 1.0, test_case, true);

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
            assembler_sc.assemble(msh, cl, lc_A, lc_B, rhs_A, rhs_B);
        else
            assembler.assemble(msh, cl, lc_A, lc_B, rhs_A, rhs_B);
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

    auto uT_l2_gp  = std::make_shared< gnuplot_output_object<RealType> >("fictdom_uT_norm.dat");
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
            locdata_vel = assembler_sc.take_velocity(msh, cl, sol);
            locdata_p   = assembler_sc.take_pressure(msh, cl, sol);
        }
        else
        {
            locdata_vel = assembler.take_velocity(msh, cl, sol);
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
                uT_l2_gp->add_data( qp.first, std::sqrt( v(0)*v(0) + v(1)*v(1) ) );

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

    postoutput.add_object(uT_l2_gp);
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
    bool sym_grad;

    stokes_interface_method(bool sym)
        : sym_grad(sym) {}

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

        Mat gr2;
        if(sym_grad)
            gr2 = make_hho_gradrec_sym_matrix(msh, cl, hdi).second;
        else
            gr2 = make_hho_gradrec_matrix(msh, cl, hdi).second;
        Mat stab = make_hho_vector_naive_stabilization(msh, cl, hdi);
        Mat lc = kappa * (gr2 + stab);
        auto dr = make_hho_divergence_reconstruction(msh, cl, hdi);
        Mat f = make_vector_rhs(msh, cl, hdi.cell_degree(), test_case.rhs_fun);

        size_t v_size = gr2.rows();
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

    Sym_gradrec_stokes_interface_method(T eta_, T gamma_, bool sym)
        : stokes_interface_method<T,ET,testType>(sym), eta(eta_), gamma_0(gamma_) {}

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
        Mat gr2_n, gr2_p;
        if(this->sym_grad)
        {
            gr2_n = make_hho_gradrec_sym_matrix_interface
                (msh, cl, level_set_function, hdi,element_location::IN_NEGATIVE_SIDE, 0.5).second;
            gr2_p = make_hho_gradrec_sym_matrix_interface
                (msh, cl, level_set_function, hdi,element_location::IN_POSITIVE_SIDE, 0.5).second;
        }
        else
        {
            gr2_n = make_hho_gradrec_matrix_interface
                (msh, cl, level_set_function, hdi,element_location::IN_NEGATIVE_SIDE, 0.5).second;
            gr2_p = make_hho_gradrec_matrix_interface
                (msh, cl, level_set_function, hdi,element_location::IN_POSITIVE_SIDE, 0.5).second;
        }

        // stab
        Mat stab = make_hho_vector_stabilization_interface(msh, cl, level_set_function, hdi,parms);

        Mat penalty = make_hho_cut_interface_vector_penalty(msh, cl, hdi, eta).block(0,0,cbs,cbs);
        stab.block(0, 0, cbs, cbs) += parms.kappa_2 * penalty;
        stab.block(0, cbs, cbs, cbs) -= parms.kappa_2 * penalty;
        stab.block(cbs, 0, cbs, cbs) -= parms.kappa_2 * penalty;
        stab.block(cbs, cbs, cbs, cbs) += parms.kappa_2 * penalty;

        Mat lc = stab + parms.kappa_1 * gr2_n + parms.kappa_2 * gr2_p;

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
                                              const T gamma_, testType test_case, bool sym)
{
    return Sym_gradrec_stokes_interface_method<T, ET, testType>(eta_, gamma_, sym);
}






///////////////////////////////////////

template<typename Mesh, typename testType, typename meth , typename Fonction , typename Velocity>
stokes_test_info<typename Mesh::coordinate_type>
run_cuthho_interface_velocity(const Mesh& msh, size_t degree, meth method, testType test_case , Fonction & level_set_function , Velocity & velocity , bool sym_grad )
{
    using RealType = typename Mesh::coordinate_type;
    
    //auto level_set_function = test_case.level_set_;

    auto bcs_vel = test_case.bcs_vel;
    
    struct params<RealType> parms = test_case.parms;

    timecounter tc;

    bool sc = true; // static condensation


    // ************** ASSEMBLE PROBLEM **************
    hho_degree_info hdi(degree+1, degree);

    tc.tic();
    auto assembler = make_stokes_interface_assembler(msh, bcs_vel, hdi);
    auto assembler_sc = make_stokes_interface_condensed_assembler(msh, bcs_vel, hdi);
    for (auto& cl : msh.cells)
    {
        // ADD BY STE
        level_set_function.cell_assignment(cl);
        //auto test_case_cell = make_test_case_stokes_2(msh, level_set_function);
        auto prm = params<RealType>();
        prm.kappa_1 = 1.0;
        prm.kappa_2 = 1.0;
        auto test_case_cell = make_test_case_kink_velocity2(msh, level_set_function, prm, sym_grad);
        //This stuff before it was before the loop in msh.cells and test_case and not test_case_cell
        auto rhs_fun = test_case_cell.rhs_fun;
        auto sol_vel = test_case_cell.sol_vel;
        auto sol_p = test_case_cell.sol_p;
        auto vel_grad = test_case_cell.vel_grad;
        //auto bcs_vel = test_case.bcs_vel;
        auto neumann_jump = test_case_cell.neumann_jump;
        assembler_sc.set_dir_func( bcs_vel);
        
        
        
        auto contrib = method.make_contrib(msh, cl, test_case_cell, hdi);
        auto lc = contrib.first;
        auto f = contrib.second;

        if( sc )
            assembler_sc.assemble(msh, cl, lc, f);
        else
            assembler.assemble(msh, cl, lc, f);
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

    // ************** SOLVE **************
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

    // ************** POSTPROCESS **************


    postprocess_output<RealType>  postoutput;

    auto uT1_gp  = std::make_shared< gnuplot_output_object<RealType> >("interface_uT1.dat");
    auto uT2_gp  = std::make_shared< gnuplot_output_object<RealType> >("interface_uT2.dat");
    auto p_gp    = std::make_shared< gnuplot_output_object<RealType> >("interface_p.dat");

    tc.tic();
    RealType    H1_error = 0.0;
    RealType    L2_error = 0.0;
    RealType    L2_pressure_error = 0.0;
    
    
    size_t i_global = 0 ; // ADD BY STE
    for (auto& cl : msh.cells)
    {
        vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> cb(msh, cl, hdi.cell_degree());
        cell_basis<cuthho_poly_mesh<RealType>, RealType> pb(msh, cl, hdi.face_degree());
        auto cbs = cb.size();
        auto pbs = pb.size();

        
        // ADD BY STE
        level_set_function.cell_assignment(cl);
        //auto test_case_cell = make_test_case_stokes_2(msh, level_set_function);
        
        auto prm = params<RealType>();
        prm.kappa_1 = 1.0;
        prm.kappa_2 = 1.0;
        auto test_case_cell = make_test_case_kink_velocity2(msh, level_set_function,  prm , sym_grad);
        //This stuff before it was before the loop in msh.cells and test_case and not test_case_cell
        auto rhs_fun = test_case_cell.rhs_fun;
        auto sol_vel = test_case_cell.sol_vel;
        auto sol_p = test_case_cell.sol_p;
        auto vel_grad = test_case_cell.vel_grad;
        //auto bcs_vel = test_case.bcs_vel;
        auto neumann_jump = test_case_cell.neumann_jump;
        assembler_sc.set_dir_func( bcs_vel);
        
        
        Matrix<RealType, Dynamic, 1> vel_locdata_n, vel_locdata_p, vel_locdata;
        Matrix<RealType, Dynamic, 1> P_locdata_n, P_locdata_p, P_locdata;
        Matrix<RealType, Dynamic, 1> vel_cell_dofs_n, vel_cell_dofs_p, vel_cell_dofs;

        if (location(msh, cl) == element_location::ON_INTERFACE)
        {
            if( sc )
            {
                vel_locdata_n = assembler_sc.take_velocity(msh, cl, sol, element_location::IN_NEGATIVE_SIDE);
                vel_locdata_p = assembler_sc.take_velocity(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
                P_locdata_n = assembler_sc.take_pressure(msh, cl, sol, element_location::IN_NEGATIVE_SIDE);
                P_locdata_p = assembler_sc.take_pressure(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
            }
            else
            {
                vel_locdata_n = assembler.take_velocity(msh, cl, sol, element_location::IN_NEGATIVE_SIDE);
                vel_locdata_p = assembler.take_velocity(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
                P_locdata_n = assembler.take_pressure(msh,cl, sol, element_location::IN_NEGATIVE_SIDE);
                P_locdata_p = assembler.take_pressure(msh,cl, sol, element_location::IN_POSITIVE_SIDE);
            }

            vel_cell_dofs_n = vel_locdata_n.head(cbs);
            vel_cell_dofs_p = vel_locdata_p.head(cbs);


            
            // Uploading velocity field by STE
            auto Lagrange_nodes_Qk = equidistriduted_nodes_ordered_bis<RealType,Mesh> (msh,cl,velocity.degree_FEM);
            size_t i_local = 0;
            for ( const auto & ln_Qk : Lagrange_nodes_Qk)
            {
                if( level_set_function(ln_Qk) > 0 )
                {
                    auto phi_HHO = cb.eval_basis( ln_Qk );
                    auto vel = phi_HHO.transpose() * vel_cell_dofs_p;
                    velocity.sol_HHO.first(i_local,i_global) = vel(0);
                    velocity.sol_HHO.second(i_local,i_global) = vel(1);
                }
                else
                {
                    auto phi_HHO = cb.eval_basis( ln_Qk );
                    auto vel = phi_HHO.transpose() * vel_cell_dofs_n;
                    velocity.sol_HHO.first(i_local,i_global) = vel(0);
                    velocity.sol_HHO.second(i_local,i_global) = vel(1);
                    //  velocity.first(i_local,i_global) = cell_dofs_n.dot( phi_HHO );
                    //  velocity.second(i_local,i_global) = 0; // elliptic case is scalar
                }
            }
            
            
            
            auto qps_n = integrate(msh, cl, 2*hdi.cell_degree(), element_location::IN_NEGATIVE_SIDE);
            for (auto& qp : qps_n)
            {
                // Compute H1-error //
                auto t_dphi = cb.eval_gradients( qp.first );
                Matrix<RealType, 2, 2> grad = Matrix<RealType, 2, 2>::Zero();

                for (size_t i = 1; i < cbs; i++ )
                    grad += vel_cell_dofs_n(i) * t_dphi[i].block(0, 0, 2, 2);

                Matrix<RealType, 2, 2> grad_diff = vel_grad(qp.first) - grad;
                H1_error += qp.second * inner_product(grad_diff , grad_diff);


                 // Compute L2-error //
                auto t_phi = cb.eval_basis( qp.first );
                auto v = t_phi.transpose() * vel_cell_dofs_n;
                Matrix<RealType, 2, 1> sol_diff = sol_vel(qp.first) - v;
                L2_error += qp.second * sol_diff.dot(sol_diff);

                uT1_gp->add_data( qp.first, v(0) );
                uT2_gp->add_data( qp.first, v(1) );

                // L2 - pressure - error //
                auto p_phi = pb.eval_basis( qp.first );
                RealType p_num = p_phi.dot(P_locdata_n);
                RealType p_diff = test_case_cell.sol_p( qp.first ) - p_num; // era test_case STE
                L2_pressure_error += qp.second * p_diff * p_diff;

                p_gp->add_data( qp.first, p_num );
            }

            auto qps_p = integrate(msh, cl, 2*hdi.cell_degree(), element_location::IN_POSITIVE_SIDE);
            for (auto& qp : qps_p)
            {
                // Compute H1-error //
                auto t_dphi = cb.eval_gradients( qp.first );
                Matrix<RealType, 2, 2> grad = Matrix<RealType, 2, 2>::Zero();

                for (size_t i = 1; i < cbs; i++ )
                    grad += vel_cell_dofs_p(i) * t_dphi[i].block(0, 0, 2, 2);

                Matrix<RealType, 2, 2> grad_diff = vel_grad(qp.first) - grad;
                H1_error += qp.second * inner_product(grad_diff , grad_diff);

                // Compute L2-error //
                auto t_phi = cb.eval_basis( qp.first );
                auto v = t_phi.transpose() * vel_cell_dofs_p;
                Matrix<RealType, 2, 1> sol_diff = sol_vel(qp.first) - v;
                L2_error += qp.second * sol_diff.dot(sol_diff);

                uT1_gp->add_data( qp.first, v(0) );
                uT2_gp->add_data( qp.first, v(1) );

                 // L2 - pressure - error //
                auto p_phi = pb.eval_basis( qp.first );
                RealType p_num = p_phi.dot(P_locdata_p);
                RealType p_diff = test_case_cell.sol_p( qp.first ) - p_num; // era test_case STE
                L2_pressure_error += qp.second * p_diff * p_diff;

                p_gp->add_data( qp.first, p_num );
            }
        }
        else
        {
            if( sc )
            {
                vel_locdata = assembler_sc.take_velocity(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
                P_locdata = assembler_sc.take_pressure(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
            }
            else
            {
                vel_locdata = assembler.take_velocity(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
                P_locdata = assembler.take_pressure(msh,cl, sol, element_location::IN_POSITIVE_SIDE);
            }
            vel_cell_dofs = vel_locdata.head(cbs);

            
            
            
            auto Lagrange_nodes_Qk = equidistriduted_nodes_ordered_bis<RealType,Mesh> (msh,cl,velocity.degree_FEM);
            size_t i_local = 0;
            for ( const auto & ln_Qk : Lagrange_nodes_Qk)
            {
                auto phi_HHO = cb.eval_basis( ln_Qk );
                auto vel = phi_HHO.transpose() * vel_cell_dofs;
                velocity.sol_HHO.first(i_local,i_global) = vel(0);
                velocity.sol_HHO.second(i_local,i_global) = vel(1);
            }
            
            
            
            
            
            
            auto qps = integrate(msh, cl, 2*hdi.cell_degree());
            for (auto& qp : qps)
            {
                // Compute H1-error //
                auto t_dphi = cb.eval_gradients( qp.first );
                Matrix<RealType, 2, 2> grad = Matrix<RealType, 2, 2>::Zero();

                for (size_t i = 1; i < cbs; i++ )
                    grad += vel_cell_dofs(i) * t_dphi[i].block(0, 0, 2, 2);

                Matrix<RealType, 2, 2> grad_diff = vel_grad(qp.first) - grad;
                H1_error += qp.second * inner_product(grad_diff , grad_diff);

                // Compute L2-error //
                auto t_phi = cb.eval_basis( qp.first );
                auto v = t_phi.transpose() * vel_cell_dofs;
                Matrix<RealType, 2, 1> sol_diff = sol_vel(qp.first) - v;
                L2_error += qp.second * sol_diff.dot(sol_diff);

                uT1_gp->add_data( qp.first, v(0) );
                uT2_gp->add_data( qp.first, v(1) );

                // L2 - pressure - error //
                auto p_phi = pb.eval_basis( qp.first );
                RealType p_num = p_phi.dot(P_locdata);
                RealType p_diff = test_case_cell.sol_p( qp.first ) - p_num; // era test_case STE
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
        // Matrix<RealType, Dynamic, Dynamic> Mat;
        if (sc)
            Mat = assembler_sc.LHS;
        else
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
    auto assembler = make_stokes_interface_assembler(msh, bcs_vel, hdi);
    auto assembler_sc = make_stokes_interface_condensed_assembler(msh, bcs_vel, hdi);
    for (auto& cl : msh.cells)
    {
        auto contrib = method.make_contrib(msh, cl, test_case, hdi);
        auto lc = contrib.first;
        auto f = contrib.second;

        if( sc )
            assembler_sc.assemble(msh, cl, lc, f);
        else
            assembler.assemble(msh, cl, lc, f);
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
            if( sc )
            {
                vel_locdata_n = assembler_sc.take_velocity(msh, cl, sol, element_location::IN_NEGATIVE_SIDE);
                vel_locdata_p = assembler_sc.take_velocity(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
                P_locdata_n = assembler_sc.take_pressure(msh, cl, sol, element_location::IN_NEGATIVE_SIDE);
                P_locdata_p = assembler_sc.take_pressure(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
            }
            else
            {
                vel_locdata_n = assembler.take_velocity(msh, cl, sol, element_location::IN_NEGATIVE_SIDE);
                vel_locdata_p = assembler.take_velocity(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
                P_locdata_n = assembler.take_pressure(msh,cl, sol, element_location::IN_NEGATIVE_SIDE);
                P_locdata_p = assembler.take_pressure(msh,cl, sol, element_location::IN_POSITIVE_SIDE);
            }

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
            if( sc )
            {
                vel_locdata = assembler_sc.take_velocity(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
                P_locdata = assembler_sc.take_pressure(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
            }
            else
            {
                vel_locdata = assembler.take_velocity(msh, cl, sol, element_location::IN_POSITIVE_SIDE);
                P_locdata = assembler.take_pressure(msh,cl, sol, element_location::IN_POSITIVE_SIDE);
            }
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
        // Matrix<RealType, Dynamic, Dynamic> Mat;
        if (sc)
            Mat = assembler_sc.LHS;
        else
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
    // mesh_sizes.push_back(128);
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

            if(1)
            {
                // auto test_case = make_test_case_stokes_1(msh, level_set_function);
                auto test_case = make_test_case_stokes_2(msh, level_set_function);
                TI = run_cuthho_fictdom(msh, k, test_case);
                // auto method = make_sym_gradrec_stokes_interface_method(msh, 1.0, 0.0, test_case, true);
                // TI = run_cuthho_interface(msh, k, method, test_case);
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

#if 0
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

    // auto test_case = make_test_case_stokes_1(msh, level_set_function);
    auto test_case = make_test_case_stokes_2(msh, level_set_function);

    auto method = make_sym_gradrec_stokes_interface_method(msh, 1.0, 0.0, test_case, true);

    if (solve_interface)
        run_cuthho_interface(msh, degree, method, test_case);

    if (solve_fictdom)
        run_cuthho_fictdom(msh, degree, test_case);


    return 0;
}
#endif





// Velocity Field from HHO STOKES
#if 1
int main(int argc, char **argv)
{
    using RealType = double;

    size_t degree           = 0;
    size_t int_refsteps     = 4;
    size_t degree_FEM       = 0;

    bool dump_debug         = false;
    bool solve_interface    = false;
    bool solve_fictdom      = false;
    bool agglomeration      = false;

    mesh_init_params<RealType> mip;
    mip.Nx = 5;
    mip.Ny = 5;
    RealType d = 0.5;
    size_t T_N = 0;
    /* k <deg>:     method degree
     * g<deg>:  method FEM degree
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
    while ( (ch = getopt(argc, argv, "k:q:c:M:N:r:T:ifDAd")) != -1 )
    {
        switch(ch)
        {
            case 'k':
                degree = atoi(optarg);
                break;

            case 'q':
                degree_FEM = atoi(optarg);
                break;
                
            case 'c':
                d = 1./atoi(optarg);
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
            
            case 'T':
                T_N = atoi(optarg);
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
    typedef cuthho_poly_mesh<RealType> Mesh;
    offset_definition(msh);
    tc.toc();
    std::cout << bold << yellow << "Mesh generation: " << tc << " seconds" << reset << std::endl;
       
    
    
    
    /************** FINITE ELEMENT INITIALIZATION **************/
    auto fe_data = Finite_Element<RealType,Mesh>( msh , degree_FEM , mip ) ;
     typedef Finite_Element<RealType,Mesh> FiniteSpace;
    
    /************** ANALYTIC LEVEL SET FUNCTION  **************/
    RealType radius = 1.0/9.0; // I PUT 6.0, IT WAS 1.0/3.0
    RealType x_centre = 0.45;
    RealType y_centre = 0.5;
    auto level_set_function_anal = circle_level_set<RealType>(radius, x_centre, y_centre);
    //auto level_set_function = circle_level_set<RealType>(radius, 0.5, 0.5);
    //auto level_set_function = flower_level_set<RealType>(0.31, 0.5, 0.5, 4, 0.04);
    // auto level_set_function = line_level_set<RealType>(0.5);

     std::cout << bold << yellow << "Initial Analytic Area Circle: "<< M_PI*radius*radius << " , and circonference: "<< 2*M_PI*radius << reset << std::endl;
    
    
    // auto level_set_function_anal = flower_level_set<RealType>(0.31, 0.5, 0.5, 4, 0.04);//(0.11, 0.1, 0.1, 4, 0.04);
    //  auto level_set_function_anal=square_level_set<RealType>(0.6,0.4,0.4,0.6);
    
    typedef RealType T;
    typedef  circle_level_set<T> Fonction;
    //typedef  flower_level_set<T> Fonction; //DA RIMETTERE POI PER LA DISCRETA
    //typedef line_level_set<T> Fonction;
    // typedef square_level_set<T> Fonction;
    
     /**************  VELOCITY FIELD   **************/ // JUST TO CHECK IMPLEMENTATION
    /*
    T u_0 = 1.00 ;
    T u_1 = 0.00 ;
    auto u = linear_velocity_field<RealType>(0,u_0,0,u_1); // analytic velocity (0,0,0,0)
    typedef linear_velocity_field<RealType> Velocity;
    std::cout<<yellow<<bold<<"--> USING LAGRANGIAN HIGH ORDER VELOCITY FIELD"<<reset<<std::endl;
    auto u_projected = projection_velocity_high_order <Mesh,Velocity,FiniteSpace,T> (fe_data , u , msh);
    */
    
    
    /************** LEVEL SET FUNCTION DISCRETISATION **************/
    std::cout<<"degree FEM "<<degree_FEM<<std::endl;
    T C = 0.2;  // Imposed max value at the boundary
    //T r0 = radius + (x_centre - radius + 2.0/16.0)/2.0;
    //T C = r0*r0 - radius*radius ;
    
    auto level_set_function = L2projected_level_set_high_order< Mesh , Fonction , FiniteSpace , T > (fe_data , level_set_function_anal , msh);
    
    //level_set_function.smooth_cut_off( C , x_centre , y_centre , radius );
    level_set_function.cut_off( C );
    
    testing_level_set(msh,level_set_function,level_set_function_anal);
     
    // Initiliatisation data for time routine
    auto crr_mesh =  Current_Mesh<Mesh>(msh);
    
    /// Initialisation area , mass -> I calculate them just the first and last time step.
    T initial_area = 0. , initial_mass = 0.;
    T area_previous_time = 0. , mass_previous_time = 0. , dt = 0. ;
    
    /// DATA CHECK INITIALISATION
    T d_a = 0. ; // It is related to area-> calculated just first and last time
    T error_normal_global = 0. ;
    T error_normal_local = 0. ;
    T perimeter_initial = 0. ;
    T centre_mass_x_inital = 0. , centre_mass_y_inital = 0. ;
    
    
    /// ERROR FEM (transport problem) wrt ANALYTICAL SOLUTION
    //T l1_L1_error = 0 , l2_L1_error = 0 , linf_L1_error = 0;
    //T l1_L2_error = 0 ,  l2_L2_error = 0 , linf_L2_error = 0;
    //T l1_Linf_error = 0 ,  l2_Linf_error = 0 , linf_Linf_error = 0;
    
    //T l1_W11_error = 0 , l2_W11_error = 0 , linf_W11_error = 0;
    //T l1_W12_error = 0 ,  l2_W12_error = 0 , linf_W12_error = 0;
    //T l1_W1inf_error = 0 ,  l2_W1inf_error = 0 , linf_W1inf_error = 0;
    
    
    for (size_t time_step = 0; time_step<=T_N; time_step++)
    {
    
        // ************** Re-Initialization mesh **************
        crr_mesh.current_mesh = msh;
        Mesh msh_i =  crr_mesh.current_mesh;
        offset_definition(msh_i);
        


    
        
        //************ DO cutHHO MESH PROCESSING **************
        tc.tic();
        detect_node_position2(msh_i, level_set_function); // In cuthho_geom
        detect_cut_faces2(msh_i, level_set_function); // In cuthho_geom
        
        if (agglomeration)
        {
            detect_cut_cells2(msh_i, level_set_function); // In cuthho_geom
            detect_cell_agglo_set(msh_i, level_set_function); // Non serve modificarla
            make_neighbors_info_cartesian(msh_i); // Non serve modificarla
            refine_interface2(msh_i, level_set_function, int_refsteps); // IN cuthho_geom
            make_agglomeration(msh_i, level_set_function); // Non serve modificarla
        }
        else
        {
            //move_nodes(msh_i, level_set_function);
            //detect_cut_faces2(msh_i, level_set_function); //do it again to update intersection points
            detect_cut_cells2(msh_i, level_set_function);
            refine_interface2(msh_i, level_set_function, int_refsteps);
        }
       
        tc.toc();
        std::cout << bold << yellow << "cutHHO-specific mesh preprocessing: " << tc << " seconds" << reset << std::endl;

        if (dump_debug)
        {
            dump_mesh(msh_i);
            output_mesh_info(msh_i, level_set_function);
        }
   
        // IN cuthho_export..Points/Nodes don't change-> it's fast
        if(time_step == 0)
            output_mesh_info2_pre_FEM(msh_i, level_set_function); // IN cuthho_export
            //output_mesh_info2(msh_i, level_set_function);
        
        typedef L2projected_level_set_high_order< Mesh , Fonction , FiniteSpace , T > Level_Set;
        auto ls_cell = LS_cell_L2proj_high_order< T , Mesh , Level_Set, Fonction , FiniteSpace >(level_set_function,msh_i);
        
        // CALCULATION OF AREA AND MASS AT TIME STEP t^n
        // CALCULATION ALSO OF CENTRE OF MASS
        
        /// DATA CHECK INITIALISATION
        T area0 = 0. , mass0 = 0. , global_mass0 = 0. ;
        T diff_area = 0. , diff_mass = 0. ;
        T error_normal_global0 = 0. ;
        T diff_interface_normal = 0. ;
        size_t counter_interface_pts = 0;
        T divergence_error_inital = 0. ;
        T centre_mass0_x = 0. , centre_mass0_y = 0. ;
        T perimeter0 = 0.;
        postprocess_output<double> postoutput_div2;
        auto test_interface_divergence2  = std::make_shared< gnuplot_output_object<double> >("divergence_interface_initial.dat");
        
        
        for(auto& cl : msh_i.cells)
        {
            ls_cell.cell_assignment(cl);
            if( location(msh_i, cl) == element_location::IN_NEGATIVE_SIDE || location(msh_i, cl) == element_location::ON_INTERFACE )
            {
                T partial_area = measure( msh_i, cl, element_location::IN_NEGATIVE_SIDE);
                area0 += partial_area;
                auto qps = integrate( msh_i , cl , 2*degree_FEM+1 , element_location::IN_NEGATIVE_SIDE);
                for(auto& qp:qps){
                    mass0 += qp.second * ls_cell(qp.first);
                    centre_mass0_x += qp.second * qp.first.x() ;
                    centre_mass0_y += qp.second * qp.first.y() ;
                }
            }
            if(cl.user_data.location == element_location::ON_INTERFACE)
            {
                // NORMAL ERROR ON THE INTERFACE POINTS: diff_interface_normal
                for(auto& interface_point : cl.user_data.interface)
                {
                    counter_interface_pts++;
                    T val = ls_cell.divergence( interface_point );
                    divergence_error_inital += pow((std::abs(val) - 1.0/radius),2) ;
                    std::cout<<"Divergence value in pt = "<<interface_point<<" is "<<val<<std::endl;
                    test_interface_divergence2->add_data(interface_point,val);
                }
                   
                for(auto interface_point = cl.user_data.interface.begin() ; interface_point < cl.user_data.interface.end() -1 ; interface_point++ )
                {
                    perimeter0 += ( *(interface_point+1) - *interface_point ).to_vector().norm();
                       
                }
                   
               
            
            }
            
        }
        
        postoutput_div2.add_object(test_interface_divergence2);
        postoutput_div2.write();
        
        divergence_error_inital /= counter_interface_pts;
        divergence_error_inital = sqrt(divergence_error_inital);
        std::cout<<yellow<<bold<<"The l2 error of the DIVERGENCE at the INTERFACE, at time "<<reset<< dt*time_step <<" is " << divergence_error_inital <<std::endl;
        std::cout<<"number of interface points is "<<counter_interface_pts<<std::endl;
        
        
        std::cout << bold << yellow << "Area at time step: "<<time_step*dt<<" is "<<reset<< area0  << reset << std::endl;
        std::cout << bold << yellow << "Internal mass at time step: "<<time_step*dt<<" is "<<reset<< mass0  << reset << std::endl;
        //std::cout<<bold<<yellow << "GLOBAL Mass at time step: "<<time_step<<" is "<<reset<<global_mass0<< reset << std::endl;
         std::cout << bold << yellow << "CENTRE OF MASS at time step: " <<reset <<time_step*dt<<" is "<<" ( "<< centre_mass0_x/area0  << " , "<< centre_mass0_y/area0 <<" ). " << reset << std::endl;
        
         d_a = sqrt(4.0*area0/M_PI) ;
        
        std::cout<<yellow<<bold<<"The PERIMETER, at time "<<reset<< dt*time_step <<" is " << perimeter0 <<std::endl;
        
        std::cout<<yellow<<bold<<"The CIRCULARITY, at time "<<reset<< dt*time_step <<" is " << M_PI*d_a/perimeter0 <<std::endl;
        
        
        
        if(time_step == 0 ){
            initial_area  = area0 ;
            initial_mass = mass0 ;
            centre_mass_x_inital = centre_mass0_x ;
            centre_mass_y_inital = centre_mass0_y ;
            perimeter_initial = perimeter0 ;
        }
        
        if(time_step > 0 )
        {
            diff_area = (area0 - initial_area)/initial_area ;
            diff_mass = (std::abs((mass0 - initial_mass)))/(std::abs(initial_mass)) ;
            std::cout << bold << yellow << "Normalised difference in Area (new - old)/old at time step: "<<time_step*dt<<" is "<<reset<< diff_area  << reset << std::endl;
            std::cout << bold << yellow << "Difference in internal MASS |new - old|/|old| at time step: "<<time_step*dt<<" is "<<reset<< diff_mass  << reset << std::endl;
            std::cout << bold << yellow << "TRANSLATION OF THE CENTRE OF MASS at time step: " <<reset <<time_step*dt<<" is "<<" ( " << centre_mass0_x/area0 - centre_mass_x_inital/initial_area <<" , " << centre_mass0_y/area0 - centre_mass_y_inital/initial_area<<" ). " << std::endl;
           
            std::cout<< bold << yellow<<"NORMALISED DIFFERENCE PERIMETER, at time "<<reset<< dt*(time_step) <<" is " << (perimeter0 - perimeter_initial)/perimeter_initial <<std::endl;
                   
        }
        
        // auto test_case = make_test_case_stokes_1(msh, level_set_function);
        // auto test_case = make_test_case_stokes_2(msh, ls_cell); //level_set_function);
        
/// DA AGGIUNGERE UNA VOLTA SISTEMATO IL CODICE

        bool sym_grad = TRUE;
        
        auto prm = params<T>();
        prm.kappa_1 = 1.0;
        prm.kappa_2 = 1.0;
        
        auto test_case = make_test_case_kink_velocity2(msh, ls_cell,  prm , sym_grad);
        
        auto method = make_sym_gradrec_stokes_interface_method(msh, 1.0, 0.0, test_case, sym_grad);

        
        auto u_projected = velocity_high_order <Mesh,FiniteSpace,T> (fe_data , msh);
       
        
    
        if(solve_interface){
            run_cuthho_interface_velocity(msh_i, degree, method,test_case, ls_cell , u_projected ,sym_grad );
            //run_cuthho_interface(msh, degree, method, test_case);
        }
        if (solve_fictdom)
            run_cuthho_fictdom(msh_i, degree, test_case);
        

        
        
        /*********************** FEM -  PROCESSING **************************/
        /// ORA HO SMOOTH OPERATOR! USE L2 PROJECTION
        std::cout<<yellow<<bold<<"------------------>>>>NOTICE: TO PASS FROM HHO TO FEM IM USING SMOOTH OPERATOR. IMPLEMENT L2 PROJECTION!!!"<<reset<<std::endl;
        u_projected.converting_into_FE_formulation(u_projected.sol_HHO);
        
        T eps = 1.0 ; // factor to be inside CFL stability zone
        //T dt1 = time_step_CFL( u , mip , eps ); // OLD IMPLEMENTATION
        T dt1 = time_step_CFL_new( u_projected , mip , eps );
        //std::cout<<"dt1 is "<<dt1<<std::endl;
        T dt2 = 1*1e-3;
        dt = std::min(dt1 , dt2);
        std::cout<<"I USE dt = "<<dt<<" AND HEURISTIC CFL IS "<<dt1<<std::endl;
        
        
         run_FEM_BERNSTEIN_CORRECT( level_set_function.msh , fe_data , level_set_function , u_projected , dt);
        
      
        /// OLD IMPLEMENTATION
        //for(size_t j=1; j<4 ; j++)
        //run_FEM_levelset( level_set_function.msh,degree_FEM,level_set_function,u,dt,mip);
       
        /************** FEM -  POST-PROCESSING **************/
        if( (T_N - time_step)==0 )
        {
            // Uploading mesh data to check out differences in mass and areas
            crr_mesh.current_mesh = msh;
            Mesh msh_i2 =  crr_mesh.current_mesh;
            offset_definition(msh_i2);
            tc.tic();
            detect_node_position2(msh_i2, level_set_function); // In cuthho_geom
            detect_cut_faces2(msh_i2, level_set_function); // In cuthho_geom
        
       
            if (agglomeration)
            {
                // std::cout<<"i m here 1"<<std::endl;
                detect_cut_cells2(msh_i2, level_set_function); // In cuthho_geom
           
                detect_cell_agglo_set(msh_i2, level_set_function); // Non serve modificarla
           
                make_neighbors_info_cartesian(msh_i2); // Non serve modificarla
                refine_interface2(msh_i2, level_set_function, int_refsteps); // IN cuthho_geom
                make_agglomeration(msh_i2, level_set_function); // Non serve modificarla
            }
            else
            {
                //move_nodes(msh_i2, level_set_function);
                //detect_cut_faces2(msh_i2, level_set_function); //do it again to update intersection     points
                detect_cut_cells2(msh_i2, level_set_function);
                refine_interface2(msh_i2, level_set_function, int_refsteps);
            }
        
            tc.toc();
            std::cout << bold << yellow << "cutHHO-specific mesh preprocessing: " << tc << " seconds" << reset << std::endl;

            if (dump_debug)
            {
                dump_mesh(msh_i2);
                output_mesh_info(msh_i2, level_set_function);
            }
        
            output_mesh_info2(msh_i2, level_set_function); // IN cuthho_export..Points/Nodes don't change
           
            // Uploading level set
            ls_cell.level_set = level_set_function;
            ls_cell.agglo_msh = msh_i2;
        
            /// DATA CHECK INITIALISATION
            postprocess_output<double> postoutput_div;
            auto test_interface_divergence  = std::make_shared< gnuplot_output_object<double> >("divergence_interface_u_const.dat");
           
            T mass_fin = 0. , area_fin = 0. ;
            T centre_mass_x = 0. , centre_mass_y = 0. ;
            counter_interface_pts = 0;
            T perimeter = 0. ;
            T divergence_error_fin = 0. ;
            
            
            for(auto& cl : msh_i2.cells)
            {
                ls_cell.cell_assignment(cl);
                
                if( (location(msh_i2, cl) == element_location::IN_NEGATIVE_SIDE) || (location(msh_i2, cl) == element_location::ON_INTERFACE) )
                {
                    
                    T partial_area = measure( msh_i2, cl, element_location::IN_NEGATIVE_SIDE);
                 
                    area_fin += partial_area;
                    
                   
                    auto qps_fin = integrate( msh_i2 , cl , 2*degree_FEM+1 , element_location::IN_NEGATIVE_SIDE);
                   
                    for(auto& qp:qps_fin){
                        mass_fin += qp.second * ls_cell(qp.first);
                        centre_mass_x += qp.second * qp.first.x() ;
                        centre_mass_y += qp.second * qp.first.y() ;
                    }
                   
                }
                if(cl.user_data.location == element_location::ON_INTERFACE)
                {
                    for(auto& interface_point : cl.user_data.interface)
                    {
                        counter_interface_pts++;
                        T val = ls_cell.divergence( interface_point );
                        std::cout<<"Divergence value in pt = "<<interface_point<<" is "<<val<<std::endl;
                        
                        divergence_error_fin += pow((std::abs(val) - 1.0/radius),2);
                    
                        test_interface_divergence->add_data(interface_point,val);
                     
                    }
                
                
                    
                    for(auto interface_point = cl.user_data.interface.begin() ; interface_point < cl.user_data.interface.end() -1 ; interface_point++ )
                    {
                        perimeter += ( *(interface_point+1) - *interface_point ).to_vector().norm();
                        
                    }
                }
                
            }
            
            postoutput_div.add_object(test_interface_divergence);
            postoutput_div.write();
           
            divergence_error_fin /= counter_interface_pts;
            divergence_error_fin = sqrt(divergence_error_fin);
            std::cout<<yellow<<bold<<"The l2 error of the DIVERGENCE at the INTERFACE, at time "<<reset<< dt*time_step <<" is " << divergence_error_fin <<std::endl;
            
            std::cout<<"number of interface points is "<<counter_interface_pts<<std::endl;
            
            std::cout<< bold << yellow<<"The PERIMETER, at time "<<reset<< dt*(time_step+1) <<" is " << perimeter <<std::endl;
            
            std::cout<<"perimeter = "<< perimeter << " AND  perimeter0 =  "<<perimeter_initial<<std::endl;
            std::cout<< bold << yellow<<"NORMALISED DIFFERENCE PERIMETER, at time "<<reset<< dt*(time_step+1) <<" is " << (perimeter - perimeter_initial)/perimeter_initial <<std::endl;
            
            d_a = sqrt(4.0*area_fin/M_PI) ;
            
             std::cout<< bold << yellow<<"The CIRCULARITY, at time "<< dt*(time_step+1)<<reset <<" is " << M_PI*d_a/perimeter <<std::endl;
            
            std::cout << bold << yellow << "Area at time step: " <<(time_step+1)*dt<<" is "<< reset<< area_fin << std::endl;
            std::cout << bold << yellow << "Internal mass at time step: "<<(time_step+1)*dt<<" is "<<reset<< mass_fin << reset << std::endl;
            
            std::cout << bold << yellow << "NORMALISED Difference in AREA AT TIME "<<(T_N+1)*(dt)<<" IS "<< reset<< (area_fin - initial_area)/initial_area << std::endl;
            std::cout << bold << yellow << "NORMALISED Difference in INTERNAL MASS AT TIME "<<(T_N+1)*(dt)<<" IS "<< reset<< (std::abs(mass_fin - initial_mass))/(std::abs( initial_mass )) << std::endl;
            std::cout << bold << yellow << "CENTRE OF MASS at time step: "<<(T_N+1)*dt<<" is "<<reset<<" ( " << centre_mass_x/area_fin <<" , " << centre_mass_y/area_fin<<" ). " << std::endl;
            std::cout << bold << yellow << "TRANSLATION OF THE CENTRE OF MASS at time step: " <<reset <<(T_N+1)*dt<<" is "<<" ( " << centre_mass_x/area_fin - centre_mass_x_inital/initial_area <<" , " << centre_mass_y/area_fin - centre_mass_y_inital/initial_area<<" ). " << std::endl;
            
            /// ERROR CHECK WITH ANALYTICAL SOLUTION
            /*
            T x_deviation = u_0*(time_step+1)*dt;
            T y_deviation = u_1*(time_step+1)*dt;
            auto analytic_level_set_post_FE = circle_level_set<RealType>(radius, x_centre + x_deviation, y_centre + y_deviation );
            auto level_set_final = L2projected_level_set_high_order < Mesh , Fonction , FiniteSpace , T > (fe_data , analytic_level_set_post_FE , msh);
            level_set_final.cut_off( C );
            testing_level_set2(msh,level_set_function,level_set_final);
            T L1_err_par = Lp_space_error_FEM( level_set_final , level_set_function,msh,degree_FEM,1.0);
                   
            T L2_err_par = Lp_space_error_FEM( level_set_final , level_set_function,msh,degree_FEM,2.0);
                    
            T Linf_err_par = Linf_error_FEM( level_set_final , level_set_function , msh,degree_FEM);
                    
            std::cout << bold << yellow << "L1-space error at time "<<(time_step+1)*(dt)<<" IS "<< reset<< L1_err_par << std::endl;
            std::cout << bold << yellow << "L2-space error at time "<<(time_step+1)*(dt)<<" IS "<< reset<< L2_err_par << std::endl;
            std::cout << bold << yellow << "Linf-space error at time "<<(time_step+1)*(dt)<<" IS "<< reset<< Linf_err_par << std::endl;
            
            */
            
           
        } // END OF T = FINAL TIME
       
        
    } // End of the temporal loop
    
    
    return 0;
}
#endif
