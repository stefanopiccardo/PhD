/*
 *       /\        Guillaume Delay (C) 2019
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



using namespace Eigen;

#include "core/core"
#include "core/solvers"
#include "dataio/silo_io.hpp"

#include "methods/hho"
#include "methods/cuthho"


////////////// TEST ON THE INTEGRATION METHOD
//// write in a file the list of the barycenter used to integrate over cut cells
template<typename Mesh, typename Function>
void
test_bar_integration(Mesh& msh, const Function& level_set_function)
{
    // open files
    std::ofstream neg_file("bar_neg.3D", std::ios::out | std::ios::trunc);
    std::ofstream pos_file("bar_pos.3D", std::ios::out | std::ios::trunc);

    if( !(neg_file && pos_file) )
        std::cerr << "bar_file has not been opened" << std::endl;

    neg_file << "X   Y   Z   val" << std::endl;
    pos_file << "X   Y   Z   val" << std::endl;

    for (auto cl : msh.cells)
    {
        if( location(msh, cl) != element_location::ON_INTERFACE )
            continue;

        // compute the barycenter used for negative and positive integrations
        // neg
        auto tp_n = collect_triangulation_points(msh, cl, element_location::IN_NEGATIVE_SIDE );
        // auto bar_n = barycenter(tp_n);
        auto bar_n = tesselation_center(msh, cl, element_location::IN_NEGATIVE_SIDE);
        neg_file << bar_n[0] << "   " <<  bar_n[1] << "   0.0     0.0" << std::endl;

        // pos
        auto tp_p = collect_triangulation_points(msh, cl, element_location::IN_POSITIVE_SIDE );
        // auto bar_p = barycenter(tp_p);
        auto bar_p = tesselation_center(msh, cl, element_location::IN_POSITIVE_SIDE);
        pos_file << bar_p[0] << "   " <<  bar_p[1] << "   0.0     0.0" << std::endl;
    }

    // close files
    neg_file.close();
    pos_file.close();
}


////////////// TEST ON THE INTEGRATION METHOD
//// compare usual integration and composite integration
void
test_comp_integration()
{
    using T = double;

    // use a 10x10 mesh with a square interface
    mesh_init_params<T> mip;
    mip.Nx = 10;
    mip.Ny = 10;
    cuthho_poly_mesh<T> msh(mip);
    size_t int_refsteps     = 1;

    auto level_set_function = square_level_set<T>(0.77, 0.23, 0.23, 0.77);

    detect_node_position(msh, level_set_function);
    detect_cut_faces(msh, level_set_function);
    detect_cut_cells(msh, level_set_function);
    detect_cell_agglo_set(msh, level_set_function);
    make_neighbors_info_cartesian(msh);
    refine_interface(msh, level_set_function, int_refsteps);

    output_mesh_info(msh, level_set_function);

    auto cl1 = msh.cells[22];
    auto cl2 = msh.cells[23];
    auto cl3 = msh.cells[32];
    auto cl4 = msh.cells[33];

    auto where = element_location::IN_NEGATIVE_SIDE;
    size_t k = 0;

    for(size_t i = 0; i < 2; i++)
    {
        if(i == 0)
        {
            std::cout << yellow << "negative side" << std::endl;
            where = element_location::IN_NEGATIVE_SIDE;
        }
        else
        {
            std::cout << blue << "positive side" << std::endl;
            where = element_location::IN_POSITIVE_SIDE;
        }

        // compute area of cl1 + cl2
        std::cout << "cl1" << std::endl;
        T area = 0.0;
        auto qps_n1 = integrate(msh, cl1, k, where);
        for (auto& qp : qps_n1)
        {
            area += qp.second;

            std::cout << qp.first << "   " << qp.second << std::endl;
        }
        std::cout << "cl2" << std::endl;
        auto qps_n2 = integrate(msh, cl2, k, where);
        for (auto& qp : qps_n2)
        {
            area += qp.second;

            std::cout << qp.first << "   " << qp.second << std::endl;
        }
        std::cout << "cl3" << std::endl;
        auto qps_n3 = integrate(msh, cl3, k, where);
        for (auto& qp : qps_n3)
        {
            area += qp.second;

            std::cout << qp.first << "   " << qp.second << std::endl;
        }
        std::cout << "cl4" << std::endl;
        auto qps_n4 = integrate(msh, cl4, k, where);
        for (auto& qp : qps_n4)
        {
            area += qp.second;

            std::cout << qp.first << "   " << qp.second << std::endl;
        }
        std::cout << "area = " << area << std::endl;

        // merge cells and compute the new area
        auto MC1 = merge_cells(msh, cl1, cl2);
        auto MC2 = merge_cells(msh, MC1.first, cl3);
        auto MC = merge_cells(msh, MC2.first, cl4);
        auto cl = MC.first;
        area = 0.0;
        std::cout << "cl" << std::endl;
        auto qps = integrate(msh, cl, k, where);
        for (auto& qp : qps)
        {
            area += qp.second;

            std::cout << qp.first << "   " << qp.second << std::endl;
        }
        std::cout << "area after merging cells = " << area << std::endl;
    }
}


//////////////  TEST_AGGLO
/// test the merging procedure
// this routine has been written for a 10x10 mesh
template<typename Mesh, typename Function>
void
test_agglo(Mesh& msh, const Function& level_set_function)
{
    // merge cell 0 with itself -> OK
    // merge_cells(msh, msh.cells[0], msh.cells[0]);

    // merge cells 0 and 2 -> OK
    // merge_cells(msh, msh.cells[0], msh.cells[2]);
    
    // merge cells 0 and 1
    typename Mesh::cell_type cl1 = msh.cells[0];
    typename Mesh::cell_type cl2 = msh.cells[1];
    merge_cells(msh, cl1, cl2);

    cl1 = msh.cells[33];
    cl2 = msh.cells[43];
    merge_cells(msh, cl1, cl2);


    cl1 = msh.cells[26];
    cl2 = msh.cells[36];
    merge_cells(msh, cl1, cl2);

    cl1 = msh.cells[73];
    cl2 = msh.cells[72];
    merge_cells(msh, cl1, cl2);

    cl1 = msh.cells[0];
    cl2 = msh.cells[1];
    merge_cells(msh, cl1, cl2);

    cl1 = msh.cells[0];
    cl2 = msh.cells[9];
    merge_cells(msh, cl1, cl2);
    
    output_mesh_info(msh, level_set_function);

    plot_quadrature_points(msh,3);
}


//////////////  PLOT_CELLS
/// plot the cells with gnuplot
//
template<typename Mesh, typename Function>
void
plot_cells(Mesh& msh, const Function& level_set_function)
{
    using T = typename Mesh::coordinate_type;

    postprocess_output<T>  postoutput;
    auto cells_gp  = std::make_shared< gnuplot_output_object<T> >("cells.dat");

    size_t count = 0;

    for (auto cl : msh.cells)
    {
        auto qps = integrate(msh, cl, 1);
        for (auto& qp : qps)
        {
            cells_gp->add_data( qp.first, count );
        }
        count++;
    }
    postoutput.add_object(cells_gp);
    postoutput.write();
}


//////////////  TEST_AGGLO_DIAG
/// test the diagonal merging procedure
// this routine has been written for a 10x10 mesh
template<typename Mesh, typename Function>
void
test_agglo_diag(Mesh& msh, const Function& level_set_function)
{
    std::vector<typename Mesh::cell_type> new_cells, removed_cells;

    typename Mesh::cell_type cl1 = msh.cells[1];
    typename Mesh::cell_type cl2 = msh.cells[12];
    auto cl = merge_cells(msh, cl1, cl2).first;
    cl.user_data.highlight = true;
    new_cells.push_back( cl );
    removed_cells.push_back( cl1 );
    removed_cells.push_back( cl2 );

    cl1 = msh.cells[76];
    cl2 = msh.cells[65];
    cl = merge_cells(msh, cl1, cl2).first;
    cl.user_data.highlight = true;
    new_cells.push_back( cl );
    removed_cells.push_back( cl1 );
    removed_cells.push_back( cl2 );

    // remove the agglomerated cells
    typename std::vector<typename Mesh::cell_type>::iterator it_RC;
    for(it_RC = removed_cells.begin(); it_RC != removed_cells.end(); it_RC++) {
        msh.cells.erase(std::remove(begin(msh.cells), end(msh.cells), *it_RC ), end(msh.cells));
    }

    // add new cells
    typename std::vector<typename Mesh::cell_type>::iterator it_NC;
    for(it_NC = new_cells.begin(); it_NC != new_cells.end(); it_NC++) {
        msh.cells.push_back(*it_NC);
    }

    // sort the new list of cells
    std::sort(msh.cells.begin(), msh.cells.end());

    // output the mesh obtained
    output_mesh_info(msh, level_set_function);
}

#if 0
int main(int argc, char **argv)
{
    test_comp_integration();
    return 0;
}
#endif
#if 1
int main(int argc, char **argv)
{
    using RealType = double;
    
    size_t int_refsteps     = 4;

    bool agglomeration      = true;

    mesh_init_params<RealType> mip;
    mip.Nx = 5;
    mip.Ny = 5;

    /* M <num>:     number of cells in x direction
     * N <num>:     number of cells in y direction
     * r <num>:     number of interface refinement steps
     *
     * D:           no agglomeration
     * A:           use agglomeration to solve bad cuts 
     */

    int ch;
    while ( (ch = getopt(argc, argv, "M:N:r:ifDAd")) != -1 )
    {
        switch(ch)
        {
            case 'M':
                mip.Nx = atoi(optarg);
                break;

            case 'N':
                mip.Ny = atoi(optarg);
                break;

            case 'r':
                int_refsteps = atoi(optarg);
                break;

            case 'D':
                agglomeration = false;
                break;

            case 'A':
                agglomeration = true;
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
    // auto level_set_function = circle_level_set<RealType>(radius, 0.5, 0.5);
    // auto level_set_function = line_level_set<RealType>(0.5);
    // auto level_set_function = square_level_set<RealType>(1.0, 0.0, 0.0, 1.0);
    auto level_set_function = square_level_set<RealType>(0.76, 0.24, 0.24, 0.76);
    // auto level_set_function = square_level_set<RealType>(1.05, -0.05, -0.05, 1.05);
    // auto level_set_function = flower_level_set<RealType>(radius, 0.5-0.03, 0.5-0.04, 12, 0.015);
    // auto level_set_function = flower_level_set<RealType>(0.31, 0.5, 0.5, 4, 0.04);
    /************** DO cutHHO MESH PROCESSING **************/

    tc.tic();
    detect_node_position(msh, level_set_function);
    detect_cut_faces(msh, level_set_function);
    detect_cut_cells(msh, level_set_function);

    if (agglomeration)
    {
        detect_cell_agglo_set(msh, level_set_function);
        make_neighbors_info_cartesian(msh);
        // make_neighbors_info(msh);
        refine_interface(msh, level_set_function, int_refsteps);
        // test_agglo(msh, level_set_function);
        // test_agglo_diag(msh, level_set_function);
        // output_mesh_info(msh, level_set_function);
        make_agglomeration(msh, level_set_function);
        // output_mesh_info(msh, level_set_function);
        // test_bar_integration(msh, level_set_function);
        // plot_cells(msh, level_set_function);
    }
    else
    {
        move_nodes(msh, level_set_function);
        detect_cut_faces(msh, level_set_function); //do it again to update intersection points
        detect_cut_cells(msh, level_set_function);
        refine_interface(msh, level_set_function, int_refsteps);
    }

    output_mesh_info(msh, level_set_function);
    
    tc.toc();
    std::cout << bold << yellow << "cutHHO-specific mesh preprocessing: " << tc << " seconds" << reset << std::endl;


    return 0;
}
#endif
