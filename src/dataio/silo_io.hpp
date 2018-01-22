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

#include <unistd.h>
#include <silo.h>

enum variable_centering_t
{
    nodal_variable_t,
    zonal_variable_t
};


class silo_database
{
    DBfile          *m_siloDb;

public:
    silo_database()
        : m_siloDb(nullptr)
    {}

    bool create(const std::string& db_name)
    {
        m_siloDb = DBCreate(db_name.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
        if (m_siloDb)
            return true;

        std::cout << "Error creating database" << std::endl;
        return false;
    }

    bool open(const std::string& db_name)
    {
        m_siloDb = DBOpen(db_name.c_str(), DB_PDB, DB_APPEND);
        if (m_siloDb)
            return true;

        std::cout << "Error opening database" << std::endl;
        return false;
    }

    bool close()
    {
        if (m_siloDb)
            DBClose(m_siloDb);
        m_siloDb = NULL;
        return true;
    }

    ~silo_database()
    {
        if (m_siloDb)
            DBClose(m_siloDb);
    }

    template<typename Mesh>
    bool add_mesh(const Mesh& msh, const std::string& name)
    {
        using T = typename Mesh::coordinate_type;
        static_assert(std::is_same<T, double>::value, "Only double for now");

        std::vector<T> x_coords, y_coords;
        x_coords.reserve(msh.points.size());
        y_coords.reserve(msh.points.size());

        for (auto itor = msh.points.begin(); itor != msh.points.end(); itor++)
        {
            auto pt = *itor;
            x_coords.push_back(pt.x());
            y_coords.push_back(pt.y());
        }

        T *coords[] = {x_coords.data(), y_coords.data()};

        std::vector<int>    nodelist;

        for (auto& cl : msh.cells)
        {
            for (auto& ptid : cl.ptids)
                nodelist.push_back(ptid+1);
        }

        int lnodelist       = nodelist.size();
        int nshapetypes     = 1;
        int shapesize       = 4;
        int shapecount      = msh.cells.size();


        int nnodes = msh.points.size();
        int nzones = msh.cells.size();
        int ndims = 2;

        DBPutZonelist(m_siloDb, "zonelist", nzones, ndims, nodelist.data(), lnodelist,
            1, &shapesize, &shapecount, nshapetypes);

        DBPutUcdmesh(m_siloDb, name.c_str(), ndims, NULL, coords, nnodes, nzones,
            "zonelist", NULL, DB_DOUBLE, NULL);

        return true;
    }

    template<typename T>
    bool add_variable(const std::string& mesh_name,
                      const std::string& var_name,
                      T* data,
                      size_t data_len,
                      variable_centering_t centering)
    {
        static_assert(std::is_same<T,double>::value, "Only double for now");

        if (!m_siloDb)
        {
            std::cout << "Silo database not opened" << std::endl;
            return false;
        }

        if (centering == zonal_variable_t)
        {
            DBPutUcdvar1(m_siloDb, var_name.c_str(), mesh_name.c_str(),
                         data, data_len, NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
            return true;
        }

        if (centering == nodal_variable_t)
        {
            DBPutUcdvar1(m_siloDb, var_name.c_str(), mesh_name.c_str(),
                         data, data_len, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
            return true;
        }

        return false;
    }

};
