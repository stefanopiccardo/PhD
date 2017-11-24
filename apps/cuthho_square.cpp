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

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/SparseExtra>

#include <unistd.h>

#include <silo.h>

using namespace Eigen;

#include "point.hpp"

template<typename T>
struct interface_box {

};


class degree_info {

    size_t cell_deg, face_deg;

public:
    degree_info()
        : cell_deg(1), face_deg(1)
    {}

    degree_info(size_t k, int l)
    {
        if ( l < -1 || l > 1)
            l = 0;

        if (l == -1 && k == 0)
            l = 0;

        cell_deg = k+l;
        face_deg = k;
    }

    size_t cell_degree() const
    {
        return cell_deg;
    }

    size_t face_degree() const
    {
        return face_deg;
    }
};

template<typename T>
bool
conjugated_gradient(const Eigen::SparseMatrix<T>& A,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
                    Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
{
    size_t                      N = A.cols();
    size_t                      iter = 0;
    T                           nr, nr0;
    T                           alpha, beta, rho;

    Eigen::Matrix<T, Eigen::Dynamic, 1> d(N), r(N), r0(N), y(N);
    x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(N);


    r0 = d = r = b - A*x;
    nr = nr0 = r.norm();

    std::ofstream ofs("cg_nopre_convergence.txt");

    while ( nr/nr0 > 1e-8 && iter < 40000 && nr/nr0 < 10000 )
    {
        std::cout << "                                                 \r";
        std::cout << " -> Iteration " << iter << ", rr = ";
        std::cout << nr/nr0 << "\b\r";
        std::cout.flush();

        ofs << nr/nr0 << std::endl;
        y = A*d;
        rho = r.dot(r);
        alpha = rho/d.dot(y);
        x = x + alpha * d;
        r = r - alpha * y;
        beta = r.dot(r)/rho;
        d = r + beta * d;

        nr = r.norm();
        iter++;
    }

    ofs << nr/nr0 << std::endl;
    ofs.close();

    std::cout << " -> Iteration " << iter << ", rr = " << nr/nr0 << std::endl;

    return true;
}

/*****************************************************************************
 *   Mesh stuff
 *****************************************************************************/
struct cell {
    std::array<size_t, 4>   ptids;

    bool cut_by_interface;

    cell()
        : cut_by_interface(false)
    {}

    bool operator<(const cell& other) const
    {
        return (this->ptids < other.ptids);
    }

    bool operator==(const cell& other) const
    {
        return (this->ptids == other.ptids);
    }
};

std::ostream&
operator<<(std::ostream& os, const cell& cl)
{
    os << "Cell: " << cl.ptids[0] << " " << cl.ptids[1];
    os << " " << cl.ptids[2] << " " << cl.ptids[3];
    return os;
}

enum class boundary
{
    NONE,
    DIRICHLET,
    NEUMANN,
    ROBIN
};

struct face {
    std::array<size_t, 2>   ptids;
    bool                    is_boundary;
    boundary                bndtype;

    face()
        : is_boundary(false),
          bndtype(boundary::NONE)
    {}

    bool operator<(const face& other) const
    {
        return (this->ptids < other.ptids);
    }

    bool operator==(const face& other) const
    {
        return (this->ptids == other.ptids);
    }
};

std::ostream&
operator<<(std::ostream& os, const face& fc)
{
    os << "Face: " << fc.ptids[0] << " " << fc.ptids[1];
    if (fc.is_boundary && fc.bndtype == boundary::NONE)
        os << " [Boundary, unspecified kind]";
    if (fc.is_boundary && fc.bndtype == boundary::DIRICHLET)
        os << " [Dirichlet boundary]";
    if (fc.is_boundary && fc.bndtype == boundary::NEUMANN)
        os << " [Neumann boundary]";
    if (fc.is_boundary && fc.bndtype == boundary::ROBIN)
        os << " [Robin boundary]";
    return os;
}

template<typename T>
struct mesh_init_params {
    T           min_x, max_x;
    T           min_y, max_y;
    size_t      Nx, Ny;

    mesh_init_params()
        : min_x(0.0), max_x(1.0),
          min_y(0.0), max_y(1.0),
          Nx(4), Ny(4)
    {}

    T hx() const {
        return (max_x - min_x)/Nx;
    }

    T hy() const {
        return (max_y - min_y)/Ny;
    }
};

template<typename T>
struct mesh {

    typedef point<T,2>          point_type;
    typedef cell                cell_type;
    typedef face                face_type;

    std::vector<point_type>     points;
    std::vector<face_type>      faces;
    std::vector<cell_type>      cells;

    mesh() : mesh( mesh_init_params<T>() )
    {}

    mesh(const mesh_init_params<T>& parms)
    {
        auto hx = parms.hx();
        auto hy = parms.hy();

        size_t numpoints = (parms.Nx + 1) * (parms.Ny + 1);
        points.reserve(numpoints);

        for (size_t j = 0; j < parms.Ny+1; j++)
        {
            for(size_t i = 0; i < parms.Nx+1; i++)
            {
                auto px = parms.min_x + i*hx;
                auto py = parms.min_y + j*hy;
                point_type pt(px, py);
                points.push_back(pt);
            }
        }

        for (size_t j = 0; j < parms.Ny; j++)
        {
            for (size_t i = 0; i < parms.Nx; i++)
            {
                auto pt0_idx = j*(parms.Nx+1) + i;
                auto pt1_idx = pt0_idx + 1;
                auto pt2_idx = pt0_idx + parms.Nx + 2;
                auto pt3_idx = pt0_idx + parms.Nx + 1;

                cell cl;
                cl.ptids = {{pt0_idx, pt1_idx, pt2_idx, pt3_idx}};
                cells.push_back(cl);

                face f0;
                f0.ptids = {pt0_idx, pt1_idx};
                if (j == 0) f0.is_boundary = true;
                faces.push_back(f0);

                face f1;
                f1.ptids = {pt1_idx, pt2_idx};
                if (i == parms.Nx-1) f1.is_boundary = true;
                faces.push_back(f1);

                face f2;
                f2.ptids = {pt3_idx, pt2_idx};
                if (j == parms.Ny-1) f2.is_boundary = true;
                faces.push_back(f2);

                face f3;
                f3.ptids = {pt0_idx, pt3_idx};
                if (i == 0) f3.is_boundary = true;
                faces.push_back(f3);

            }
        }

        std::sort(cells.begin(), cells.end());
        std::sort(faces.begin(), faces.end());
        faces.erase( std::unique(faces.begin(), faces.end()), faces.end() );
    }
};

/*****************************************************************************
 *   SILO stuff
 *****************************************************************************/

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

    template<typename T>
    bool add_mesh(const mesh<T>& msh, const std::string& name)
    {
        static_assert(std::is_same<T,double>::value, "Only double for now");

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
                      const T* data,
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


/*****************************************************************************
 *   Geometry
 *****************************************************************************/
template<typename T>
size_t
offset(const mesh<T>& msh, const typename mesh<T>::cell_type cl)
{
    auto itor = std::lower_bound(msh.cells.begin(), msh.cells.end(), cl);
    if ( itor == msh.cells.end() )
        throw std::logic_error("Cell not found: this is likely a bug.");

    return std::distance(msh.cells.begin(), itor);
}

template<typename T>
size_t
offset(const mesh<T>& msh, const typename mesh<T>::face_type fc)
{
    auto itor = std::lower_bound(msh.faces.begin(), msh.faces.end(), fc);
    if ( itor == msh.faces.end() )
        throw std::logic_error("Face not found: this is likely a bug.");

    return std::distance(msh.faces.begin(), itor);
}

template<typename T>
std::array< typename mesh<T>::point_type, 4 >
points(const mesh<T>& msh, const typename mesh<T>::cell_type cl)
{
    std::array< typename mesh<T>::point_type, 4 > ret;

    auto ptid2pt = [&](size_t ptid) -> auto {
        return msh.points.at(ptid);
    };

    std::transform( cl.ptids.begin(), cl.ptids.end(), ret.begin(), ptid2pt );

    return ret;
}

template<typename T>
std::array< typename mesh<T>::point_type, 2 >
points(const mesh<T>& msh, const typename mesh<T>::face_type fc)
{
    std::array< typename mesh<T>::point_type, 2 > ret;

    auto ptid2pt = [&](size_t ptid) -> auto {
        return msh.points.at(ptid);
    };

    std::transform( fc.ptids.begin(), fc.ptids.end(), ret.begin(), ptid2pt );

    return ret;
}

template<typename T>
std::array< typename mesh<T>::face_type, 4 >
faces(const mesh<T>& msh, const typename mesh<T>::cell_type cl)
{
    typedef typename mesh<T>::face_type face_type;
    std::array< typename mesh<T>::face_type, 4 > ret;

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

template<typename T>
typename mesh<T>::point_type
barycenter(const mesh<T>& msh, const typename mesh<T>::cell_type& cl)
{
    auto pts = points(msh, cl);
    return (pts[0] + pts[1] + pts[2] + pts[3])/4.0;
}

template<typename T>
typename mesh<T>::point_type
barycenter(const mesh<T>& msh, const typename mesh<T>::face_type& fc)
{
    auto pts = points(msh, fc);
    return (pts[0] + pts[1])/2.0;
}

template<typename T>
T
diameter(const mesh<T>& msh, const typename mesh<T>::cell_type& cl)
{
    auto p0 = msh.points.at( cl.ptids[0] );
    auto p1 = msh.points.at( cl.ptids[2] );

    return (p1-p0).to_vector().norm();
}

template<typename T>
T
diameter(const mesh<T>& msh, const typename mesh<T>::face_type& fc)
{
    auto p0 = msh.points.at( fc.ptids[0] );
    auto p1 = msh.points.at( fc.ptids[1] );

    return (p1-p0).to_vector().norm();
}

template<typename T>
T
measure(const mesh<T>& msh, const typename mesh<T>::cell_type& cl)
{
    auto p0 = msh.points.at( cl.ptids[0] );
    auto p2 = msh.points.at( cl.ptids[2] );

    return (p2.x() - p0.x()) * (p2.y() - p0.y());
}

template<typename T>
T
measure(const mesh<T>& msh, const typename mesh<T>::face_type& fc)
{
    auto p0 = msh.points.at( fc.ptids[0] );
    auto p1 = msh.points.at( fc.ptids[1] );

    return (p1-p0).to_vector().norm();
}

template<typename T>
std::array< Matrix<T,2,1>, 4 >
normals(const mesh<T>& msh, const typename mesh<T>::cell_type& cl)
{
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
std::vector< typename mesh<T>::point_type >
make_test_points(const mesh<T>& msh, const typename mesh<T>::cell_type& cl)
{
    const size_t N = 10;

    std::vector< typename mesh<T>::point_type > ret;

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
            typename mesh<T>::point_type pt(px, py);
            ret.push_back(pt);
        }
    }

    return ret;
}

template<typename T>
std::vector< typename mesh<T>::point_type >
make_test_points(const mesh<T>& msh, const typename mesh<T>::face_type& fc)
{
    const size_t N = 10;

    std::vector< typename mesh<T>::point_type > ret;

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

template<typename T>
T iexp_pow(T x, size_t n)
{
    if (n == 0)
        return 1;

    T y = 1;
    while (n > 1)
    {
        if (n % 2 == 0)
        {
            x = x * x;
            n = n / 2;
        }
        else
        {
            y = x * y;
            x = x * x;
            n = (n-1)/2;
        }
    }

    return x*y;
}

/*****************************************************************************
 *   Bases
 *****************************************************************************/
template<typename CT, typename VT>
class cell_basis
{
    typedef typename mesh<CT>::point_type   point_type;

    point_type  cell_bar;
    CT          cell_h;
    size_t      basis_degree, basis_size;

public:
    cell_basis(const mesh<CT>& msh, const typename mesh<CT>::cell_type& cl, size_t degree)
    {
        cell_bar        = barycenter(msh, cl);
        cell_h          = diameter(msh, cl);
        basis_degree    = degree;
        basis_size      = (basis_degree+2)*(basis_degree+1)/2;
    }

    Matrix<VT, Dynamic, 1>
    eval_basis(const point_type& pt)
    {
        Matrix<VT, Dynamic, 1> ret = Matrix<VT, Dynamic, 1>::Zero(basis_size);

        auto bx = (pt.x() - cell_bar.x()) / (0.5*cell_h);
        auto by = (pt.y() - cell_bar.y()) / (0.5*cell_h);

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t i = 0; i <= k; i++)
            {
                auto pow_x = k-i;
                auto pow_y = i;
                auto bv = iexp_pow(bx, pow_x) * iexp_pow(by, pow_y);
                ret(pos++) = bv;
            }
        }

        assert(pos == basis_size);

        return ret;
    }

    Matrix<VT, Dynamic, 2>
    eval_gradients(const point_type& pt)
    {
        Matrix<VT, Dynamic, 2> ret = Matrix<VT, Dynamic, 2>::Zero(basis_size, 2);

        auto bx = (pt.x() - cell_bar.x()) / (0.5*cell_h);
        auto by = (pt.y() - cell_bar.y()) / (0.5*cell_h);
        auto ih = 2.0/cell_h;

        size_t pos = 0;
        for (size_t k = 0; k <= basis_degree; k++)
        {
            for (size_t i = 0; i <= k; i++)
            {
                auto pow_x = k-i;
                auto pow_y = i;
                auto px = iexp_pow(bx, pow_x);
                auto py = iexp_pow(by, pow_y);
                auto dx = (pow_x == 0) ? 0 : pow_x*ih*iexp_pow(bx, pow_x-1);
                auto dy = (pow_y == 0) ? 0 : pow_y*ih*iexp_pow(by, pow_y-1);

                ret(pos,0) = dx*py;
                ret(pos,1) = px*dy;
                pos++;
            }
        }

        assert(pos == basis_size);

        return ret;
    }

    size_t size() const
    {
        return basis_size;
    }

    static size_t size(size_t degree)
    {
        return (degree+2)*(degree+1)/2;
    }
};

#if 0
template<template<typename, typename> class Basis, typename CT, typename VT>
class trial_function
{
    const Basis<CT, VT>& basis;

public:
    trial_function(const Basis<CT, VT>& b)
        : basis(b)
    {}
};

template<template<typename, typename> class Basis, typename CT, typename VT>
auto trial(const Basis<CT, VT>& b)
{
    return trial_function<Basis, CT, VT>(b);
}

template<template<typename, typename> class Basis, typename CT, typename VT>
class test_function
{
    const Basis<CT, VT>& basis;

public:
    test_function(const Basis<CT, VT>& b)
        : basis(b)
    {}
};

template<template<typename, typename> class Basis, typename CT, typename VT>
auto test(const Basis<CT, VT>& b)
{
    return test_function<Basis, CT, VT>(b);
}

template<template <typename, typename> class BasisL, template <typename, typename> class BasisR, typename CT, typename VT>
int
operator,(const trial_function<BasisL,CT,VT> trial, const test_function<BasisR,CT,VT>& test)
{
    std::cout << "\\o/  \\o/  \\o/" << std::endl;
    return 42;
}
#endif

template<typename CT, typename VT>
class face_basis
{
    typedef typename mesh<CT>::point_type   point_type;

    point_type  face_bar;
    point_type  base;
    CT          face_h;
    size_t      basis_degree, basis_size;

public:
    face_basis(const mesh<CT>& msh, const typename mesh<CT>::face_type& fc, size_t degree)
    {
        face_bar        = barycenter(msh, fc);
        face_h          = diameter(msh, fc);
        basis_degree    = degree;
        basis_size      = degree+1;

        auto pts = points(msh, fc);
        base = face_bar - pts[0];
    }

    Matrix<VT, Dynamic, 1>
    eval_basis(const point_type& pt)
    {
        Matrix<VT, Dynamic, 1> ret = Matrix<VT, Dynamic, 1>::Zero(basis_size);

        auto v = base.to_vector();
        auto t = (pt - face_bar).to_vector();
        auto dot = v.dot(t);
        auto ep = 4.0*dot/(face_h*face_h);

        for (size_t i = 0; i <= basis_degree; i++)
        {
            auto bv = iexp_pow(ep, i);
            ret(i) = bv;
        }
        return ret;
    }

    size_t size() const
    {
        return basis_size;
    }

    static size_t size(size_t degree)
    {
        return degree+1;
    }
};

/*****************************************************************************
 *   Quadrature rules
 *****************************************************************************/
template<typename T>
std::vector<std::pair<point<T,1>, T>>
edge_quadrature(size_t doe)
{
    std::vector<std::pair<point<T,1>, T>> ret;

    using namespace Eigen;

    if (doe%2 == 0)
        doe++;

    size_t num_nodes = (doe+1)/2;

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

template<typename T>
std::vector<std::pair<point<T,2>, T>>
integrate(const mesh<T>& msh, const typename mesh<T>::cell_type& cl, size_t degree)
{
    typedef typename mesh<T>::point_type    point_type;

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

template<typename T>
std::vector<std::pair<point<T,2>, T>>
integrate(const mesh<T>& msh, const typename mesh<T>::face_type& fc, size_t degree)
{
    typedef typename mesh<T>::point_type    point_type;

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

/*****************************************************************************
 *   Utilities
 *****************************************************************************/
template<typename T>
Matrix<T, Dynamic, Dynamic>
make_mass_matrix(const mesh<T>& msh, const typename mesh<T>::cell_type& cl, size_t degree)
{
    cell_basis<T,T> cb(msh, cl, degree);
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

template<typename T>
Matrix<T, Dynamic, Dynamic>
make_mass_matrix(const mesh<T>& msh, const typename mesh<T>::face_type& fc, size_t degree)
{
    face_basis<T,T> fb(msh, fc, degree);
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

template<typename T, typename Function>
Matrix<T, Dynamic, 1>
make_rhs(const mesh<T>& msh, const typename mesh<T>::cell_type& cl,
         size_t degree, const Function& f)
{
    cell_basis<T,T> cb(msh, cl, degree);
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

template<typename T, typename Function>
Matrix<T, Dynamic, 1>
make_rhs(const mesh<T>& msh, const typename mesh<T>::face_type& fc,
         size_t degree, const Function& f)
{
    face_basis<T,T> fb(msh, fc, degree);
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

template<typename T, typename Function>
Matrix<T, Dynamic, 1>
project_function(const mesh<T>& msh, const typename mesh<T>::cell_type& cl,
                 size_t degree, const Function& f)
{
    auto cbs = cell_basis<T,T>::size(degree);
    auto fbs = face_basis<T,T>::size(degree);

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

/*****************************************************************************
 *   HHO stuff
 *****************************************************************************/
template<typename T>
std::pair<Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>>
make_hho_laplacian(const mesh<T>& msh, const typename mesh<T>::cell_type& cl, size_t degree)
{
    cell_basis<T,T>     cb(msh, cl, degree+1);

    auto rbs = cell_basis<T,T>::size(degree+1);
    auto cbs = cell_basis<T,T>::size(degree);
    auto fbs = face_basis<T,T>::size(degree);

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
        face_basis<T,T> fb(msh, fc, degree);
        auto qps_f = integrate(msh, fc, 2*degree);
        for (auto& qp : qps_f)
        {
            Matrix<T, Dynamic, 1> c_phi = cb.eval_basis(qp.first);
            c_phi = c_phi.head(cbs);
            Matrix<T, Dynamic, 2> c_dphi = cb.eval_gradients(qp.first);
            c_dphi = c_dphi.block(1, 0, rbs-1, 2);
            Matrix<T, Dynamic, 1> f_phi = fb.eval_basis(qp.first);
            gr_rhs.block(0, cbs+i*fbs, rbs-1, fbs) += qp.second * (c_dphi * n) * f_phi.transpose();
            gr_rhs.block(0, 0, rbs-1, cbs) -= qp.second * (c_dphi * n) * c_phi.transpose();
        }
    }

    Matrix<T, Dynamic, Dynamic> oper = gr_lhs.llt().solve(gr_rhs);
    Matrix<T, Dynamic, Dynamic> data = gr_rhs.transpose() * oper;

    return std::make_pair(oper, data);
}

template<typename T>
std::pair<Matrix<T, Dynamic, Dynamic>, Matrix<T, Dynamic, Dynamic>>
make_operator(const mesh<T>& msh, const typename mesh<T>::cell_type& cl, size_t degree)
{
    return make_hho_laplacian(msh, cl, degree);
}

template<typename T>
Matrix<T, Dynamic, Dynamic>
make_hho_naive_stabilization(const mesh<T>& msh, const typename mesh<T>::cell_type& cl, size_t degree)
{
    auto cbs = cell_basis<T,T>::size(degree);
    auto fbs = face_basis<T,T>::size(degree);

    auto fcs = faces(msh, cl);

    Matrix<T, Dynamic, Dynamic> data = Matrix<T, Dynamic, Dynamic>::Zero(cbs+4*fbs, cbs+4*fbs);
    Matrix<T, Dynamic, Dynamic> If = Matrix<T, Dynamic, Dynamic>::Identity(fbs, fbs);

    cell_basis<T,T> cb(msh, cl, degree);
    auto h = measure(msh, cl);

    for (size_t i = 0; i < 4; i++)
    {
        auto fc = fcs[i];
        face_basis<T,T> fb(msh, fc, degree);

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

        data += oper.transpose() * mass * oper / h;
    }

    return data;
}

template<typename T>
Matrix<T, Dynamic, Dynamic>
make_stabilization(const mesh<T>& msh, const typename mesh<T>::cell_type& cl, size_t degree)
{
    return make_hho_naive_stabilization(msh, cl, degree);
}

/*****************************************************************************
 *   Test stuff
 *****************************************************************************/
template<typename T>
void
plot_basis_functions(const mesh<T>& msh)
{
    std::ofstream c_ofs("cell_basis_check.dat");

    for (auto cl : msh.cells)
    {
        cell_basis<double, double> cb(msh, cl, 2);

        auto tps = make_test_points(msh, cl);

        for (auto& tp : tps)
        {
            c_ofs << tp.x() << " " << tp.y() << " ";

            auto vals = cb.eval_basis(tp);
            for(size_t i = 0; i < cb.size(); i++)
                c_ofs << vals(i) << " ";

            c_ofs << std::endl;
        }
    }

    c_ofs.close();

    std::ofstream f_ofs("face_basis_check.dat");

    for (auto fc : msh.faces)
    {
        face_basis<double, double> fb(msh, fc, 2);

        auto tps = make_test_points(msh, fc);

        for (auto& tp : tps)
        {
            f_ofs << tp.x() << " " << tp.y() << " ";

            auto vals = fb.eval_basis(tp);
            for(size_t i = 0; i < fb.size(); i++)
                f_ofs << vals(i) << " ";

            f_ofs << std::endl;
        }
    }

    f_ofs.close();
}

template<typename T>
void
plot_quadrature_points(const mesh<T>& msh, size_t degree)
{
    std::ofstream c_ofs("cell_quadrature_check.dat");

    for (auto& cl : msh.cells)
    {
        auto qps = integrate(msh, cl, degree);

        for (auto& qp : qps)
        {
            c_ofs << qp.first.x() << " " << qp.first.y();
            c_ofs << " " << qp.second << std::endl;
        }
    }

    c_ofs.close();

    std::ofstream f_ofs("face_quadrature_check.dat");

    for (auto& fc : msh.faces)
    {
        auto qps = integrate(msh, fc, degree);

        for (auto& qp : qps)
        {
            f_ofs << qp.first.x() << " " << qp.first.y();
            f_ofs << " " << qp.second << std::endl;
        }
    }

    f_ofs.close();
}


template<typename T>
void
test_mass_matrices(const mesh<T>& msh, size_t degree)
{
    auto rhs_fun = [](const typename mesh<T>::point_type& pt) -> T {
        return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
    };

    std::ofstream c_ofs("cell_mass_check.dat");

    cell_basis<T, T>::print_structure(degree);

    for (auto& cl : msh.cells)
    {
        Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, degree);
        Matrix<T, Dynamic, 1> rhs = make_rhs(msh, cl, degree, rhs_fun);
        Matrix<T, Dynamic, 1> sol = mass.llt().solve(rhs);

        cell_basis<T,T> cb(msh, cl, degree);

        auto tps = make_test_points(msh, cl);
        for (auto& tp : tps)
        {
            auto phi = cb.eval_basis(tp);
            auto val = sol.dot(phi);
            c_ofs << tp.x() << " " << tp.y() << " " << val << std::endl;
        }

    }

    c_ofs.close();


    std::ofstream f_ofs("face_mass_check.dat");

    for (auto& fc : msh.faces)
    {
        Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, degree);
        Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, degree, rhs_fun);
        Matrix<T, Dynamic, 1> sol = mass.llt().solve(rhs);

        face_basis<T,T> fb(msh, fc, degree);

        auto tps = make_test_points(msh, fc);
        for (auto& tp : tps)
        {
            auto phi = fb.eval_basis(tp);
            auto val = sol.dot(phi);
            f_ofs << tp.x() << " " << tp.y() << " " << val << std::endl;
        }

    }

    f_ofs.close();
}


template<typename T>
class assembler
{
    std::vector<size_t>                 compress_table;
    std::vector<size_t>                 expand_table;
    //std::vector<Matrix<T, Dynamic, 1>>  dirichlet_data;
    size_t                              degree;

    std::vector< Triplet<T> >           triplets;

    class assembly_index
    {
        size_t  idx;
        bool    assem;

    public:
        assembly_index(size_t i, bool as)
            : idx(i), assem(as)
        {}

        operator size_t() const
        {
            if (!assem)
                throw std::logic_error("Invalid assembly_index");

            return idx;
        }

        bool assemble() const
        {
            return assem;
        }

        friend std::ostream& operator<<(std::ostream& os, const assembly_index& as)
        {
            os << "(" << as.idx << "," << as.assem << ")";
            return os;
        }
    };

public:

    SparseMatrix<T>         LHS;
    Matrix<T, Dynamic, 1>   RHS;

    assembler(const mesh<T>& msh, size_t deg)
        : degree(deg)
    {
        auto is_dirichlet = [&](const typename mesh<T>::face_type& fc) -> bool {
            return fc.is_boundary && fc.bndtype == boundary::DIRICHLET;
        };

        auto num_all_faces = msh.faces.size();
        auto num_dirichlet_faces = std::count_if(msh.faces.begin(), msh.faces.end(), is_dirichlet);
        auto num_other_faces = num_all_faces - num_dirichlet_faces;

        compress_table.resize( num_all_faces );
        expand_table.resize( num_other_faces );
        //dirichlet_data.resize( num_dirichlet_faces );

        size_t compressed_offset = 0;
        for (size_t i = 0; i < num_all_faces; i++)
        {
            auto fc = msh.faces[i];
            if ( !is_dirichlet(fc) )
            {
                compress_table.at(i) = compressed_offset;
                expand_table.at(compressed_offset) = i;
                compressed_offset++;
            }
        }

        auto cbs = cell_basis<T,T>::size(degree);
        auto fbs = face_basis<T,T>::size(degree);
        auto system_size = cbs * msh.cells.size() + fbs * num_other_faces;

        LHS = SparseMatrix<T>( system_size, system_size );
        RHS = Matrix<T, Dynamic, 1>::Zero( system_size );
    }

    void dump_tables() const
    {
        std::cout << "Compress table: " << std::endl;
        for (size_t i = 0; i < compress_table.size(); i++)
            std::cout << i << " -> " << compress_table.at(i) << std::endl;
    }

    template<typename Function>
    void
    assemble(const mesh<T>& msh, const typename mesh<T>::cell_type& cl,
             const Matrix<T, Dynamic, Dynamic>& lhs, const Matrix<T, Dynamic, 1>& rhs,
             const Function& dirichlet_bf)
    {
        auto cbs = cell_basis<T,T>::size(degree);
        auto fbs = face_basis<T,T>::size(degree);

        std::vector<assembly_index> asm_map;
        asm_map.reserve(cbs + 4*fbs);

        auto cell_offset        = offset(msh, cl);
        auto cell_LHS_offset    = cell_offset * cbs;

        for (size_t i = 0; i < cbs; i++)
            asm_map.push_back( assembly_index(cell_LHS_offset+i, true) );

        Matrix<T, Dynamic, 1> dirichlet_data = Matrix<T, Dynamic, 1>::Zero(cbs + 4*fbs);

        auto fcs = faces(msh, cl);
        for (size_t face_i = 0; face_i < 4; face_i++)
        {
            auto fc = fcs[face_i];
            auto face_offset = offset(msh, fc);
            auto face_LHS_offset = cbs * msh.cells.size() + compress_table.at(face_offset)*fbs;

            bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET;

            for (size_t i = 0; i < fbs; i++)
                asm_map.push_back( assembly_index(face_LHS_offset+i, !dirichlet) );

            if (dirichlet)
            {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, degree);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, degree, dirichlet_bf);
                dirichlet_data.block(cbs+face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }
        }

        assert( asm_map.size() == lhs.rows() && asm_map.size() == lhs.cols() );

        for (size_t i = 0; i < lhs.rows(); i++)
        {
            if (!asm_map[i].assemble())
                continue;

            for (size_t j = 0; j < lhs.cols(); j++)
            {
                if ( asm_map[j].assemble() )
                    triplets.push_back( Triplet<T>(asm_map[i], asm_map[j], lhs(i,j)) );
                else
                    RHS(asm_map[i]) -= lhs(i,j)*dirichlet_data(j);
            }
        }

        RHS.block(cell_LHS_offset, 0, cbs, 1) += rhs.block(0, 0, cbs, 1);
    } // assemble()

    template<typename Function>
    Matrix<T, Dynamic, 1>
    take_local_data(const mesh<T>& msh, const typename mesh<T>::cell_type& cl,
    const Matrix<T, Dynamic, 1>& solution, const Function& dirichlet_bf)
    {
        auto cbs = cell_basis<T,T>::size(degree);
        auto fbs = face_basis<T,T>::size(degree);

        auto cell_offset        = offset(msh, cl);
        auto cell_SOL_offset    = cell_offset * cbs;

        Matrix<T, Dynamic, 1> ret = Matrix<T, Dynamic, 1>::Zero(cbs + 4*fbs);
        ret.block(0, 0, cbs, 1) = solution.block(cell_SOL_offset, 0, cbs, 1);

        auto fcs = faces(msh, cl);
        for (size_t face_i = 0; face_i < 4; face_i++)
        {
            auto fc = fcs[face_i];

            bool dirichlet = fc.is_boundary && fc.bndtype == boundary::DIRICHLET;

            if (dirichlet)
            {
                Matrix<T, Dynamic, Dynamic> mass = make_mass_matrix(msh, fc, degree);
                Matrix<T, Dynamic, 1> rhs = make_rhs(msh, fc, degree, dirichlet_bf);
                ret.block(cbs+face_i*fbs, 0, fbs, 1) = mass.llt().solve(rhs);
            }
            else
            {
                auto face_offset = offset(msh, fc);
                auto face_SOL_offset = cbs * msh.cells.size() + compress_table.at(face_offset)*fbs;
                ret.block(cbs+face_i*fbs, 0, fbs, 1) = solution.block(face_SOL_offset, 0, fbs, 1);
            }
        }

        return ret;
    }

    void finalize(void)
    {
        LHS.setFromTriplets( triplets.begin(), triplets.end() );
        triplets.clear();
    }
};


template<typename T>
auto make_assembler(const mesh<T>& msh, size_t deg)
{
    return assembler<T>(msh, deg);
}




#if 0
auto rhs_fun = [](const typename mesh<RealType>::point_type& pt) -> RealType {
    return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
};

std::ofstream ofs("gr_check.dat");

auto assembler = make_assembler(msh, degree);

for (auto& cl : msh.cells)
{
    Matrix<RealType, Dynamic, Dynamic> gr = make_hho_laplacian(msh, cl, degree);
    Matrix<RealType, Dynamic, Dynamic> stab = make_hho_naive_stabilization(msh, cl, degree);
    Matrix<RealType, Dynamic, 1> f = project_function(msh, cl, degree, rhs_fun);

    Matrix<RealType, Dynamic, 1> g = gr*f;

    cell_basis<RealType,RealType> rb(msh, cl, degree+1);

    auto tps = make_test_points(msh, cl);
    for (auto& tp : tps)
    {
        Matrix<RealType, Dynamic, 1> phi = rb.eval_basis(tp);
        auto val = g.dot( phi.tail(rb.size()-1) ) + f(0);
        ofs << tp.x() << " " << tp.y() << " " << val << std::endl;
    }


    assembler.assemble(msh, cl, stab, f, rhs_fun);

}

assembler.finalize();

ofs.close();
#endif

template<typename T>
void dump_mesh(const mesh<T>& msh)
{
    std::ofstream ofs("mesh_dump.m");
    size_t i = 0;
    for (auto& fc : msh.faces)
    {
        auto pts = points(msh, fc);
        if (fc.is_boundary)
            ofs << "line([" << pts[0].x() << ", " << pts[1].x() << "], [" << pts[0].y() << ", " << pts[1].y() << "], 'Color', 'r');" << std::endl;
        else
            ofs << "line([" << pts[0].x() << ", " << pts[1].x() << "], [" << pts[0].y() << ", " << pts[1].y() << "]);" << std::endl;

        auto bar = barycenter(msh, fc);
        ofs << "text(" << bar.x() << ", " << bar.y() << ", '" << i << "');" << std::endl;

        i++;
    }
    ofs.close();
}


template<typename T>
void add_interface_edges(const mesh<T>& msh)
{
    for (auto& cl : msh.cells)
    {
        auto fcs = faces(msh, cl);
        for (auto& fc : fcs)
        {
            auto pts = points(msh, fc);
            auto len = (pts[1] - pts[0]).to_vector().norm();
        }
    }
}


int main(int argc, char **argv)
{
    using RealType = double;
    size_t degree = 0;

    mesh_init_params<RealType> mip;
    mip.Nx = 5;
    mip.Ny = 5;

    int ch;
    while ( (ch = getopt(argc, argv, "k:M:N:")) != -1 )
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

            case '?':
            default:
                std::cout << "wrong arguments" << std::endl;
                exit(1);
        }
    }

    argc -= optind;
    argv += optind;

    mesh<RealType> msh(mip);

    dump_mesh(msh);

    silo_database silo;
    silo.create("cuthho_square.silo");
    silo.add_mesh(msh, "mesh");



    auto level_set_function = [](const typename mesh<RealType>::point_type pt) -> RealType {
        auto x = pt.x();
        auto y = pt.y();
        auto alpha = 0.5;
        auto beta = 0.5;

        return (x-alpha)*(x-alpha) + (y-beta)*(y-beta) - 0.15;
    };

    detect_cut_cells(msh, level_set_function);

    std::vector<RealType> cut_cell_markers;
    for (auto& cl : msh.cells)
    {
        if (cl.cut_by_interface)
            cut_cell_markers.push_back(1.0);
        else
            cut_cell_markers.push_back(0.0);
    }
    silo.add_variable("mesh", "cut_cells", cut_cell_markers.data(), cut_cell_markers.size(), zonal_variable_t);


    std::vector<RealType> level_set_vals;
    for (auto& pt : msh.points)
        level_set_vals.push_back( level_set_function(pt) );
    silo.add_variable("mesh", "level_set", level_set_vals.data(), level_set_vals.size(), nodal_variable_t);



    auto rhs_fun = [](const typename mesh<RealType>::point_type& pt) -> RealType {
        return 2.0 * M_PI * M_PI * std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
        //return pt.x() * pt.x() * pt.x();
    };

    auto sol_fun = [](const typename mesh<RealType>::point_type& pt) -> RealType {
        return std::sin(M_PI*pt.x()) * std::sin(M_PI*pt.y());
        //return pt.x();
    };

    for (auto& fc : msh.faces)
    {
        if (fc.is_boundary)
            fc.bndtype = boundary::DIRICHLET;
    }

    auto assembler = make_assembler(msh, degree);
    for (auto& cl : msh.cells)
    {
        auto gr = make_operator(msh, cl, degree);
        Matrix<RealType, Dynamic, Dynamic> stab = make_stabilization(msh, cl, degree);
        Matrix<RealType, Dynamic, Dynamic> lc = gr.second + stab;
        //Matrix<RealType, Dynamic, 1> ft = project_function(msh, cl, degree, rhs_fun);
        Matrix<RealType, Dynamic, 1> f = make_rhs(msh, cl, degree, rhs_fun);
        assembler.assemble(msh, cl, lc, f, sol_fun);
    }

    assembler.finalize();

#if 0
    std::ofstream mat_ofs("matrix.dat");
    for (int k=0; k<assembler.LHS.outerSize(); ++k)
      for (SparseMatrix<RealType>::InnerIterator it(assembler.LHS,k); it; ++it)
      {

        mat_ofs << it.row() << " " << it.col() << " " << it.value() << std::endl;
      }
    mat_ofs.close();
#endif


    std::cout << "System unknowns: " << assembler.LHS.rows() << std::endl;
    std::cout << "Cells: " << msh.cells.size() << std::endl;
    std::cout << "Faces: " << msh.faces.size() << std::endl;

//#if 0
    SparseLU<SparseMatrix<RealType>>  solver;

    solver.analyzePattern(assembler.LHS);
    solver.factorize(assembler.LHS);
    Matrix<RealType, Dynamic, 1> sol = solver.solve(assembler.RHS);
//#endif
#if 0
    Matrix<RealType, Dynamic, 1> sol = Matrix<RealType, Dynamic, 1>::Zero(assembler.RHS.rows());
    conjugated_gradient(assembler.LHS, assembler.RHS, sol);
#endif

    std::ofstream ofs("solution.dat");
    RealType error_int = 0.0;
    RealType error_mm = 0.0;

    std::vector<RealType>   error_int_percell;
    std::vector<RealType>   error_mm_percell;
    std::vector<RealType>   solution_val;
    std::vector<RealType>   stab_val;

    size_t cell_i = 0;
    for (auto& cl : msh.cells)
    {
        cell_basis<RealType,RealType> cb(msh, cl, degree);
        auto cbs = cb.size();

        Matrix<RealType, Dynamic, 1> cdofs = sol.block(cell_i*cbs, 0, cbs, 1);

        auto tps = make_test_points(msh, cl);
        for (auto& tp : tps)
        {
            Matrix<RealType, Dynamic, 1> phi = cb.eval_basis(tp);
            auto val = cdofs.dot( phi );
            ofs << tp.x() << " " << tp.y() << " " << val << std::endl;
        }

        RealType e_int_c = 0.0;
        auto qps = integrate(msh, cl, 2*degree);
        for (auto& qp : qps)
        {
            Matrix<RealType, Dynamic, 1> phi = cb.eval_basis(qp.first);
            auto val = cdofs.dot( phi );
            auto real_val = sol_fun( qp.first );
            e_int_c += qp.second*(real_val - val)*(real_val - val);
        }
        error_int_percell.push_back(e_int_c);
        error_int += e_int_c;

        Matrix<RealType, Dynamic, Dynamic> mass = make_mass_matrix(msh, cl, degree);
        Matrix<RealType, Dynamic, 1> rhs = make_rhs(msh, cl, degree, sol_fun);
        Matrix<RealType, Dynamic, 1> real_dofs = mass.llt().solve(rhs);
        Matrix<RealType, Dynamic, 1> diff = real_dofs - cdofs;
        RealType e_mm_c = diff.dot(mass*diff);
        error_mm_percell.push_back(e_mm_c);
        error_mm += e_mm_c;

        {
            auto bar = barycenter(msh, cl);
            Matrix<RealType, Dynamic, 1> phi = cb.eval_basis(bar);
            auto val = cdofs.dot( phi );
            solution_val.push_back(val);
        }

        Matrix<RealType, Dynamic, Dynamic> stab = make_stabilization(msh, cl, degree);
        Matrix<RealType, Dynamic, 1> c_alldofs = assembler.take_local_data(msh, cl, sol, sol_fun);
        auto sv = c_alldofs.dot(stab*c_alldofs);
        stab_val.push_back(sv);

        cell_i++;
    }
    ofs.close();

    silo.add_variable("mesh", "err_int_c", error_int_percell.data(), error_int_percell.size(), zonal_variable_t);
    silo.add_variable("mesh", "err_mm_c", error_mm_percell.data(), error_mm_percell.size(), zonal_variable_t);
    silo.add_variable("mesh", "u", solution_val.data(), solution_val.size(), zonal_variable_t);
    silo.add_variable("mesh", "stab_val", stab_val.data(), stab_val.size(), zonal_variable_t);

    std::cout << "L2-error (int): " << sqrt(error_int) << std::endl;
    std::cout << "L2-error (mm) : " << sqrt(error_mm) << std::endl;

    return 0;
}
