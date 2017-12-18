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

#include <Eigen/Dense>

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

template<typename Mesh, typename VT>
class cell_basis
{
    typedef typename Mesh::coordinate_type  coordinate_type;
    typedef typename Mesh::point_type       point_type;

    point_type          cell_bar;
    coordinate_type     cell_h;
    size_t              basis_degree, basis_size;

public:
    cell_basis(const Mesh& msh, const typename Mesh::cell_type& cl, size_t degree)
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

template<typename Mesh, typename VT>
class face_basis
{
    typedef typename Mesh::coordinate_type  coordinate_type;
    typedef typename Mesh::point_type       point_type;

    point_type          face_bar;
    point_type          base;
    coordinate_type     face_h;
    size_t              basis_degree, basis_size;

public:
    face_basis(const Mesh& msh, const typename Mesh::face_type& fc, size_t degree)
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
