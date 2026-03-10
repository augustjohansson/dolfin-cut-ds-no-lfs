// Copyright (C) 2014 August Johansson and Anders Logg
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

#include <cassert>
#include <cmath>
#include <stdexcept>
#include "SimplexQuadrature.h"
#include "SimplexQuadratureCompression.h"
#include "SimplexQuadratureTables.h"
#include "predicates.h"

using namespace simpex;

// ---------------------------------------------------------------------------
// File-scope helpers
// ---------------------------------------------------------------------------
static void write_point(double* out, std::size_t gdim, const Point& p)
{
  out[0] = p.x();
  if (gdim > 1) out[1] = p.y();
  if (gdim > 2) out[2] = p.z();
}

static double segment_length(const Point& a, const Point& b, std::size_t gdim)
{
  const Point e = b - a;
  if (gdim == 1) return std::abs(e.x());
  if (gdim == 2) return std::sqrt(e.x()*e.x() + e.y()*e.y());
  return e.norm();
}

static double triangle_area(const Point& a, const Point& b, const Point& c,
                             std::size_t gdim)
{
  if (gdim == 2) return 0.5 * std::abs(orient2d(a, b, c));
  const Point u = b - a;
  const Point v = c - a;
  return 0.5 * u.cross(v).norm();
}

static double tetra_volume(const Point& a, const Point& b,
                            const Point& c, const Point& d)
{
  return std::abs(orient3d(a, b, c, d)) / 6.0;
}

//-----------------------------------------------------------------------------
SimplexQuadrature::SimplexQuadrature(std::size_t tdim, std::size_t order)
{
  // Create and store quadrature rule for reference simplex
  switch (tdim)
    {
    case 1:
      setup_qr_reference_interval(order);
      break;
    case 2:
      setup_qr_reference_triangle(order);
      break;
    case 3:
      setup_qr_reference_tetrahedron(order);
      break;
    default:
      throw std::runtime_error("Only implemented for topological dimension 1, 2, 3");
    }

}
//-----------------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
SimplexQuadrature::compute_quadrature_rule(const std::vector<Point>& coordinates,
					   std::size_t gdim) const
{
  std::size_t tdim = coordinates.size() - 1;

  switch (tdim)
    {
    case 0:
      // FIXME: should we return empty qr or should we have detected this earlier?
      break;
    case 1:
      return compute_quadrature_rule_interval(coordinates, gdim);
      break;
    case 2:
      return compute_quadrature_rule_triangle(coordinates, gdim);
      break;
    case 3:
      return compute_quadrature_rule_tetrahedron(coordinates, gdim);
      break;
    default:
      throw std::runtime_error("Only implemented for topological dimension 1, 2, 3");
    };

  std::pair<std::vector<double>, std::vector<double>> quadrature_rule;
  return quadrature_rule;
}
//-----------------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
SimplexQuadrature::compute_quadrature_rule_interval(const std::vector<Point>& X,
                                                    std::size_t gdim) const
{
  if (X.size() != 2) throw std::runtime_error("interval: need 2 vertices");
  const std::size_t n = _p.size();

  std::vector<double> pts(n * gdim);
  std::vector<double> w(n);

  const Point& x0 = X[0];
  const Point  e  = X[1] - X[0];

  const double L = segment_length(X[0], X[1], gdim);
  const double ref = 1.0;                              // |[0,1]|
  const double scale = L / ref;

  for (std::size_t i = 0; i < n; ++i)
    {
      const double r = _p[i][0];
      const Point p = x0 + e * r;
      write_point(&pts[i*gdim], gdim, p);
      w[i] = _w[i] * scale;
    }

  return {std::move(pts), std::move(w)};
}

//-----------------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
SimplexQuadrature::compute_quadrature_rule_triangle(const std::vector<Point>& X,
                                                    std::size_t gdim) const
{
  if (X.size() != 3) throw std::runtime_error("triangle: need 3 vertices");
  const std::size_t n = _p.size();

  std::vector<double> pts(n * gdim);
  std::vector<double> w(n);

  const Point& x0 = X[0];
  const Point  e1 = X[1] - X[0];
  const Point  e2 = X[2] - X[0];

  const double A = triangle_area(X[0], X[1], X[2], gdim);
  const double ref = 0.5; // reference triangle area
  const double scale = A / ref;

  for (std::size_t i = 0; i < n; ++i)
    {
      const double r = _p[i][0];
      const double s = _p[i][1];
      const Point p = x0 + e1 * r + e2 * s;
      write_point(&pts[i*gdim], gdim, p);
      w[i] = _w[i] * scale;
    }

  return {std::move(pts), std::move(w)};
}

//-----------------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
SimplexQuadrature::compute_quadrature_rule_tetrahedron(const std::vector<Point>& X,
                                                       std::size_t gdim) const
{
  if (X.size() != 4) throw std::runtime_error("tetrahedron: need 4 vertices");
  if (gdim != 3) throw std::runtime_error("tetrahedron: requires gdim=3");
  const std::size_t n = _p.size();

  std::vector<double> pts(n * 3);
  std::vector<double> w(n);

  const Point& x0 = X[0];
  const Point  e1 = X[1] - X[0];
  const Point  e2 = X[2] - X[0];
  const Point  e3 = X[3] - X[0];

  const double V = tetra_volume(X[0], X[1], X[2], X[3]);
  const double ref = 1.0 / 6.0;
  const double scale = V / ref;

  for (std::size_t i = 0; i < n; ++i)
    {
      const double r = _p[i][0];
      const double s = _p[i][1];
      const double t = _p[i][2];
      const Point p = x0 + e1 * r + e2 * s + e3 * t;
      pts[i*3 + 0] = p.x();
      pts[i*3 + 1] = p.y();
      pts[i*3 + 2] = p.z();
      w[i] = _w[i] * scale;
    }

  return {std::move(pts), std::move(w)};
}

//-----------------------------------------------------------------------------
void SimplexQuadrature::setup_qr_reference_interval(std::size_t order)
{
  const auto rr = tables::interval_rule_legendre(order);
  const std::size_t n = rr.w.size();
  _p.assign(n, std::vector<double>(1));
  _w = rr.w;
  for (std::size_t i = 0; i < n; ++i)
    _p[i][0] = rr.pts[3*i + 0]; // r
}

//-----------------------------------------------------------------------------
void SimplexQuadrature::setup_qr_reference_triangle(std::size_t order)
{
  if (order > 20)
    throw std::runtime_error("Triangle quadrature supported up to order 20.");
  const auto rr = tables::triangle_rule_dunavant(order);
  const std::size_t n = rr.w.size();
  _p.assign(n, std::vector<double>(2));
  _w = rr.w;
  for (std::size_t i = 0; i < n; ++i)
  {
    _p[i][0] = rr.pts[3*i + 0]; // r
    _p[i][1] = rr.pts[3*i + 1]; // s
  }
}

//-----------------------------------------------------------------------------
void SimplexQuadrature::setup_qr_reference_tetrahedron(std::size_t order)
{
  if (order > 6)
    throw std::runtime_error("Tetrahedron quadrature supported up to order 6.");
  const auto rr = tables::tetra_rule_keast(order);
  const std::size_t n = rr.w.size();
  _p.assign(n, std::vector<double>(3));
  _w = rr.w;
  for (std::size_t i = 0; i < n; ++i)
  {
    _p[i][0] = rr.pts[3*i + 0]; // r
    _p[i][1] = rr.pts[3*i + 1]; // s
    _p[i][2] = rr.pts[3*i + 2]; // t
  }
}
//-----------------------------------------------------------------------------
