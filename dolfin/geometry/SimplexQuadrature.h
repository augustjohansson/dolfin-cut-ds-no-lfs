// Copyright (C) 2014 Anders Logg
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

#pragma once

#include <cstddef>
#include <vector>
#include "Point.h"
#include "SimplexQuadratureCompression.h"

namespace simpex
{

/// Quadrature rules for simplices.
class SimplexQuadrature
{
public:

  /// Construct a quadrature rule for the reference simplex of the given
  /// topological dimension and polynomial order.
  ///
  /// \param tdim  Topological dimension (1 = interval, 2 = triangle, 3 = tetrahedron).
  /// \param order Order of the quadrature rule.
  SimplexQuadrature(std::size_t tdim, std::size_t order);

  /// Compute a quadrature rule for a physical simplex.
  ///
  /// \param coordinates Vertex coordinates of the simplex.
  /// \param gdim        Geometric (embedding) dimension.
  /// \returns A pair (points, weights), where points is a flat array of length
  ///          n * gdim and weights has length n.
  std::pair<std::vector<double>, std::vector<double>>
  compute_quadrature_rule(const std::vector<Point>& coordinates,
                          std::size_t gdim) const;

  /// Compute a quadrature rule for a physical interval.
  std::pair<std::vector<double>, std::vector<double>>
  compute_quadrature_rule_interval(const std::vector<Point>& coordinates,
                                   std::size_t gdim) const;

  /// Compute a quadrature rule for a physical triangle.
  std::pair<std::vector<double>, std::vector<double>>
  compute_quadrature_rule_triangle(const std::vector<Point>& coordinates,
                                   std::size_t gdim) const;

  /// Compute a quadrature rule for a physical tetrahedron.
  std::pair<std::vector<double>, std::vector<double>>
  compute_quadrature_rule_tetrahedron(const std::vector<Point>& coordinates,
                                      std::size_t gdim) const;

  /// Compress a quadrature rule.  Delegates to compress_quadrature_rule()
  /// in SimplexQuadratureCompression.h.
  ///
  /// \param[in,out] qr              Quadrature rule to compress.
  /// \param[in]     gdim            Geometric dimension.
  /// \param[in]     quadrature_order Polynomial order integrated exactly.
  /// \param[in]     method           Compression algorithm to use.
  /// \returns Indices of the retained points, or empty if no compression was made.
  static std::vector<std::size_t>
  compress(std::pair<std::vector<double>, std::vector<double>>& qr,
           std::size_t gdim,
           std::size_t quadrature_order,
           compression::CompressionMethod method
               = compression::CompressionMethod::legacy)
  {
    return compression::compress(qr, gdim, quadrature_order, method);
  }

private:

  void setup_qr_reference_interval(std::size_t order);
  void setup_qr_reference_triangle(std::size_t order);
  void setup_qr_reference_tetrahedron(std::size_t order);

  // Reference quadrature points and weights
  std::vector<std::vector<double>> _p;
  std::vector<double>              _w;

};

} // namespace simpex

// ============================================================================
// dolfin::SimplexQuadrature
//
// Extends simpex::SimplexQuadrature with a convenience overload that accepts
// a dolfin::Cell directly (used in MultiMesh.cpp).
// ============================================================================

namespace dolfin
{
  // Forward declaration to avoid circular includes.
  class Cell;

  /// SimplexQuadrature with dolfin Cell support.
  class SimplexQuadrature : public simpex::SimplexQuadrature
  {
  public:
    using simpex::SimplexQuadrature::SimplexQuadrature;

    /// Compute quadrature rule for a dolfin Cell.
    std::pair<std::vector<double>, std::vector<double>>
    compute_quadrature_rule(const Cell& cell) const;
  };
}

