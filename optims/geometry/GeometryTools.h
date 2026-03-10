// Copyright (C) 2017 Anders Logg
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

#ifndef __GEOMETRY_TOOLS_H
#define __GEOMETRY_TOOLS_H

#include "Point.h"
#include "predicates.h"

namespace simpex
{

/// Utility functions for computational geometry.
class GeometryTools
{
public:

  /// Compute the numerically stable cross product (a - c) x (b - c).
  /// See Shewchuk, Lecture Notes on Geometric Robustness.
  static inline Point cross_product(const Point& a, const Point& b, const Point& c)
  {
    double ayz[2] = {a.y(), a.z()};
    double byz[2] = {b.y(), b.z()};
    double cyz[2] = {c.y(), c.z()};
    double azx[2] = {a.z(), a.x()};
    double bzx[2] = {b.z(), b.x()};
    double czx[2] = {c.z(), c.x()};
    double axy[2] = {a.x(), a.y()};
    double bxy[2] = {b.x(), b.y()};
    double cxy[2] = {c.x(), c.y()};
    return Point(_orient2d(ayz, byz, cyz),
                 _orient2d(azx, bzx, czx),
                 _orient2d(axy, bxy, cxy));
  }

};

} // namespace simpex

#endif
