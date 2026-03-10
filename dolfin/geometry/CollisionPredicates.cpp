// Copyright (C) 2014-2017 Anders Logg, August Johansson and Benjamin Kehlet
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
//
// First added:  2014-02-03
// Last changed: 2017-10-09

#include <cassert>
#include <stdexcept>
#include "predicates.h"
#include "Point.h"
#include "CollisionPredicates.h"
#include "GeometryTools.h"

#include "CGALExactArithmetic.h"

using namespace simpex;

//-----------------------------------------------------------------------------
// Low-level collision detection predicates
//-----------------------------------------------------------------------------
bool CollisionPredicates::collides_segment_point(const Point& p0,
                                                 const Point& p1,
                                                 const Point& point,
                                                 std::size_t gdim)
{
  switch (gdim)
    {
    case 1:
      return collides_segment_point_1d(p0[0], p1[0], point[0]);
    case 2:
      return collides_segment_point_2d(p0, p1, point);
    case 3:
      return collides_segment_point_3d(p0, p1, point);
    default:
      throw std::runtime_error("Unknown dimension (only implemented for dimension 2 and 3");
    }
  return false;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::collides_segment_segment(const Point& p0,
                                                   const Point& p1,
                                                   const Point& q0,
                                                   const Point& q1,
                                                   std::size_t gdim)
{
  switch (gdim)
    {
    case 1:
      return collides_segment_segment_1d(p0[0], p1[0], q0[0], q1[0]);
    case 2:
      return collides_segment_segment_2d(p0, p1, q0, q1);
    case 3:
      return collides_segment_segment_3d(p0, p1, q0, q1);
    default:
      throw std::runtime_error("Unknown dimension (Implemented for dimension 1, 2 and 3)");
    }
  return false;
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_triangle_point(const Point& p0,
                                                  const Point& p1,
                                                  const Point& p2,
                                                  const Point& point,
                                                  std::size_t gdim)
{
  switch (gdim)
    {
    case 2:
      return collides_triangle_point_2d(p0, p1, p2, point);
    case 3:
      return collides_triangle_point_3d(p0, p1, p2, point);
    default:
      throw std::runtime_error("Implemented only for dimension 2 and 3.");
    }
  return false;
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_triangle_segment(const Point& p0,
                                                    const Point& p1,
                                                    const Point& p2,
                                                    const Point& q0,
                                                    const Point& q1,
                                                    std::size_t gdim)
{
  switch (gdim)
    {
    case 2:
      return collides_triangle_segment_2d(p0, p1, p2, q0, q1);
    case 3:
      return collides_triangle_segment_3d(p0, p1, p2, q0, q1);
    default:
      throw std::runtime_error("Implmented only for dimension 2 and 3.");
    }
  return false;
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_triangle_triangle(const Point& p0,
                                                     const Point& p1,
                                                     const Point& p2,
                                                     const Point& q0,
                                                     const Point& q1,
                                                     const Point& q2,
                                                     std::size_t gdim)
{
  switch (gdim)
    {
    case 2:
      return collides_triangle_triangle_2d(p0, p1, p2, q0, q1, q2);
    case 3:
      return collides_triangle_triangle_3d(p0, p1, p2, q0, q1, q2);
    default:
      throw std::runtime_error("Implmented only for dimension 2 and 3.");
    }
  return false;
}



//--- Low-level collision detection predicates ---
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_segment_point_1d(double p0,
						    double p1,
						    double point)
{
  // FIXME: Skip CGAL for now
  return _collides_segment_point_1d(p0, p1, point);
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_segment_point_2d(const Point& p0,
						    const Point& p1,
						    const Point& point)
{
  return CHECK_CGAL(_collides_segment_point_2d(p0, p1, point),
		    cgal_collides_segment_point_2d(p0, p1, point));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_segment_point_3d(const Point& p0,
						    const Point& p1,
						    const Point& point)
{
  return CHECK_CGAL(_collides_segment_point_3d(p0, p1, point),
		    cgal_collides_segment_point_3d(p0, p1, point));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_segment_segment_1d(double p0,
						      double p1,
						      double q0,
						      double q1)
{
  return _collides_segment_segment_1d(p0, p1, q0, q1);
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_segment_segment_2d(const Point& p0,
						      const Point& p1,
						      const Point& q0,
						      const Point& q1)
{
  return CHECK_CGAL(_collides_segment_segment_2d(p0, p1, q0, q1),
		    cgal_collides_segment_segment_2d(p0, p1, q0, q1));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_segment_segment_3d(const Point& p0,
						      const Point& p1,
						      const Point& q0,
						      const Point& q1)
{
  return CHECK_CGAL(_collides_segment_segment_3d(p0, p1, q0, q1),
		    cgal_collides_segment_segment_3d(p0, p1, q0, q1));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_triangle_point_2d(const Point& p0,
						     const Point& p1,
						     const Point& p2,
						     const Point& point)
{
  return CHECK_CGAL(_collides_triangle_point_2d(p0, p1, p2, point),
		    cgal_collides_triangle_point_2d(p0, p1, p2, point));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_triangle_point_3d(const Point& p0,
						     const Point& p1,
						     const Point& p2,
						     const Point& point)
{
  return CHECK_CGAL(_collides_triangle_point_3d(p0, p1, p2, point),
		    cgal_collides_triangle_point_3d(p0, p1, p2, point));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_triangle_segment_2d(const Point& p0,
						       const Point& p1,
						       const Point& p2,
						       const Point& q0,
						       const Point& q1)
{
  return CHECK_CGAL(_collides_triangle_segment_2d(p0, p1, p2, q0, q1),
		    cgal_collides_triangle_segment_2d(p0, p1, p2, q0, q1));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_triangle_segment_3d(const Point& p0,
						       const Point& p1,
						       const Point& p2,
						       const Point& q0,
						       const Point& q1)
{
  return CHECK_CGAL(_collides_triangle_segment_3d(p0, p1, p2, q0, q1),
		    cgal_collides_triangle_segment_3d(p0, p1, p2, q0, q1));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_triangle_triangle_2d(const Point& p0,
							const Point& p1,
							const Point& p2,
							const Point& q0,
							const Point& q1,
							const Point& q2)
{
  return CHECK_CGAL(_collides_triangle_triangle_2d(p0, p1, p2, q0, q1, q2),
		    cgal_collides_triangle_triangle_2d(p0, p1, p2, q0, q1, q2));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_triangle_triangle_3d(const Point& p0,
							const Point& p1,
							const Point& p2,
							const Point& q0,
							const Point& q1,
							const Point& q2)
{
  return CHECK_CGAL(_collides_triangle_triangle_3d(p0, p1, p2, q0, q1, q2),
		    cgal_collides_triangle_triangle_3d(p0, p1, p2, q0, q1, q2));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_tetrahedron_point_3d(const Point& p0,
							const Point& p1,
							const Point& p2,
							const Point& p3,
							const Point& point)
{
  return CHECK_CGAL(_collides_tetrahedron_point_3d(p0, p1, p2, p3, point),
		    cgal_collides_tetrahedron_point_3d(p0, p1, p2, p3, point));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_tetrahedron_segment_3d(const Point& p0,
							  const Point& p1,
							  const Point& p2,
							  const Point& p3,
							  const Point& q0,
							  const Point& q1)
{
  return CHECK_CGAL(_collides_tetrahedron_segment_3d(p0, p1, p2, p3, q0, q1),
		    cgal_collides_tetrahedron_segment_3d(p0, p1, p2, p3, q0, q1));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_tetrahedron_triangle_3d(const Point& p0,
							   const Point& p1,
							   const Point& p2,
							   const Point& p3,
							   const Point& q0,
							   const Point& q1,
							   const Point& q2)
{
  return CHECK_CGAL(_collides_tetrahedron_triangle_3d(p0, p1, p2, p3, q0, q1, q2),
		    cgal_collides_tetrahedron_triangle_3d(p0, p1, p2, p3, q0, q1, q2));
}
//------------------------------------------------------------------------------
bool CollisionPredicates::collides_tetrahedron_tetrahedron_3d(const Point& p0,
							      const Point& p1,
							      const Point& p2,
							      const Point& p3,
							      const Point& q0,
							      const Point& q1,
							      const Point& q2,
							      const Point& q3)
{
  return CHECK_CGAL(_collides_tetrahedron_tetrahedron_3d(p0, p1, p2, p3, q0, q1, q2, q3),
		    cgal_collides_tetrahedron_tetrahedron_3d(p0, p1, p2, p3, q0, q1, q2, q3));
}
//-----------------------------------------------------------------------------
// Implementation of private members
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_segment_point_1d(double p0,
                                                     double p1,
                                                     double point)
{
  if (p0 > p1)
    std::swap(p0, p1);
  return p0 <= point and point <= p1;
}
//-----------------------------------------------------------------------------
// Returns true iff p lies within the bounding box of segment [a,b] (2D check).
static inline bool on_segment_bbox_2d(const Point& a, const Point& b,
                                       const Point& p) noexcept
{
  return (std::min(a.x(), b.x()) <= p.x() && p.x() <= std::max(a.x(), b.x()) &&
          std::min(a.y(), b.y()) <= p.y() && p.y() <= std::max(a.y(), b.y()));
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_segment_point_2d(const Point& p0,
                                                     const Point& p1,
                                                     const Point& point)
{
  if (orient2d(p0, p1, point) != 0.0) return false;
  return on_segment_bbox_2d(p0, p1, point);
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_segment_point_3d(const Point& p0,
                                                     const Point& p1,
                                                     const Point& point)
{
  if (point == p0 || point == p1)
    return true;

  // Project to 2D planes in turn; point is collinear iff all three are zero.
  // xy-plane
  const double a0[2] = {p0.x(),    p0.y()   };
  const double a1[2] = {p1.x(),    p1.y()   };
  const double a2[2] = {point.x(), point.y()};
  if (_orient2d(a0, a1, a2) != 0.0) return false;

  // xz-plane
  const double b0[2] = {p0.x(),    p0.z()   };
  const double b1[2] = {p1.x(),    p1.z()   };
  const double b2[2] = {point.x(), point.z()};
  if (_orient2d(b0, b1, b2) != 0.0) return false;

  // yz-plane
  const double c0[2] = {p0.y(),    p0.z()   };
  const double c1[2] = {p1.y(),    p1.z()   };
  const double c2[2] = {point.y(), point.z()};
  if (_orient2d(c0, c1, c2) != 0.0) return false;

  // Point is collinear: check it lies between p0 and p1.
  const double length = (p0 - p1).squared_norm();
  return (point - p0).squared_norm() <= length &&
         (point - p1).squared_norm() <= length;
}
//------------------------------------------------------------------------------
bool CollisionPredicates::_collides_segment_segment_1d(double p0,
                                                       double p1,
                                                       double q0,
                                                       double q1)
{
  // Get range
  const double a0 = std::min(p0, p1);
  const double b0 = std::max(p0, p1);
  const double a1 = std::min(q0, q1);
  const double b1 = std::max(q0, q1);

  // Check for collision
  const double dx = std::min(b0 - a0, b1 - a1);
  return b1 >= a0 - dx && a1 <= b0 + dx;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_segment_point_2d_with_hint(const Point& p0,
                                                                const Point& p1,
                                                                const Point& point,
                                                                double op01_pt)
{
  // op01_pt = orient2d(p0, p1, point) is pre-computed.
  if (op01_pt != 0.0) return false;
  return on_segment_bbox_2d(p0, p1, point);
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_segment_segment_2d_with_hint(const Point& p0,
                                                                  const Point& p1,
                                                                  const Point& q0,
                                                                  const Point& q1,
                                                                  double pq0,
                                                                  double pq1,
                                                                  double qp0,
                                                                  double qp1)
{
  // pq0 = orient2d(p0,p1,q0)  pq1 = orient2d(p0,p1,q1)
  // qp0 = orient2d(q0,q1,p0)  qp1 = orient2d(q0,q1,p1)

  // Proper (non-degenerate) crossing: each segment straddles the other.
  if (((pq0 > 0.0) != (pq1 > 0.0)) && ((qp0 > 0.0) != (qp1 > 0.0)))
    return true;

  // Collinear / touching: at least one endpoint lies on the other segment.
  if (pq0 == 0.0 && on_segment_bbox_2d(p0, p1, q0)) return true;
  if (pq1 == 0.0 && on_segment_bbox_2d(p0, p1, q1)) return true;
  if (qp0 == 0.0 && on_segment_bbox_2d(q0, q1, p0)) return true;
  if (qp1 == 0.0 && on_segment_bbox_2d(q0, q1, p1)) return true;

  return false;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_segment_segment_2d(const Point& p0,
                                                       const Point& p1,
                                                       const Point& q0,
                                                       const Point& q1)
{
  const double o1 = orient2d(p0, p1, q0);
  const double o2 = orient2d(p0, p1, q1);
  const double o3 = orient2d(q0, q1, p0);
  const double o4 = orient2d(q0, q1, p1);

  // Proper intersection: each segment straddles the other's line.
  if (((o1 > 0.0) != (o2 > 0.0)) && ((o3 > 0.0) != (o4 > 0.0)))
    return true;

  // Collinear / touching: check bbox containment.
  if (o1 == 0.0 && on_segment_bbox_2d(p0, p1, q0)) return true;
  if (o2 == 0.0 && on_segment_bbox_2d(p0, p1, q1)) return true;
  if (o3 == 0.0 && on_segment_bbox_2d(q0, q1, p0)) return true;
  if (o4 == 0.0 && on_segment_bbox_2d(q0, q1, p1)) return true;

  return false;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_segment_segment_3d(const Point& p0,
                                                       const Point& p1,
                                                       const Point& q0,
                                                       const Point& q1)
{
  // Vertex collisions
  if (p0 == q0 || p0 == q1 || p1 == q0 || p1 == q1)
    return true;

  if (collides_segment_point_3d(p0, p1, q0) or
      collides_segment_point_3d(p0, p1, q1) or
      collides_segment_point_3d(q0, q1, p0) or
      collides_segment_point_3d(q0, q1, p1))
    {
      return true;
    }

  // Determinant must be zero
  const double det = orient3d(p0, p1, q0, q1);

  if (det < 0.0 or det > 0.0)
    return false;

  // Now we know that the segments are in the same plane. This means
  // that they can be parallel, or even collinear.

  // Check for collinearity
  const Point u = GeometryTools::cross_product(p0, p1, q0);
  if (u[0] == 0.0 and u[1] == 0.0 and u[2] == 0.0)
    {
      const Point v = GeometryTools::cross_product(p0, p1, q1);
      if (v[0] == 0.0 and v[1] == 0.0 and v[2] == 0.0)
	{
	  // Now we know that the segments are collinear
	  if ((p0-q0).squared_norm() <= (q1-q0).squared_norm() and
	      (p0-q1).squared_norm() <= (q0-q1).squared_norm())
	    return true;

	  if ((p1-q0).squared_norm() <= (q1-q0).squared_norm() and
	      (p1-q1).squared_norm() <= (q0-q1).squared_norm())
	    return true;

	  if ((q0-p0).squared_norm() <= (p1-p0).squared_norm() and
	      (q0-p1).squared_norm() <= (p0-p1).squared_norm())
	    return true;

	  if ((q1-p0).squared_norm() <= (p1-p0).squared_norm() and
	      (q1-p1).squared_norm() <= (p0-p1).squared_norm())
	    return true;
	}
    }

  // Segments are not collinear, but in the same plane
  // Try to reduce to 2d by elimination

  for (std::size_t d = 0; d < 3; ++d)
    {
      if (p0[d] == p1[d] and p0[d] == q0[d] and p0[d] == q1[d])
	{
	  const std::array<std::array<std::size_t, 2>, 3> dims = {{ {{1, 2}},
								    {{0, 2}},
								    {{0, 1}} }};
	  const Point p0_2d(p0[dims[d][0]], p0[dims[d][1]]);
	  const Point p1_2d(p1[dims[d][0]], p1[dims[d][1]]);
	  const Point q0_2d(q0[dims[d][0]], q0[dims[d][1]]);
	  const Point q1_2d(q1[dims[d][0]], q1[dims[d][1]]);

	  return collides_segment_segment_2d(p0_2d, p1_2d, q0_2d, q1_2d);
	}
    }

  return false;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_triangle_point_2d_with_hint(const Point& p0,
                                                                 const Point& p1,
                                                                 const Point& p2,
                                                                 const Point& point,
                                                                 double ref,
                                                                 double o01,
                                                                 double o12,
                                                                 double o20)
{
  // ref = orient2d(p0,p1,p2), o01 = orient2d(p0,p1,point),
  // o12 = orient2d(p1,p2,point), o20 = orient2d(p2,p0,point)
  if (ref > 0.0)
    return o12 >= 0.0 && o20 >= 0.0 && o01 >= 0.0;
  if (ref < 0.0)
    return o12 <= 0.0 && o20 <= 0.0 && o01 <= 0.0;
  // Degenerate (collinear) triangle: test edges
  return ((o01 == 0.0 &&
           collides_segment_point_1d(p0[0], p1[0], point[0]) &&
           collides_segment_point_1d(p0[1], p1[1], point[1])) ||
          (o12 == 0.0 &&
           collides_segment_point_1d(p1[0], p2[0], point[0]) &&
           collides_segment_point_1d(p1[1], p2[1], point[1])) ||
          (o20 == 0.0 &&
           collides_segment_point_1d(p2[0], p0[0], point[0]) &&
           collides_segment_point_1d(p2[1], p0[1], point[1])));
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_triangle_point_2d(const Point& p0,
                                                      const Point& p1,
                                                      const Point& p2,
                                                      const Point& point)
{
  const double ref = orient2d(p0, p1, p2);
  const double o01 = orient2d(p0, p1, point);
  const double o12 = orient2d(p1, p2, point);
  const double o20 = orient2d(p2, p0, point);
  return _collides_triangle_point_2d_with_hint(p0, p1, p2, point, ref, o01, o12, o20);
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_triangle_point_3d_in_plane(const Point& p0,
                                                                const Point& p1,
                                                                const Point& p2,
                                                                const Point& point)
{
  // Assumes orient3d(p0,p1,p2,point) == 0.0 -- point is in the triangle's plane.
  // Project to 2D by dropping the coordinate corresponding to the largest
  // component of the triangle normal.  This is exact with Shewchuk predicates.
  const Point n = GeometryTools::cross_product(p0, p1, p2);
  const double nx = std::abs(n[0]);
  const double ny = std::abs(n[1]);
  const double nz = std::abs(n[2]);

  // Choose the two coordinates to keep (drop the dominant one).
  std::size_t keep0, keep1;  // coordinate indices to keep after projection
  if (nz >= nx && nz >= ny)      { keep0 = 0; keep1 = 1; } // drop z
  else if (ny >= nx && ny >= nz) { keep0 = 0; keep1 = 2; } // drop y
  else                           { keep0 = 1; keep1 = 2; } // drop x

  // Apply orient2d in the projected plane.
  const double a0[2] = {p0[keep0],    p0[keep1]   };
  const double a1[2] = {p1[keep0],    p1[keep1]   };
  const double a2[2] = {p2[keep0],    p2[keep1]   };
  const double qp[2] = {point[keep0], point[keep1]};

  const double ref2d = _orient2d(a0, a1, a2);
  const double e01   = _orient2d(a0, a1, qp);
  const double e12   = _orient2d(a1, a2, qp);
  const double e20   = _orient2d(a2, a0, qp);

  if (ref2d > 0.0) return e01 >= 0.0 && e12 >= 0.0 && e20 >= 0.0;
  if (ref2d < 0.0) return e01 <= 0.0 && e12 <= 0.0 && e20 <= 0.0;
  // Degenerate projected triangle (co-projected): fall back to vertex check.
  return (p0 == point || p1 == point || p2 == point);
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_triangle_point_3d(const Point& p0,
                                                      const Point& p1,
                                                      const Point& p2,
                                                      const Point& point)
{
  if (p0 == point or p1 == point or p2 == point)
    return true;

  const double tet_det = orient3d(p0, p1, p2, point);

  if (tet_det < 0.0 or tet_det > 0.0)
    return false;

  return _collides_triangle_point_3d_in_plane(p0, p1, p2, point);
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_triangle_segment_2d(const Point& p0,
                                                        const Point& p1,
                                                        const Point& p2,
                                                        const Point& q0,
                                                        const Point& q1)
{
  // Pre-compute orient2d values for all edge-point combinations.
  const double ref   = orient2d(p0, p1, p2);
  const double o01_q0 = orient2d(p0, p1, q0);
  const double o01_q1 = orient2d(p0, p1, q1);
  const double o12_q0 = orient2d(p1, p2, q0);
  const double o12_q1 = orient2d(p1, p2, q1);
  const double o20_q0 = orient2d(p2, p0, q0);
  const double o20_q1 = orient2d(p2, p0, q1);

  // Check if q0 or q1 are inside the triangle.
  if (_collides_triangle_point_2d_with_hint(p0, p1, p2, q0, ref, o01_q0, o12_q0, o20_q0))
    return true;
  if (_collides_triangle_point_2d_with_hint(p0, p1, p2, q1, ref, o01_q1, o12_q1, o20_q1))
    return true;

  // Check if any triangle edge is crossed by the segment.
  // orient2d(pi, pj, q_k) and orient2d(q0, q1, p_k) are already computed.
  // orient2d(p0,p2,q_k) = -orient2d(p2,p0,q_k) = -o20_q_k

  const double oq_p0  = orient2d(q0, q1, p0);
  const double oq_p1  = orient2d(q0, q1, p1);
  const double oq_p2  = orient2d(q0, q1, p2);

  if (_collides_segment_segment_2d_with_hint(p0, p1, q0, q1, o01_q0, o01_q1, oq_p0, oq_p1))
    return true;
  if (_collides_segment_segment_2d_with_hint(p0, p2, q0, q1, -o20_q0, -o20_q1, oq_p0, oq_p2))
    return true;
  if (_collides_segment_segment_2d_with_hint(p1, p2, q0, q1, o12_q0, o12_q1, oq_p1, oq_p2))
    return true;

  return false;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_triangle_segment_3d_with_hint(const Point& r,
                                                                   const Point& s,
                                                                   const Point& t,
                                                                   const Point& a,
                                                                   const Point& b,
                                                                   double rsta,
                                                                   double rstb)
{
  // rsta = orient3d(r,s,t,a) and rstb = orient3d(r,s,t,b) are pre-computed.

  // a and b on the same strict side of plane rst -> no intersection
  if ((rsta < 0.0 && rstb < 0.0) || (rsta > 0.0 && rstb > 0.0))
    return false;

  // Check endpoints: if in the plane (orient==0), test if inside triangle
  if (rsta == 0.0 && _collides_triangle_point_3d_in_plane(r, s, t, a))
    return true;
  if (rstb == 0.0 && _collides_triangle_point_3d_in_plane(r, s, t, b))
    return true;

  if (rsta == 0.0 && rstb == 0.0)
  {
    // Both on the plane but neither inside the triangle: check edge-edge crossings
    if (collides_segment_segment_3d(r, s, a, b)) return true;
    if (collides_segment_segment_3d(r, t, a, b)) return true;
    if (collides_segment_segment_3d(s, t, a, b)) return true;
    return false;
  }

  // a and b on strictly different sides: use tetrahedra orientation test.
  // Orient so that the "a" side is always the positive side (no copy needed).
  const Point& A = (rsta < 0.0) ? b : a;
  const Point& B = (rsta < 0.0) ? a : b;
  if (orient3d(r, A, s, B) < 0) return false;
  if (orient3d(s, A, t, B) < 0) return false;
  if (orient3d(t, A, r, B) < 0) return false;
  return true;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_triangle_segment_3d(const Point& r,
                                                        const Point& s,
                                                        const Point& t,
                                                        const Point& a,
                                                        const Point& b)
{
  const double rsta = orient3d(r, s, t, a);
  const double rstb = orient3d(r, s, t, b);
  return _collides_triangle_segment_3d_with_hint(r, s, t, a, b, rsta, rstb);
}
//------------------------------------------------------------------------------
bool CollisionPredicates::_collides_triangle_triangle_2d(const Point& p0,
                                                         const Point& p1,
                                                         const Point& p2,
                                                         const Point& q0,
                                                         const Point& q1,
                                                         const Point& q2)
{
  // Pre-compute orient2d of each triangle's edges applied to each vertex of
  // the other triangle.  These 18 values replace all redundant calls.
  const double ref_p = orient2d(p0, p1, p2);
  const double ref_q = orient2d(q0, q1, q2);

  // orient2d(p-edges) applied to q-vertices:
  // fp01_q_i = orient2d(p0,p1, q_i), etc.
  const double fp01_q0 = orient2d(p0, p1, q0), fp01_q1 = orient2d(p0, p1, q1), fp01_q2 = orient2d(p0, p1, q2);
  const double fp12_q0 = orient2d(p1, p2, q0), fp12_q1 = orient2d(p1, p2, q1), fp12_q2 = orient2d(p1, p2, q2);
  const double fp20_q0 = orient2d(p2, p0, q0), fp20_q1 = orient2d(p2, p0, q1), fp20_q2 = orient2d(p2, p0, q2);

  // orient2d(q-edges) applied to p-vertices:
  const double fq01_p0 = orient2d(q0, q1, p0), fq01_p1 = orient2d(q0, q1, p1), fq01_p2 = orient2d(q0, q1, p2);
  const double fq12_p0 = orient2d(q1, q2, p0), fq12_p1 = orient2d(q1, q2, p1), fq12_p2 = orient2d(q1, q2, p2);
  const double fq20_p0 = orient2d(q2, q0, p0), fq20_p1 = orient2d(q2, q0, p1), fq20_p2 = orient2d(q2, q0, p2);

  // Check each q-vertex inside triangle p.
  if (_collides_triangle_point_2d_with_hint(p0,p1,p2, q0, ref_p, fp01_q0, fp12_q0, fp20_q0)) return true;
  if (_collides_triangle_point_2d_with_hint(p0,p1,p2, q1, ref_p, fp01_q1, fp12_q1, fp20_q1)) return true;
  if (_collides_triangle_point_2d_with_hint(p0,p1,p2, q2, ref_p, fp01_q2, fp12_q2, fp20_q2)) return true;

  // Check each p-vertex inside triangle q.
  if (_collides_triangle_point_2d_with_hint(q0,q1,q2, p0, ref_q, fq01_p0, fq12_p0, fq20_p0)) return true;
  if (_collides_triangle_point_2d_with_hint(q0,q1,q2, p1, ref_q, fq01_p1, fq12_p1, fq20_p1)) return true;
  if (_collides_triangle_point_2d_with_hint(q0,q1,q2, p2, ref_q, fq01_p2, fq12_p2, fq20_p2)) return true;

  // Check edge-edge collisions using the pre-computed orientations.
  // fp01_qi = orient2d(p0,p1, q_i); fp12_qi = orient2d(p1,p2, q_i);
  // fp20_qi = orient2d(p2,p0, q_i)  [note: orient2d(p0,p2,.) = -fp20_.]
  // fq01_pi = orient2d(q0,q1, p_i); fq12_pi = orient2d(q1,q2, p_i);
  // fq20_pi = orient2d(q2,q0, p_i)  [note: orient2d(q0,q2,.) = -fq20_.]
  // Edges of p: (p0,p1), (p1,p2), (p2,p0)
  // Edges of q: (q0,q1), (q1,q2), (q2,q0)
  // p-edge (p0,p1) vs q-edges
  if (_collides_segment_segment_2d_with_hint(p0,p1, q0,q1, fp01_q0, fp01_q1, fq01_p0, fq01_p1)) return true;
  if (_collides_segment_segment_2d_with_hint(p0,p1, q1,q2, fp01_q1, fp01_q2, fq12_p0, fq12_p1)) return true;
  if (_collides_segment_segment_2d_with_hint(p0,p1, q2,q0, fp01_q2, fp01_q0, -fq20_p0, -fq20_p1)) return true;
  // p-edge (p1,p2) vs q-edges
  if (_collides_segment_segment_2d_with_hint(p1,p2, q0,q1, fp12_q0, fp12_q1, fq01_p1, fq01_p2)) return true;
  if (_collides_segment_segment_2d_with_hint(p1,p2, q1,q2, fp12_q1, fp12_q2, fq12_p1, fq12_p2)) return true;
  if (_collides_segment_segment_2d_with_hint(p1,p2, q2,q0, fp12_q2, fp12_q0, -fq20_p1, -fq20_p2)) return true;
  // p-edge (p2,p0) vs q-edges
  if (_collides_segment_segment_2d_with_hint(p2,p0, q0,q1, -fp20_q0, -fp20_q1, fq01_p2, fq01_p0)) return true;
  if (_collides_segment_segment_2d_with_hint(p2,p0, q1,q2, -fp20_q1, -fp20_q2, fq12_p2, fq12_p0)) return true;
  if (_collides_segment_segment_2d_with_hint(p2,p0, q2,q0, -fp20_q2, -fp20_q0, -fq20_p2, -fq20_p0)) return true;

  return false;
}
bool CollisionPredicates::_collides_triangle_triangle_3d(const Point& p0,
                                                         const Point& p1,
                                                         const Point& p2,
                                                         const Point& q0,
                                                         const Point& q1,
                                                         const Point& q2)
{
  // Pre-compute orient3d of the p-plane for each q-vertex and vice versa.
  // This eliminates all redundant orient3d calls in the edge-face tests.
  const double pq0 = orient3d(p0, p1, p2, q0);
  const double pq1 = orient3d(p0, p1, p2, q1);
  const double pq2 = orient3d(p0, p1, p2, q2);

  const double qp0 = orient3d(q0, q1, q2, p0);
  const double qp1 = orient3d(q0, q1, q2, p1);
  const double qp2 = orient3d(q0, q1, q2, p2);

  // opposite_sides: returns true iff the face-edge test is worth running
  // (at least one endpoint on the plane, or endpoints on opposite sides).
  auto opposite_sides = [](double s0, double s1) noexcept {
    return (s0 == 0.0) || (s1 == 0.0) || ((s0 > 0.0) != (s1 > 0.0));
  };

  // Check each edge of q against the p-face, using pre-computed hints.
  if (opposite_sides(pq0, pq1) &&
      _collides_triangle_segment_3d_with_hint(p0, p1, p2, q0, q1, pq0, pq1)) return true;
  if (opposite_sides(pq1, pq2) &&
      _collides_triangle_segment_3d_with_hint(p0, p1, p2, q1, q2, pq1, pq2)) return true;
  if (opposite_sides(pq2, pq0) &&
      _collides_triangle_segment_3d_with_hint(p0, p1, p2, q2, q0, pq2, pq0)) return true;

  // Check each edge of p against the q-face, using pre-computed hints.
  if (opposite_sides(qp0, qp1) &&
      _collides_triangle_segment_3d_with_hint(q0, q1, q2, p0, p1, qp0, qp1)) return true;
  if (opposite_sides(qp1, qp2) &&
      _collides_triangle_segment_3d_with_hint(q0, q1, q2, p1, p2, qp1, qp2)) return true;
  if (opposite_sides(qp2, qp0) &&
      _collides_triangle_segment_3d_with_hint(q0, q1, q2, p2, p0, qp2, qp0)) return true;

  return false;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_tetrahedron_point_3d(const Point& p0,
                                                         const Point& p1,
                                                         const Point& p2,
                                                         const Point& p3,
                                                         const Point& point)
{
  const double ref = orient3d(p0, p1, p2, p3);

  if (ref > 0.0)
    {
      return (orient3d(p0, p1, p2, point) >= 0.0 and
	      orient3d(p0, p3, p1, point) >= 0.0 and
	      orient3d(p0, p2, p3, point) >= 0.0 and
	      orient3d(p1, p3, p2, point) >= 0.0);
    }
  else if (ref < 0.0)
    {
      return (orient3d(p0, p1, p2, point) <= 0.0 and
	      orient3d(p0, p3, p1, point) <= 0.0 and
	      orient3d(p0, p2, p3, point) <= 0.0 and
	      orient3d(p1, p3, p2, point) <= 0.0);
    }
  else
    {
      throw std::runtime_error("Not implemented for degenerate tetrahedron");
    }

  return false;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_tetrahedron_segment_3d(const Point& p0,
                                                           const Point& p1,
                                                           const Point& p2,
                                                           const Point& p3,
                                                           const Point& q0,
                                                           const Point& q1)
{
  const double ref = orient3d(p0, p1, p2, p3);
  const bool pos = (ref > 0.0);

  auto inside_from_fo = [&](double fo0, double fo1, double fo2, double fo3) -> bool
  {
    if (pos) return (fo0 >= 0.0 && fo1 >= 0.0 && fo2 >= 0.0 && fo3 >= 0.0);
    else     return (fo0 <= 0.0 && fo1 <= 0.0 && fo2 <= 0.0 && fo3 <= 0.0);
  };

  // Phase 1: q0
  const double fo0_q0 = orient3d(p0, p1, p2, q0);
  const double fo1_q0 = orient3d(p0, p3, p1, q0);
  const double fo2_q0 = orient3d(p0, p2, p3, q0);
  const double fo3_q0 = orient3d(p1, p3, p2, q0);
  if (inside_from_fo(fo0_q0, fo1_q0, fo2_q0, fo3_q0)) return true;

  // Phase 2: q1
  const double fo0_q1 = orient3d(p0, p1, p2, q1);
  const double fo1_q1 = orient3d(p0, p3, p1, q1);
  const double fo2_q1 = orient3d(p0, p2, p3, q1);
  const double fo3_q1 = orient3d(p1, p3, p2, q1);
  if (inside_from_fo(fo0_q1, fo1_q1, fo2_q1, fo3_q1)) return true;

  // Phase 3: face crossings -- skip if both endpoints are on the same strict side.
  auto opp = [](double s0, double s1) noexcept {
    return (s0 == 0.0) || (s1 == 0.0) || ((s0 > 0.0) != (s1 > 0.0));
  };
  if (opp(-fo3_q0, -fo3_q1) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q0,q1, -fo3_q0, -fo3_q1)) return true;
  if (opp( fo2_q0,  fo2_q1) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q0,q1,  fo2_q0,  fo2_q1)) return true;
  if (opp(-fo1_q0, -fo1_q1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q0,q1, -fo1_q0, -fo1_q1)) return true;
  if (opp( fo0_q0,  fo0_q1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q0,q1,  fo0_q0,  fo0_q1)) return true;

  return false;
}
//-----------------------------------------------------------------------------
bool CollisionPredicates::_collides_tetrahedron_triangle_3d(const Point& p0,
							    const Point& p1,
							    const Point& p2,
							    const Point& p3,
							    const Point& q0,
							    const Point& q1,
							    const Point& q2)
{
  const double ref = orient3d(p0, p1, p2, p3);
  const bool pos = (ref > 0.0);

  auto inside = [&](double fo0, double fo1, double fo2, double fo3) -> bool
  {
    if (pos) return (fo0 >= 0.0 && fo1 >= 0.0 && fo2 >= 0.0 && fo3 >= 0.0);
    else     return (fo0 <= 0.0 && fo1 <= 0.0 && fo2 <= 0.0 && fo3 <= 0.0);
  };

  // Compute face-orientations for q0, q1, q2 only as needed.
  struct FO { double fo0, fo1, fo2, fo3; };

  const auto eval_q = [&](const Point& q) -> FO {
    return FO{
      orient3d(p0,p1,p2,q),
      orient3d(p0,p3,p1,q),
      orient3d(p0,p2,p3,q),
      orient3d(p1,p3,p2,q)
    };
  };

  const FO fq0 = eval_q(q0);
  if (inside(fq0.fo0, fq0.fo1, fq0.fo2, fq0.fo3)) return true;

  const FO fq1 = eval_q(q1);
  if (inside(fq1.fo0, fq1.fo1, fq1.fo2, fq1.fo3)) return true;

  const FO fq2 = eval_q(q2);
  if (inside(fq2.fo0, fq2.fo1, fq2.fo2, fq2.fo3)) return true;

  auto opp = [](double s0, double s1) noexcept {
    return (s0 == 0.0) || (s1 == 0.0) || ((s0 > 0.0) != (s1 > 0.0));
  };

  // Face (p1,p2,p3): plane-orient = -fo3_*
  if (opp(-fq0.fo3, -fq1.fo3) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q0,q1, -fq0.fo3, -fq1.fo3)) return true;
  if (opp(-fq1.fo3, -fq2.fo3) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q1,q2, -fq1.fo3, -fq2.fo3)) return true;
  if (opp(-fq0.fo3, -fq2.fo3) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q0,q2, -fq0.fo3, -fq2.fo3)) return true;

  // Face (p0,p2,p3): plane-orient =  fo2_*
  if (opp( fq0.fo2,  fq1.fo2) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q0,q1,  fq0.fo2,  fq1.fo2)) return true;
  if (opp( fq1.fo2,  fq2.fo2) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q1,q2,  fq1.fo2,  fq2.fo2)) return true;
  if (opp( fq0.fo2,  fq2.fo2) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q0,q2,  fq0.fo2,  fq2.fo2)) return true;

  // Face (p0,p1,p3): plane-orient = -fo1_*
  if (opp(-fq0.fo1, -fq1.fo1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q0,q1, -fq0.fo1, -fq1.fo1)) return true;
  if (opp(-fq1.fo1, -fq2.fo1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q1,q2, -fq1.fo1, -fq2.fo1)) return true;
  if (opp(-fq0.fo1, -fq2.fo1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q0,q2, -fq0.fo1, -fq2.fo1)) return true;

  // Face (p0,p1,p2): plane-orient =  fo0_*
  if (opp( fq0.fo0,  fq1.fo0) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q0,q1,  fq0.fo0,  fq1.fo0)) return true;
  if (opp( fq1.fo0,  fq2.fo0) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q1,q2,  fq1.fo0,  fq2.fo0)) return true;
  if (opp( fq0.fo0,  fq2.fo0) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q0,q2,  fq0.fo0,  fq2.fo0)) return true;

  // Tet-edge vs triangle-face: compute orient3d for each tet vertex w.r.t q-face.
  const double qf_p0 = orient3d(q0,q1,q2,p0);
  const double qf_p1 = orient3d(q0,q1,q2,p1);
  const double qf_p2 = orient3d(q0,q1,q2,p2);
  const double qf_p3 = orient3d(q0,q1,q2,p3);

  if (opp(qf_p0,qf_p1) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p0,p1, qf_p0,qf_p1)) return true;
  if (opp(qf_p0,qf_p2) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p0,p2, qf_p0,qf_p2)) return true;
  if (opp(qf_p0,qf_p3) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p0,p3, qf_p0,qf_p3)) return true;
  if (opp(qf_p1,qf_p2) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p1,p2, qf_p1,qf_p2)) return true;
  if (opp(qf_p1,qf_p3) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p1,p3, qf_p1,qf_p3)) return true;
  if (opp(qf_p2,qf_p3) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p2,p3, qf_p2,qf_p3)) return true;

  return false;
}


//-----------------------------------------------------------------------------

struct FO { double f0,f1,f2,f3; };

static inline bool inside_fo(double ref, const FO& fo) noexcept
{
  if (ref > 0.0) return (fo.f0>=0 && fo.f1>=0 && fo.f2>=0 && fo.f3>=0);
  else           return (fo.f0<=0 && fo.f1<=0 && fo.f2<=0 && fo.f3<=0);
}

static inline FO faces_of_p_at_q(const Point& p0,const Point& p1,const Point& p2,const Point& p3,
				 const Point& q)
{
  return FO{
    orient3d(p0,p1,p2,q),
    orient3d(p0,p3,p1,q),
    orient3d(p0,p2,p3,q),
    orient3d(p1,p3,p2,q)
  };
}

bool CollisionPredicates::_collides_tetrahedron_tetrahedron_3d(const Point& p0,
							       const Point& p1,
							       const Point& p2,
							       const Point& p3,
							       const Point& q0,
							       const Point& q1,
							       const Point& q2,
							       const Point& q3)
{
// Phase A: test q-vertices inside p (compute per vertex, early-out)
  const double ref_p = orient3d(p0,p1,p2,p3);

  const FO p_at_q0 = faces_of_p_at_q(p0,p1,p2,p3,q0); if (inside_fo(ref_p, p_at_q0)) return true;
  const FO p_at_q1 = faces_of_p_at_q(p0,p1,p2,p3,q1); if (inside_fo(ref_p, p_at_q1)) return true;
  const FO p_at_q2 = faces_of_p_at_q(p0,p1,p2,p3,q2); if (inside_fo(ref_p, p_at_q2)) return true;
  const FO p_at_q3 = faces_of_p_at_q(p0,p1,p2,p3,q3); if (inside_fo(ref_p, p_at_q3)) return true;

  // Phase B: test p-vertices inside q (same)
  const double ref_q = orient3d(q0,q1,q2,q3);

  // note: face order for q in your code: fq0=orient(q0,q1,q2,.), fq1=orient(q0,q3,q1,.), ...
  auto faces_of_q_at_p = [&](const Point& p) -> FO {
    return FO{
      orient3d(q0,q1,q2,p),
      orient3d(q0,q3,q1,p),
      orient3d(q0,q2,q3,p),
      orient3d(q1,q3,q2,p)
    };
  };

  const FO q_at_p0 = faces_of_q_at_p(p0); if (inside_fo(ref_q, q_at_p0)) return true;
  const FO q_at_p1 = faces_of_q_at_p(p1); if (inside_fo(ref_q, q_at_p1)) return true;
  const FO q_at_p2 = faces_of_q_at_p(p2); if (inside_fo(ref_q, q_at_p2)) return true;
  const FO q_at_p3 = faces_of_q_at_p(p3); if (inside_fo(ref_q, q_at_p3)) return true;

  // Phase C: face-edge tests using the FOs already computed above.
  // Skip each test when both q-edge endpoints lie strictly on the same side.
  auto opp = [](double s0, double s1) noexcept {
    return (s0 == 0.0) || (s1 == 0.0) || ((s0 > 0.0) != (s1 > 0.0));
  };

  // 4 p-faces vs 6 q-edges
  // p-face (p1,p2,p3): plane-orient for qi = -p_at_qi.f3
  if (opp(-p_at_q0.f3, -p_at_q1.f3) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q0,q1, -p_at_q0.f3,-p_at_q1.f3)) return true;
  if (opp(-p_at_q0.f3, -p_at_q2.f3) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q0,q2, -p_at_q0.f3,-p_at_q2.f3)) return true;
  if (opp(-p_at_q0.f3, -p_at_q3.f3) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q0,q3, -p_at_q0.f3,-p_at_q3.f3)) return true;
  if (opp(-p_at_q1.f3, -p_at_q2.f3) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q1,q2, -p_at_q1.f3,-p_at_q2.f3)) return true;
  if (opp(-p_at_q1.f3, -p_at_q3.f3) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q1,q3, -p_at_q1.f3,-p_at_q3.f3)) return true;
  if (opp(-p_at_q2.f3, -p_at_q3.f3) &&
      _collides_triangle_segment_3d_with_hint(p1,p2,p3, q2,q3, -p_at_q2.f3,-p_at_q3.f3)) return true;

  // p-face (p0,p2,p3): plane-orient for qi =  p_at_qi.f2
  if (opp( p_at_q0.f2,  p_at_q1.f2) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q0,q1,  p_at_q0.f2, p_at_q1.f2)) return true;
  if (opp( p_at_q0.f2,  p_at_q2.f2) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q0,q2,  p_at_q0.f2, p_at_q2.f2)) return true;
  if (opp( p_at_q0.f2,  p_at_q3.f2) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q0,q3,  p_at_q0.f2, p_at_q3.f2)) return true;
  if (opp( p_at_q1.f2,  p_at_q2.f2) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q1,q2,  p_at_q1.f2, p_at_q2.f2)) return true;
  if (opp( p_at_q1.f2,  p_at_q3.f2) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q1,q3,  p_at_q1.f2, p_at_q3.f2)) return true;
  if (opp( p_at_q2.f2,  p_at_q3.f2) &&
      _collides_triangle_segment_3d_with_hint(p0,p2,p3, q2,q3,  p_at_q2.f2, p_at_q3.f2)) return true;

  // p-face (p0,p1,p3): plane-orient for qi = -p_at_qi.f1
  if (opp(-p_at_q0.f1, -p_at_q1.f1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q0,q1, -p_at_q0.f1,-p_at_q1.f1)) return true;
  if (opp(-p_at_q0.f1, -p_at_q2.f1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q0,q2, -p_at_q0.f1,-p_at_q2.f1)) return true;
  if (opp(-p_at_q0.f1, -p_at_q3.f1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q0,q3, -p_at_q0.f1,-p_at_q3.f1)) return true;
  if (opp(-p_at_q1.f1, -p_at_q2.f1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q1,q2, -p_at_q1.f1,-p_at_q2.f1)) return true;
  if (opp(-p_at_q1.f1, -p_at_q3.f1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q1,q3, -p_at_q1.f1,-p_at_q3.f1)) return true;
  if (opp(-p_at_q2.f1, -p_at_q3.f1) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p3, q2,q3, -p_at_q2.f1,-p_at_q3.f1)) return true;

  // p-face (p0,p1,p2): plane-orient for qi =  p_at_qi.f0
  if (opp( p_at_q0.f0,  p_at_q1.f0) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q0,q1,  p_at_q0.f0, p_at_q1.f0)) return true;
  if (opp( p_at_q0.f0,  p_at_q2.f0) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q0,q2,  p_at_q0.f0, p_at_q2.f0)) return true;
  if (opp( p_at_q0.f0,  p_at_q3.f0) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q0,q3,  p_at_q0.f0, p_at_q3.f0)) return true;
  if (opp( p_at_q1.f0,  p_at_q2.f0) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q1,q2,  p_at_q1.f0, p_at_q2.f0)) return true;
  if (opp( p_at_q1.f0,  p_at_q3.f0) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q1,q3,  p_at_q1.f0, p_at_q3.f0)) return true;
  if (opp( p_at_q2.f0,  p_at_q3.f0) &&
      _collides_triangle_segment_3d_with_hint(p0,p1,p2, q2,q3,  p_at_q2.f0, p_at_q3.f0)) return true;

  // 4 q-faces vs 6 p-edges
  // q-face (q1,q2,q3): plane-orient for pi = -q_at_pi.f3
  if (opp(-q_at_p0.f3, -q_at_p1.f3) &&
      _collides_triangle_segment_3d_with_hint(q1,q2,q3, p0,p1, -q_at_p0.f3,-q_at_p1.f3)) return true;
  if (opp(-q_at_p0.f3, -q_at_p2.f3) &&
      _collides_triangle_segment_3d_with_hint(q1,q2,q3, p0,p2, -q_at_p0.f3,-q_at_p2.f3)) return true;
  if (opp(-q_at_p0.f3, -q_at_p3.f3) &&
      _collides_triangle_segment_3d_with_hint(q1,q2,q3, p0,p3, -q_at_p0.f3,-q_at_p3.f3)) return true;
  if (opp(-q_at_p1.f3, -q_at_p2.f3) &&
      _collides_triangle_segment_3d_with_hint(q1,q2,q3, p1,p2, -q_at_p1.f3,-q_at_p2.f3)) return true;
  if (opp(-q_at_p1.f3, -q_at_p3.f3) &&
      _collides_triangle_segment_3d_with_hint(q1,q2,q3, p1,p3, -q_at_p1.f3,-q_at_p3.f3)) return true;
  if (opp(-q_at_p2.f3, -q_at_p3.f3) &&
      _collides_triangle_segment_3d_with_hint(q1,q2,q3, p2,p3, -q_at_p2.f3,-q_at_p3.f3)) return true;

  // q-face (q0,q2,q3): plane-orient for pi =  q_at_pi.f2
  if (opp( q_at_p0.f2,  q_at_p1.f2) &&
      _collides_triangle_segment_3d_with_hint(q0,q2,q3, p0,p1,  q_at_p0.f2, q_at_p1.f2)) return true;
  if (opp( q_at_p0.f2,  q_at_p2.f2) &&
      _collides_triangle_segment_3d_with_hint(q0,q2,q3, p0,p2,  q_at_p0.f2, q_at_p2.f2)) return true;
  if (opp( q_at_p0.f2,  q_at_p3.f2) &&
      _collides_triangle_segment_3d_with_hint(q0,q2,q3, p0,p3,  q_at_p0.f2, q_at_p3.f2)) return true;
  if (opp( q_at_p1.f2,  q_at_p2.f2) &&
      _collides_triangle_segment_3d_with_hint(q0,q2,q3, p1,p2,  q_at_p1.f2, q_at_p2.f2)) return true;
  if (opp( q_at_p1.f2,  q_at_p3.f2) &&
      _collides_triangle_segment_3d_with_hint(q0,q2,q3, p1,p3,  q_at_p1.f2, q_at_p3.f2)) return true;
  if (opp( q_at_p2.f2,  q_at_p3.f2) &&
      _collides_triangle_segment_3d_with_hint(q0,q2,q3, p2,p3,  q_at_p2.f2, q_at_p3.f2)) return true;

  // q-face (q0,q1,q3): plane-orient for pi = -q_at_pi.f1
  if (opp(-q_at_p0.f1, -q_at_p1.f1) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q3, p0,p1, -q_at_p0.f1,-q_at_p1.f1)) return true;
  if (opp(-q_at_p0.f1, -q_at_p2.f1) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q3, p0,p2, -q_at_p0.f1,-q_at_p2.f1)) return true;
  if (opp(-q_at_p0.f1, -q_at_p3.f1) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q3, p0,p3, -q_at_p0.f1,-q_at_p3.f1)) return true;
  if (opp(-q_at_p1.f1, -q_at_p2.f1) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q3, p1,p2, -q_at_p1.f1,-q_at_p2.f1)) return true;
  if (opp(-q_at_p1.f1, -q_at_p3.f1) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q3, p1,p3, -q_at_p1.f1,-q_at_p3.f1)) return true;
  if (opp(-q_at_p2.f1, -q_at_p3.f1) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q3, p2,p3, -q_at_p2.f1,-q_at_p3.f1)) return true;

  // q-face (q0,q1,q2): plane-orient for pi =  q_at_pi.f0
  if (opp( q_at_p0.f0,  q_at_p1.f0) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p0,p1,  q_at_p0.f0, q_at_p1.f0)) return true;
  if (opp( q_at_p0.f0,  q_at_p2.f0) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p0,p2,  q_at_p0.f0, q_at_p2.f0)) return true;
  if (opp( q_at_p0.f0,  q_at_p3.f0) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p0,p3,  q_at_p0.f0, q_at_p3.f0)) return true;
  if (opp( q_at_p1.f0,  q_at_p2.f0) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p1,p2,  q_at_p1.f0, q_at_p2.f0)) return true;
  if (opp( q_at_p1.f0,  q_at_p3.f0) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p1,p3,  q_at_p1.f0, q_at_p3.f0)) return true;
  if (opp( q_at_p2.f0,  q_at_p3.f0) &&
      _collides_triangle_segment_3d_with_hint(q0,q1,q2, p2,p3,  q_at_p2.f0, q_at_p3.f0)) return true;

  return false;
}


//-----------------------------------------------------------------------------

// ============================================================================
// dolfin::CollisionPredicates -- high-level MeshEntity wrappers
// ============================================================================
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/MeshGeometry.h>

bool dolfin::CollisionPredicates::collides(const dolfin::MeshEntity& entity,
                                           const simpex::Point& point)
{
  // Only simplex meshes are supported
  const dolfin::MeshGeometry& g = entity.mesh().geometry();
  const unsigned int* v = entity.entities(0);
  const std::size_t tdim = entity.mesh().topology().dim();
  const std::size_t gdim = g.dim();

  if (tdim == 1 && gdim == 1)
    return collides_segment_point_1d(g.point(v[0])[0], g.point(v[1])[0], point[0]);

  if (tdim == 1 && gdim == 2)
    return collides_segment_point_2d(g.point(v[0]), g.point(v[1]), point);

  if (tdim == 1 && gdim == 3)
    return collides_segment_point_3d(g.point(v[0]), g.point(v[1]), point);

  if (tdim == 2 && gdim == 2)
    return collides_triangle_point_2d(g.point(v[0]), g.point(v[1]),
                                      g.point(v[2]), point);

  if (tdim == 2 && gdim == 3)
    return collides_triangle_point_3d(g.point(v[0]), g.point(v[1]),
                                      g.point(v[2]), point);

  if (tdim == 3)
    return collides_tetrahedron_point_3d(g.point(v[0]), g.point(v[1]),
                                         g.point(v[2]), g.point(v[3]), point);

  return false;
}
//-----------------------------------------------------------------------------
bool dolfin::CollisionPredicates::collides(const dolfin::MeshEntity& entity_0,
                                           const dolfin::MeshEntity& entity_1)
{
  const dolfin::MeshGeometry& g0 = entity_0.mesh().geometry();
  const dolfin::MeshGeometry& g1 = entity_1.mesh().geometry();
  const unsigned int* v0 = entity_0.entities(0);
  const unsigned int* v1 = entity_1.entities(0);
  const std::size_t d0 = entity_0.dim();
  const std::size_t d1 = entity_1.dim();
  const std::size_t gdim = g0.dim();

  if (d0 == 1 && d1 == 1)
    return collides_segment_segment(g0.point(v0[0]), g0.point(v0[1]),
                                    g1.point(v1[0]), g1.point(v1[1]), gdim);

  if (d0 == 1 && d1 == 2)
    return collides_triangle_segment(g1.point(v1[0]), g1.point(v1[1]),
                                     g1.point(v1[2]),
                                     g0.point(v0[0]), g0.point(v0[1]), gdim);

  if (d0 == 2 && d1 == 1)
    return collides_triangle_segment(g0.point(v0[0]), g0.point(v0[1]),
                                     g0.point(v0[2]),
                                     g1.point(v1[0]), g1.point(v1[1]), gdim);

  if (d0 == 2 && d1 == 2)
    return collides_triangle_triangle(g0.point(v0[0]), g0.point(v0[1]),
                                      g0.point(v0[2]),
                                      g1.point(v1[0]), g1.point(v1[1]),
                                      g1.point(v1[2]), gdim);

  if (d0 == 2 && d1 == 3)
    return collides_tetrahedron_triangle_3d(g1.point(v1[0]), g1.point(v1[1]),
                                            g1.point(v1[2]), g1.point(v1[3]),
                                            g0.point(v0[0]), g0.point(v0[1]),
                                            g0.point(v0[2]));

  if (d0 == 3 && d1 == 2)
    return collides_tetrahedron_triangle_3d(g0.point(v0[0]), g0.point(v0[1]),
                                            g0.point(v0[2]), g0.point(v0[3]),
                                            g1.point(v1[0]), g1.point(v1[1]),
                                            g1.point(v1[2]));

  if (d0 == 3 && d1 == 3)
    return collides_tetrahedron_tetrahedron_3d(g0.point(v0[0]), g0.point(v0[1]),
                                               g0.point(v0[2]), g0.point(v0[3]),
                                               g1.point(v1[0]), g1.point(v1[1]),
                                               g1.point(v1[2]), g1.point(v1[3]));

  return false;
}
//-----------------------------------------------------------------------------
