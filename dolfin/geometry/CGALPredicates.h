// CGALPredicates.h
//
// Provides well-optimized CGAL-based geometry predicates that mirror the
// Shewchuk API in CollisionPredicates.h / predicates.h.
//
// Available only when compiled with -DSIMPEX_WITH_CGAL.
// Uses CGAL::Exact_predicates_inexact_constructions_kernel (EPICK) throughout.
// All collision tests use CGAL::do_intersect -- no general Polyhedron routines.

#ifndef SIMPEX_CGAL_PREDICATES_H
#define SIMPEX_CGAL_PREDICATES_H

#ifdef SIMPEX_WITH_CGAL

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>

#include "Point.h"

namespace dolfin {

// ---------------------------------------------------------------------------
// Internal helpers: convert simpex::Point to CGAL EPICK types
// ---------------------------------------------------------------------------
namespace _cgal {
  using K    = CGAL::Exact_predicates_inexact_constructions_kernel;
  using P2   = K::Point_2;
  using P3   = K::Point_3;
  using Seg2 = K::Segment_2;
  using Seg3 = K::Segment_3;
  using Tri2 = K::Triangle_2;
  using Tri3 = K::Triangle_3;
  using Tet3 = K::Tetrahedron_3;

  inline P2   to2(const Point& p) { return P2(p.x(), p.y()); }
  inline P3   to3(const Point& p) { return P3(p.x(), p.y(), p.z()); }
  inline Seg2 seg2(const Point& a, const Point& b) { return Seg2(to2(a), to2(b)); }
  inline Seg3 seg3(const Point& a, const Point& b) { return Seg3(to3(a), to3(b)); }
  inline Tri2 tri2(const Point& a, const Point& b, const Point& c)
  { return Tri2(to2(a), to2(b), to2(c)); }
  inline Tri3 tri(const Point& a, const Point& b, const Point& c)
  { return Tri3(to3(a), to3(b), to3(c)); }
  inline Tet3 tet(const Point& a, const Point& b, const Point& c, const Point& d)
  { return Tet3(to3(a), to3(b), to3(c), to3(d)); }
}

// ---------------------------------------------------------------------------
// Orientation predicates
// ---------------------------------------------------------------------------

/// 2D orientation using CGAL EPICK.
/// Returns positive (CCW), negative (CW), or zero (collinear).
inline double cgal_orient2d(const Point& a, const Point& b, const Point& c)
{
  return static_cast<double>(
    CGAL::orientation(_cgal::to2(a), _cgal::to2(b), _cgal::to2(c)));
}

/// 3D orientation using CGAL EPICK, adjusted to match Shewchuk's sign convention.
///
/// Shewchuk's orient3d uses the **left-hand rule** (as stated in predicates.h):
///   orient3d(a,b,c,d) > 0  when a,b,c,d are in left-hand orientation.
/// CGAL's orientation uses the **right-hand rule**:
///   CGAL::orientation(a,b,c,d) = POSITIVE when d is on the right-hand side of (a,b,c).
///
/// These conventions are exactly opposite for non-degenerate points, so we negate
/// the CGAL result here to make cgal_orient3d() agree with orient3d() in sign.
inline double cgal_orient3d(const Point& a, const Point& b,
                             const Point& c, const Point& d)
{
  return -static_cast<double>(
    CGAL::orientation(_cgal::to3(a), _cgal::to3(b),
                      _cgal::to3(c), _cgal::to3(d)));
}

// ---------------------------------------------------------------------------
// Collision predicates -- all use CGAL::do_intersect (EPICK)
// ---------------------------------------------------------------------------

/// 2D segment-point collision.
inline bool cgal_collides_segment_point_2d(const Point& p0, const Point& p1, const Point& q)
{ return CGAL::do_intersect(_cgal::seg2(p0, p1), _cgal::to2(q)); }

/// 3D segment-point collision.
inline bool cgal_collides_segment_point_3d(const Point& p0, const Point& p1, const Point& q)
{ return CGAL::do_intersect(_cgal::seg3(p0, p1), _cgal::to3(q)); }

/// 2D segment-segment collision.
inline bool cgal_collides_segment_segment_2d(
    const Point& p0, const Point& p1, const Point& q0, const Point& q1)
{ return CGAL::do_intersect(_cgal::seg2(p0, p1), _cgal::seg2(q0, q1)); }

/// 3D segment-segment collision.
inline bool cgal_collides_segment_segment_3d(
    const Point& p0, const Point& p1, const Point& q0, const Point& q1)
{ return CGAL::do_intersect(_cgal::seg3(p0, p1), _cgal::seg3(q0, q1)); }

/// 2D triangle-point collision.
inline bool cgal_collides_triangle_point_2d(
    const Point& p0, const Point& p1, const Point& p2, const Point& q)
{ return CGAL::do_intersect(_cgal::tri2(p0, p1, p2), _cgal::to2(q)); }

/// 2D triangle-segment collision.
inline bool cgal_collides_triangle_segment_2d(
    const Point& p0, const Point& p1, const Point& p2,
    const Point& q0, const Point& q1)
{ return CGAL::do_intersect(_cgal::tri2(p0, p1, p2), _cgal::seg2(q0, q1)); }

/// 3D triangle-segment collision.
inline bool cgal_collides_triangle_segment_3d(
    const Point& p0, const Point& p1, const Point& p2,
    const Point& q0, const Point& q1)
{ return CGAL::do_intersect(_cgal::tri(p0, p1, p2), _cgal::seg3(q0, q1)); }

/// 2D triangle-triangle collision.
inline bool cgal_collides_triangle_triangle_2d(
    const Point& p0, const Point& p1, const Point& p2,
    const Point& q0, const Point& q1, const Point& q2)
{ return CGAL::do_intersect(_cgal::tri2(p0, p1, p2), _cgal::tri2(q0, q1, q2)); }

/// Check whether two 3D triangles intersect.
inline bool cgal_collides_triangle_triangle_3d(
    const Point& p0, const Point& p1, const Point& p2,
    const Point& q0, const Point& q1, const Point& q2)
{
  return CGAL::do_intersect(_cgal::tri(p0, p1, p2),
                            _cgal::tri(q0, q1, q2));
}

/// Check whether two 3D tetrahedra intersect.
inline bool cgal_collides_tetrahedron_tetrahedron_3d(
    const Point& p0, const Point& p1, const Point& p2, const Point& p3,
    const Point& q0, const Point& q1, const Point& q2, const Point& q3)
{
  return CGAL::do_intersect(_cgal::tet(p0, p1, p2, p3),
                            _cgal::tet(q0, q1, q2, q3));
}

/// Check whether a 3D triangle intersects a point.
inline bool cgal_collides_triangle_point_3d(
    const Point& p0, const Point& p1, const Point& p2, const Point& q)
{
  return CGAL::do_intersect(_cgal::tri(p0, p1, p2), _cgal::to3(q));
}

/// Check whether a 3D tetrahedron intersects a point.
inline bool cgal_collides_tetrahedron_point_3d(
    const Point& p0, const Point& p1, const Point& p2, const Point& p3,
    const Point& q)
{
  return CGAL::do_intersect(_cgal::tet(p0, p1, p2, p3), _cgal::to3(q));
}

/// Check whether a 3D tetrahedron intersects a segment.
inline bool cgal_collides_tetrahedron_segment_3d(
    const Point& p0, const Point& p1, const Point& p2, const Point& p3,
    const Point& q0, const Point& q1)
{
  return CGAL::do_intersect(_cgal::tet(p0, p1, p2, p3), _cgal::seg3(q0, q1));
}

/// Check whether a 3D tetrahedron intersects a triangle.
inline bool cgal_collides_tetrahedron_triangle_3d(
    const Point& p0, const Point& p1, const Point& p2, const Point& p3,
    const Point& q0, const Point& q1, const Point& q2)
{
  return CGAL::do_intersect(_cgal::tet(p0, p1, p2, p3), _cgal::tri(q0, q1, q2));
}

} // namespace dolfin

#endif // SIMPEX_WITH_CGAL
#endif // SIMPEX_CGAL_PREDICATES_H
