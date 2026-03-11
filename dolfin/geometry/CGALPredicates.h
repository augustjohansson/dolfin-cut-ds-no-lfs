// CGALPredicates.h
//
// Provides well-optimized CGAL-based geometry predicates that mirror the
// Shewchuk API in CollisionPredicates.h / predicates.h.
//
// Available only when compiled with -DDOLFIN_WITH_CGAL.
// Uses CGAL::Exact_predicates_inexact_constructions_kernel (EPICK) throughout.
// All collision tests use CGAL::do_intersect -- no general Polyhedron routines.

#ifndef DOLFIN_CGAL_PREDICATES_H
#define DOLFIN_CGAL_PREDICATES_H

#ifdef DOLFIN_WITH_CGAL

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <vector>

#if CGAL_VERSION_NR >= 1060100000
#include <variant>
namespace { template<typename T, typename V> static const T* _cgalp_get_if(const V* v) { return std::get_if<T>(v); } }
#else
#include <boost/variant.hpp>
namespace { template<typename T, typename V> static const T* _cgalp_get_if(const V* v) { return boost::get<T>(v); } }
#endif

#include "Point.h"

namespace dolfin {

// ---------------------------------------------------------------------------
// Internal helpers: convert dolfin::Point to CGAL EPICK types
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

// ---------------------------------------------------------------------------
// Degeneracy predicates
// ---------------------------------------------------------------------------

/// Check if a 2D segment is degenerate (both endpoints coincide).
inline bool is_degenerate_2d(const Point& a, const Point& b)
{ return _cgal::seg2(a, b).is_degenerate(); }

/// Check if a 3D segment is degenerate.
inline bool is_degenerate_3d(const Point& a, const Point& b)
{ return _cgal::seg3(a, b).is_degenerate(); }

/// Check if a 2D triangle is degenerate (collinear or duplicate vertices).
inline bool is_degenerate_2d(const Point& a, const Point& b, const Point& c)
{ return _cgal::tri2(a, b, c).is_degenerate(); }

/// Check if a 3D triangle is degenerate.
inline bool is_degenerate_3d(const Point& a, const Point& b, const Point& c)
{ return _cgal::tri(a, b, c).is_degenerate(); }

/// Check if a 3D tetrahedron is degenerate.
inline bool is_degenerate_3d(const Point& a, const Point& b,
                              const Point& c, const Point& d)
{ return _cgal::tet(a, b, c, d).is_degenerate(); }

// ---------------------------------------------------------------------------
// Intersection construction functions
// Uses CGAL::Exact_predicates_exact_constructions_kernel (EPECK) to ensure
// correct floating-point results when new points are constructed.
// ---------------------------------------------------------------------------

namespace _cgal_epeck {
  using EK   = CGAL::Exact_predicates_exact_constructions_kernel;
  using P3   = EK::Point_3;
  using S3   = EK::Segment_3;
  using T3   = EK::Triangle_3;
  using Tet3 = EK::Tetrahedron_3;
  using Poly = CGAL::Polyhedron_3<EK>;
  using Nef  = CGAL::Nef_polyhedron_3<EK>;

  inline P3   to3(const Point& p) { return P3(p.x(), p.y(), p.z()); }
  inline Point from3(const P3& p)
  {
    return Point(CGAL::to_double(p.x()),
                 CGAL::to_double(p.y()),
                 CGAL::to_double(p.z()));
  }
  inline S3   seg(const Point& a, const Point& b) { return S3(to3(a), to3(b)); }
  inline T3   tri(const Point& a, const Point& b, const Point& c)
  { return T3(to3(a), to3(b), to3(c)); }
}

/// Return the intersection of a 3D triangle and a 3D segment as a point set.
inline std::vector<Point>
cgal_intersection_triangle_segment_3d(const Point& p0, const Point& p1,
                                      const Point& p2,
                                      const Point& q0, const Point& q1)
{
  using namespace _cgal_epeck;
  const auto ii = CGAL::intersection(tri(p0, p1, p2), seg(q0, q1));
  if (!ii) return {};
  if (const P3* p = _cgalp_get_if<P3>(&*ii))
    return { from3(*p) };
  if (const S3* s = _cgalp_get_if<S3>(&*ii))
    return { from3(s->source()), from3(s->target()) };
  return {};
}

/// Return the vertices of the intersection of two 3D tetrahedra.
inline std::vector<Point>
cgal_intersection_tetrahedron_tetrahedron_3d(
    const Point& p0, const Point& p1, const Point& p2, const Point& p3,
    const Point& q0, const Point& q1, const Point& q2, const Point& q3)
{
  using namespace _cgal_epeck;
  Poly tet_a, tet_b;
  tet_a.make_tetrahedron(to3(p0), to3(p1), to3(p2), to3(p3));
  tet_b.make_tetrahedron(to3(q0), to3(q1), to3(q2), to3(q3));
  const Nef nef_a(tet_a), nef_b(tet_b);
  const Nef nef_i = nef_a * nef_b;
  std::vector<Point> res;
  for (auto v = nef_i.vertices_begin(); v != nef_i.vertices_end(); ++v)
    res.push_back(from3(v->point()));
  return res;
}

} // namespace dolfin

#endif // DOLFIN_WITH_CGAL
#endif // DOLFIN_CGAL_PREDICATES_H
