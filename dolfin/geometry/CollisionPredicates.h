// Copyright (C) 2014-2016 Anders Logg, August Johansson and Benjamin Kehlet
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
// Last changed: 2017-09-29

#ifndef __COLLISION_PREDICATES_H
#define __COLLISION_PREDICATES_H

#include <cstddef>

namespace simpex
{

  // Forward declarations
  class Point;
  /// This class implements algorithms for detecting pairwise
  /// collisions between mesh entities of varying dimensions.

  class CollisionPredicates
  {
  public:

    //--- Low-level collision detection predicates ---

    /// Check whether segment p0-p1 collides with point
    static bool collides_segment_point(const Point& p0,
                                       const Point& p1,
                                       const Point& point,
                                       std::size_t gdim);

    /// Check whether segment p0-p1 collides with point (1D version)
    static bool collides_segment_point_1d(double p0,
                                          double p1,
                                          double point);

    /// Check whether segment p0-p1 collides with point (2D version)
    static bool collides_segment_point_2d(const Point& p0,
                                          const Point& p1,
                                          const Point& point);

    /// Check whether segment p0-p1 collides with point (3D version)
    static bool collides_segment_point_3d(const Point& p0,
                                          const Point& p1,
                                          const Point& point);

    /// Check whether segment p0-p1 collides with segment q0-q1
    static bool collides_segment_segment(const Point& p0,
                                         const Point& p1,
                                         const Point& q0,
                                         const Point& q1,
                                         std::size_t gdim);

    /// Check whether segment p0-p1 collides with segment q0-q1 (1D version)
    static bool collides_segment_segment_1d(double p0,
                                            double p1,
                                            double q0,
                                            double q1);

    /// Check whether segment p0-p1 collides with segment q0-q1 (2D version)
    static bool collides_segment_segment_2d(const Point& p0,
                                            const Point& p1,
                                            const Point& q0,
                                            const Point& q1);

    /// Check whether segment p0-p1 collides with segment q0-q1 (3D version)
    static bool collides_segment_segment_3d(const Point& p0,
                                            const Point& p1,
                                            const Point& q0,
                                            const Point& q1);

    /// Check whether triangle p0-p1-p2 collides with point
    static bool collides_triangle_point(const Point& p0,
                                        const Point& p1,
                                        const Point& p2,
                                        const Point& point,
                                        std::size_t gdim);

    /// Check whether triangle p0-p1-p2 collides with point (2D version)
    static bool collides_triangle_point_2d(const Point& p0,
                                           const Point& p1,
                                           const Point& p2,
                                           const Point& point);

    /// Check whether triangle p0-p1-p2 collides with point (3D version)
    static bool collides_triangle_point_3d(const Point& p0,
                                           const Point& p1,
                                           const Point& p2,
                                           const Point& point);

    /// Check whether triangle p0-p1-p2 collides with segment q0-q1
    static bool collides_triangle_segment(const Point& p0,
                                          const Point& p1,
                                          const Point& p2,
                                          const Point& q0,
                                          const Point& q1,
                                          std::size_t gdim);

    /// Check whether triangle p0-p1-p2 collides with segment q0-q1 (2D version)
    static bool collides_triangle_segment_2d(const Point& p0,
                                             const Point& p1,
                                             const Point& p2,
                                             const Point& q0,
                                             const Point& q1);

    /// Check whether triangle p0-p1-p2 collides with segment q0-q1 (3D version)
    static bool collides_triangle_segment_3d(const Point& p0,
                                             const Point& p1,
                                             const Point& p2,
                                             const Point& q0,
                                             const Point& q1);

    /// Check whether triangle p0-p1-p2 collides with triangle q0-q1-q2
    static bool collides_triangle_triangle(const Point& p0,
                                           const Point& p1,
                                           const Point& p2,
                                           const Point& q0,
                                           const Point& q1,
                                           const Point& q2,
                                           std::size_t gdim);

    /// Check whether triangle p0-p1-p2 collides with triangle q0-q1-q2 (2D version)
    static bool collides_triangle_triangle_2d(const Point& p0,
                                              const Point& p1,
                                              const Point& p2,
                                              const Point& q0,
                                              const Point& q1,
                                              const Point& q2);

    /// Check whether triangle p0-p1-p2 collides with triangle q0-q1-q2 (3D version)
    static bool collides_triangle_triangle_3d(const Point& p0,
                                              const Point& p1,
                                              const Point& p2,
                                              const Point& q0,
                                              const Point& q1,
                                              const Point& q2);

    /// Check whether tetrahedron p0-p1-p2-p3 collides with point
    static bool collides_tetrahedron_point_3d(const Point& p0,
                                              const Point& p1,
                                              const Point& p2,
                                              const Point& p3,
                                              const Point& point);

    /// Check whether tetrahedron p0-p1-p2-p3 collides with segment q0-q1
    static bool collides_tetrahedron_segment_3d(const Point& p0,
                                                const Point& p1,
                                                const Point& p2,
                                                const Point& p3,
                                                const Point& q0,
                                                const Point& q1);

    /// Check whether tetrahedron p0-p1-p2-p3 collides with triangle q0-q1-q2
    static bool collides_tetrahedron_triangle_3d(const Point& p0,
                                                 const Point& p1,
                                                 const Point& p2,
                                                 const Point& p3,
                                                 const Point& q0,
                                                 const Point& q1,
                                                 const Point& q2);

    /// Check whether tetrahedron p0-p1-p2-p3 collides with tetrahedron q0-q1-q2
    static bool collides_tetrahedron_tetrahedron_3d(const Point& p0,
                                                    const Point& p1,
                                                    const Point& p2,
                                                    const Point& p3,
                                                    const Point& q0,
                                                    const Point& q1,
                                                    const Point& q2,
						    const Point& q3);

  private:

    // Implementation of collision detection predicates

    static bool _collides_segment_point_1d(double p0,
                                           double p1,
                                           double point);

    static bool _collides_segment_point_2d(const Point& p0,
                                           const Point& p1,
                                           const Point& point);

    static bool _collides_segment_point_3d(const Point& p0,
                                           const Point& p1,
                                           const Point& point);

    static bool _collides_segment_segment_1d(double p0,
                                             double p1,
                                             double q0,
                                             double q1);

    static bool _collides_segment_segment_2d(const Point& p0,
                                             const Point& p1,
                                             const Point& q0,
                                             const Point& q1);

    static bool _collides_segment_segment_3d(const Point& p0,
                                             const Point& p1,
                                             const Point& q0,
                                             const Point& q1);

    static bool _collides_triangle_point_2d(const Point& p0,
                                            const Point& p1,
                                            const Point& p2,
                                            const Point& point);

    static bool _collides_triangle_point_3d(const Point& p0,
                                            const Point& p1,
                                            const Point& p2,
                                            const Point& point);

    static bool _collides_triangle_segment_2d(const Point& p0,
                                              const Point& p1,
                                              const Point& p2,
                                              const Point& q0,
                                              const Point& q1);

    static bool _collides_triangle_segment_3d(const Point& p0,
                                              const Point& p1,
                                              const Point& p2,
                                              const Point& q0,
                                              const Point& q1);

    static bool _collides_triangle_triangle_2d(const Point& p0,
                                               const Point& p1,
                                               const Point& p2,
                                               const Point& q0,
                                               const Point& q1,
                                               const Point& q2);

    static bool _collides_triangle_triangle_3d(const Point& p0,
                                               const Point& p1,
                                               const Point& p2,
                                               const Point& q0,
                                               const Point& q1,
                                               const Point& q2);

    static bool _collides_tetrahedron_point_3d(const Point& p0,
                                               const Point& p1,
                                               const Point& p2,
                                               const Point& p3,
                                               const Point& point);

    static bool _collides_tetrahedron_segment_3d(const Point& p0,
                                                 const Point& p1,
                                                 const Point& p2,
                                                 const Point& p3,
                                                 const Point& q0,
                                                 const Point& q1);

    static bool _collides_tetrahedron_triangle_3d(const Point& p0,
                                                  const Point& p1,
                                                  const Point& p2,
                                                  const Point& p3,
                                                  const Point& q0,
                                                  const Point& q1,
                                                  const Point& q2);

    static bool _collides_tetrahedron_tetrahedron_3d(const Point& p0,
                                                     const Point& p1,
                                                     const Point& p2,
                                                     const Point& p3,
                                                     const Point& q0,
                                                     const Point& q1,
                                                     const Point& q2,
                                                     const Point& q3);

    // Helper: point-in-triangle assuming the point is already in the triangle's plane
    // (skips the orient3d coplanarity check).
    static bool _collides_triangle_point_3d_in_plane(const Point& p0,
                                                      const Point& p1,
                                                      const Point& p2,
                                                      const Point& point);

    // Helper: triangle-segment collision with pre-computed plane orientations
    // rsta = orient3d(r,s,t,a), rstb = orient3d(r,s,t,b)
    static bool _collides_triangle_segment_3d_with_hint(const Point& r,
                                                         const Point& s,
                                                         const Point& t,
                                                         const Point& a,
                                                         const Point& b,
                                                         double rsta,
                                                         double rstb);

    // Helper: 2D segment-point test with pre-computed orient2d value.
    // op01_pt = orient2d(p0, p1, point)
    static bool _collides_segment_point_2d_with_hint(const Point& p0,
                                                      const Point& p1,
                                                      const Point& point,
                                                      double op01_pt);

    // Helper: 2D segment-segment collision with all four orient2d values
    // pre-computed.
    // pq0 = orient2d(p0,p1,q0)  pq1 = orient2d(p0,p1,q1)
    // qp0 = orient2d(q0,q1,p0)  qp1 = orient2d(q0,q1,p1)
    static bool _collides_segment_segment_2d_with_hint(const Point& p0,
                                                        const Point& p1,
                                                        const Point& q0,
                                                        const Point& q1,
                                                        double pq0,
                                                        double pq1,
                                                        double qp0,
                                                        double qp1);

    // Helper: 2D triangle-point collision with pre-computed edge orientations.
    // ref   = orient2d(p0,p1,p2)
    // o01   = orient2d(p0,p1,point)
    // o12   = orient2d(p1,p2,point)
    // o20   = orient2d(p2,p0,point)
    static bool _collides_triangle_point_2d_with_hint(const Point& p0,
                                                       const Point& p1,
                                                       const Point& p2,
                                                       const Point& point,
                                                       double ref,
                                                       double o01,
                                                       double o12,
                                                       double o20);
  };

} // namespace simpex

// ============================================================================
// dolfin::CollisionPredicates
//
// Extends simpex::CollisionPredicates with high-level methods that accept
// dolfin MeshEntity objects.  Existing code using dolfin::CollisionPredicates
// therefore continues to work without modification.
// ============================================================================

namespace dolfin
{
  // Forward declaration (avoids circular #include with mesh headers).
  class MeshEntity;

  /// Collision predicates for dolfin mesh entities and low-level simplex geometry.
  class CollisionPredicates : public simpex::CollisionPredicates
  {
  public:
    //--- High-level collision detection predicates (dolfin mesh entities) ---

    /// Check whether entity collides with point.
    static bool collides(const MeshEntity& entity,
                         const simpex::Point& point);

    /// Check whether two entities collide.
    static bool collides(const MeshEntity& entity_0,
                         const MeshEntity& entity_1);
  };
}

#endif
