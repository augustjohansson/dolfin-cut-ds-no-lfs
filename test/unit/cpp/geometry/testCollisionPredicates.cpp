// Unit tests for CollisionPredicates

#include "../geometry/CollisionPredicates.h"
#include "../geometry/Point.h"
#include "catch/catch.hpp"

using namespace dolfin;

// ---------------------------------------------------------------------------
// Segment-point collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: segment-point 1D")
{
  SECTION("point inside segment")
  {
    CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 0.5));
  }
  SECTION("point at left endpoint")
  {
    CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 0.0));
  }
  SECTION("point at right endpoint")
  {
    CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 1.0));
  }
  SECTION("point outside left")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, -0.5));
  }
  SECTION("point outside right")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 1.5));
  }
}

TEST_CASE("CollisionPredicates: segment-point 2D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);

  SECTION("point on segment")
  {
    CHECK(CollisionPredicates::collides_segment_point_2d(p0, p1, Point(0.5, 0.0, 0.0)));
  }
  SECTION("point at endpoint")
  {
    CHECK(CollisionPredicates::collides_segment_point_2d(p0, p1, p0));
    CHECK(CollisionPredicates::collides_segment_point_2d(p0, p1, p1));
  }
  SECTION("point off segment")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_2d(p0, p1, Point(0.5, 0.1, 0.0)));
  }
  SECTION("point collinear outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_2d(p0, p1, Point(2.0, 0.0, 0.0)));
  }
  SECTION("diagonal segment midpoint")
  {
    CHECK(CollisionPredicates::collides_segment_point_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 0.0), Point(0.5, 0.5, 0.0)));
  }
}

TEST_CASE("CollisionPredicates: segment-point 3D")
{
  SECTION("point on segment")
  {
    CHECK(CollisionPredicates::collides_segment_point_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0), Point(0.5, 0.5, 0.5)));
  }
  SECTION("point off segment")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_point_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.5, 0.0, 0.1)));
  }
  SECTION("point at endpoint")
  {
    const Point p0(0.0, 0.0, 0.0);
    const Point p1(1.0, 0.0, 0.0);
    CHECK(CollisionPredicates::collides_segment_point_3d(p0, p1, p0));
    CHECK(CollisionPredicates::collides_segment_point_3d(p0, p1, p1));
  }
}

// ---------------------------------------------------------------------------
// Segment-segment collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: segment-segment 1D")
{
  SECTION("overlapping")
  {
    CHECK(CollisionPredicates::collides_segment_segment_1d(0.0, 1.0, 0.5, 1.5));
  }
  SECTION("touching at endpoint")
  {
    CHECK(CollisionPredicates::collides_segment_segment_1d(0.0, 1.0, 1.0, 2.0));
  }
  SECTION("non-overlapping")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_segment_1d(0.0, 1.0, 2.5, 3.5));
  }
  SECTION("one inside other")
  {
    CHECK(CollisionPredicates::collides_segment_segment_1d(0.0, 2.0, 0.5, 1.5));
  }
  SECTION("identical")
  {
    CHECK(CollisionPredicates::collides_segment_segment_1d(0.0, 1.0, 0.0, 1.0));
  }
}

TEST_CASE("CollisionPredicates: segment-segment 2D")
{
  SECTION("crossing")
  {
    CHECK(CollisionPredicates::collides_segment_segment_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 0.0),
      Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0)));
  }
  SECTION("parallel non-overlapping")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_segment_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(1.0, 1.0, 0.0)));
  }
  SECTION("collinear overlapping")
  {
    CHECK(CollisionPredicates::collides_segment_segment_2d(
      Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0),
      Point(1.0, 0.0, 0.0), Point(3.0, 0.0, 0.0)));
  }
  SECTION("T-junction")
  {
    CHECK(CollisionPredicates::collides_segment_segment_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.5, 0.0, 0.0), Point(0.5, 1.0, 0.0)));
  }
}

TEST_CASE("CollisionPredicates: segment-segment 3D")
{
  SECTION("crossing in same plane")
  {
    CHECK(CollisionPredicates::collides_segment_segment_3d(
      Point(0.0, 0.5, 0.0), Point(1.0, 0.5, 0.0),
      Point(0.5, 0.0, 0.0), Point(0.5, 1.0, 0.0)));
  }
  SECTION("skew segments")
  {
    CHECK_FALSE(CollisionPredicates::collides_segment_segment_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 0.0, 1.0), Point(0.0, 1.0, 1.0)));
  }
}

// ---------------------------------------------------------------------------
// Triangle-point collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: triangle-point 2D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);

  SECTION("point inside")
  {
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, Point(0.25, 0.25, 0.0)));
  }
  SECTION("point on vertex")
  {
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, p0));
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, p1));
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, p2));
  }
  SECTION("point on edge")
  {
    CHECK(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, Point(0.5, 0.0, 0.0)));
  }
  SECTION("point outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, Point(1.0, 1.0, 0.0)));
  }
  SECTION("point far outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_point_2d(p0, p1, p2, Point(5.0, 5.0, 0.0)));
  }
}

TEST_CASE("CollisionPredicates: triangle-point 3D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);

  SECTION("point on triangle plane")
  {
    CHECK(CollisionPredicates::collides_triangle_point_3d(p0, p1, p2, Point(0.25, 0.25, 0.0)));
  }
  SECTION("point above triangle")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_point_3d(p0, p1, p2, Point(0.25, 0.25, 1.0)));
  }
}

// ---------------------------------------------------------------------------
// Triangle-triangle collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: triangle-triangle 2D")
{
  SECTION("overlapping triangles")
  {
    CHECK(CollisionPredicates::collides_triangle_triangle_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
      Point(0.25, 0.25, 0.0), Point(1.25, 0.25, 0.0), Point(0.25, 1.25, 0.0)));
  }
  SECTION("non-overlapping triangles")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_triangle_2d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
      Point(2.0, 0.0, 0.0), Point(3.0, 0.0, 0.0), Point(2.0, 1.0, 0.0)));
  }
}

TEST_CASE("CollisionPredicates: triangle-triangle 3D")
{
  SECTION("intersecting triangles")
  {
    CHECK(CollisionPredicates::collides_triangle_triangle_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
      Point(0.25, 0.25, -0.5), Point(0.25, 0.25, 0.5), Point(0.5, 0.5, 0.0)));
  }
}

// ---------------------------------------------------------------------------
// Tetrahedron-point collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: tetrahedron-point 3D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);
  const Point p3(0.0, 0.0, 1.0);

  SECTION("point inside tetrahedron")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_point_3d(
      p0, p1, p2, p3, Point(0.1, 0.1, 0.1)));
  }
  SECTION("point at vertex")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_point_3d(p0, p1, p2, p3, p0));
  }
  SECTION("point outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_point_3d(
      p0, p1, p2, p3, Point(1.0, 1.0, 1.0)));
  }
}

// ---------------------------------------------------------------------------
// Triangle-segment collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: triangle-segment 2D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);

  SECTION("segment crosses triangle")
  {
    CHECK(CollisionPredicates::collides_triangle_segment_2d(
      p0, p1, p2, Point(0.25, -0.5, 0.0), Point(0.25, 0.5, 0.0)));
  }
  SECTION("segment inside triangle")
  {
    CHECK(CollisionPredicates::collides_triangle_segment_2d(
      p0, p1, p2, Point(0.1, 0.1, 0.0), Point(0.2, 0.1, 0.0)));
  }
  SECTION("segment endpoint on triangle vertex")
  {
    CHECK(CollisionPredicates::collides_triangle_segment_2d(
      p0, p1, p2, Point(0.0, 0.0, 0.0), Point(-1.0, -1.0, 0.0)));
  }
  SECTION("segment entirely outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_segment_2d(
      p0, p1, p2, Point(2.0, 0.0, 0.0), Point(3.0, 0.0, 0.0)));
  }
  SECTION("segment parallel outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_segment_2d(
      p0, p1, p2, Point(0.0, 2.0, 0.0), Point(1.0, 2.0, 0.0)));
  }
  SECTION("segment on triangle edge")
  {
    CHECK(CollisionPredicates::collides_triangle_segment_2d(
      p0, p1, p2, Point(0.1, 0.0, 0.0), Point(0.9, 0.0, 0.0)));
  }
}

TEST_CASE("CollisionPredicates: triangle-segment 3D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);

  SECTION("segment pierces triangle")
  {
    CHECK(CollisionPredicates::collides_triangle_segment_3d(
      p0, p1, p2, Point(0.25, 0.25, -1.0), Point(0.25, 0.25, 1.0)));
  }
  SECTION("segment endpoint in triangle plane, inside")
  {
    CHECK(CollisionPredicates::collides_triangle_segment_3d(
      p0, p1, p2, Point(0.25, 0.25, 0.0), Point(0.25, 0.25, 1.0)));
  }
  SECTION("segment on same side of plane (no intersection)")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_segment_3d(
      p0, p1, p2, Point(0.25, 0.25, 1.0), Point(0.25, 0.25, 2.0)));
  }
  SECTION("segment pierces plane outside triangle")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_segment_3d(
      p0, p1, p2, Point(2.0, 2.0, -1.0), Point(2.0, 2.0, 1.0)));
  }
  SECTION("segment grazes triangle vertex")
  {
    CHECK(CollisionPredicates::collides_triangle_segment_3d(
      p0, p1, p2, Point(0.0, 0.0, -1.0), Point(0.0, 0.0, 1.0)));
  }
}

// ---------------------------------------------------------------------------
// Triangle-triangle 3D: additional cases
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: triangle-triangle 3D extra")
{
  SECTION("non-intersecting triangles (parallel planes)")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_triangle_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
      Point(0.0, 0.0, 1.0), Point(1.0, 0.0, 1.0), Point(0.0, 1.0, 1.0)));
  }
  SECTION("non-intersecting triangles (same plane, apart)")
  {
    CHECK_FALSE(CollisionPredicates::collides_triangle_triangle_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
      Point(3.0, 0.0, 0.0), Point(4.0, 0.0, 0.0), Point(3.0, 1.0, 0.0)));
  }
  SECTION("triangles sharing an edge")
  {
    CHECK(CollisionPredicates::collides_triangle_triangle_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0),
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, -1.0, 0.0)));
  }
}

// ---------------------------------------------------------------------------
// Tetrahedron-segment collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: tetrahedron-segment 3D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);
  const Point p3(0.0, 0.0, 1.0);

  SECTION("segment inside tetrahedron")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_segment_3d(
      p0, p1, p2, p3, Point(0.1, 0.1, 0.1), Point(0.2, 0.1, 0.1)));
  }
  SECTION("segment pierces tetrahedron face")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_segment_3d(
      p0, p1, p2, p3, Point(0.1, 0.1, -1.0), Point(0.1, 0.1, 0.5)));
  }
  SECTION("segment endpoint at tetrahedron vertex")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_segment_3d(
      p0, p1, p2, p3, Point(0.0, 0.0, 0.0), Point(-1.0, 0.0, 0.0)));
  }
  SECTION("segment entirely outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_segment_3d(
      p0, p1, p2, p3, Point(2.0, 0.0, 0.0), Point(3.0, 0.0, 0.0)));
  }
  SECTION("segment passes beside tetrahedron")
  {
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_segment_3d(
      p0, p1, p2, p3, Point(1.0, 1.0, 0.0), Point(1.0, 1.0, 1.0)));
  }
}

// ---------------------------------------------------------------------------
// Tetrahedron-triangle collisions
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: tetrahedron-triangle 3D")
{
  const Point p0(0.0, 0.0, 0.0);
  const Point p1(1.0, 0.0, 0.0);
  const Point p2(0.0, 1.0, 0.0);
  const Point p3(0.0, 0.0, 1.0);

  SECTION("triangle inside tetrahedron")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_triangle_3d(
      p0, p1, p2, p3,
      Point(0.1, 0.1, 0.1), Point(0.2, 0.1, 0.1), Point(0.1, 0.2, 0.1)));
  }
  SECTION("triangle cuts through tetrahedron face")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_triangle_3d(
      p0, p1, p2, p3,
      Point(0.2, 0.2, -0.5), Point(0.2, 0.2, 0.5), Point(0.4, 0.1, 0.0)));
  }
  SECTION("triangle entirely outside")
  {
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_triangle_3d(
      p0, p1, p2, p3,
      Point(2.0, 0.0, 0.0), Point(3.0, 0.0, 0.0), Point(2.0, 1.0, 0.0)));
  }
  SECTION("triangle on tetrahedron face (shared face)")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_triangle_3d(
      p0, p1, p2, p3, p0, p1, p2));
  }
}

// ---------------------------------------------------------------------------
// Near-degenerate / near-touching cases
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: near-touching (1e-14 separation)")
{
  SECTION("tet-point just outside (1e-14)")
  {
    const Point p0(0.0, 0.0, 0.0);
    const Point p1(1.0, 0.0, 0.0);
    const Point p2(0.0, 1.0, 0.0);
    const Point p3(0.0, 0.0, 1.0);
    // Point just outside the x=0 face
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_point_3d(
      p0, p1, p2, p3, Point(-1e-14, 0.0, 0.0)));
  }
  SECTION("tri-seg 3D just missing (segment does not reach plane)")
  {
    const Point p0(0.0, 0.0, 0.0);
    const Point p1(1.0, 0.0, 0.0);
    const Point p2(0.0, 1.0, 0.0);
    // Segment above the triangle plane, does not cross it
    CHECK_FALSE(CollisionPredicates::collides_triangle_segment_3d(
      p0, p1, p2,
      Point(0.25, 0.25, 1e-14), Point(0.25, 0.25, 1.0)));
  }
}

TEST_CASE("CollisionPredicates: tetrahedron-tetrahedron 3D")
{
  SECTION("overlapping tetrahedra")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0),
      Point(0.1, 0.1, 0.1), Point(1.1, 0.1, 0.1),
      Point(0.1, 1.1, 0.1), Point(0.1, 0.1, 1.1)));
  }
  SECTION("non-overlapping tetrahedra")
  {
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0),
      Point(3.0, 0.0, 0.0), Point(4.0, 0.0, 0.0),
      Point(3.0, 1.0, 0.0), Point(3.0, 0.0, 1.0)));
  }
}


// ---------------------------------------------------------------------------
// Magnitude tests: entities of size ~1e-13 and ~1e13
// ---------------------------------------------------------------------------

TEST_CASE("CollisionPredicates: magnitude tests ~1e-13")
{
  const double s = 1e-13;

  SECTION("tet-point inside (small scale)")
  {
    // Unit tet scaled by s
    CHECK(CollisionPredicates::collides_tetrahedron_point_3d(
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0),
      Point(0.0, s, 0.0),   Point(0.0, 0.0, s),
      Point(s*0.1, s*0.1, s*0.1)));
  }

  SECTION("tet-point outside (small scale)")
  {
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_point_3d(
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0),
      Point(0.0, s, 0.0),   Point(0.0, 0.0, s),
      Point(s*2.0, 0.0, 0.0)));
  }

  SECTION("tri-tri intersecting (small scale)")
  {
    CHECK(CollisionPredicates::collides_triangle_triangle_3d(
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0), Point(0.0, s, 0.0),
      Point(s*0.25, s*0.25, -s), Point(s*0.25, s*0.25, s), Point(s*0.3, s*0.1, 0.0)));
  }

  SECTION("tet-tet intersecting (small scale)")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0),
      Point(0.0, s, 0.0),   Point(0.0, 0.0, s),
      Point(s*0.1, s*0.1, s*0.1), Point(s*1.1, s*0.1, s*0.1),
      Point(s*0.1, s*1.1, s*0.1), Point(s*0.1, s*0.1, s*1.1)));
  }
}

TEST_CASE("CollisionPredicates: magnitude tests ~1e13")
{
  const double s = 1e13;

  SECTION("tet-point inside (large scale)")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_point_3d(
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0),
      Point(0.0, s, 0.0),   Point(0.0, 0.0, s),
      Point(s*0.1, s*0.1, s*0.1)));
  }

  SECTION("tet-point outside (large scale)")
  {
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_point_3d(
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0),
      Point(0.0, s, 0.0),   Point(0.0, 0.0, s),
      Point(s*2.0, 0.0, 0.0)));
  }

  SECTION("tri-tri intersecting (large scale)")
  {
    CHECK(CollisionPredicates::collides_triangle_triangle_3d(
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0), Point(0.0, s, 0.0),
      Point(s*0.25, s*0.25, -s), Point(s*0.25, s*0.25, s), Point(s*0.3, s*0.1, 0.0)));
  }

  SECTION("tet-tet intersecting (large scale)")
  {
    CHECK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0),
      Point(0.0, s, 0.0),   Point(0.0, 0.0, s),
      Point(s*0.1, s*0.1, s*0.1), Point(s*1.1, s*0.1, s*0.1),
      Point(s*0.1, s*1.1, s*0.1), Point(s*0.1, s*0.1, s*1.1)));
  }
}

TEST_CASE("CollisionPredicates: mixed magnitude (1e-13 vs 1e13)")
{
  SECTION("seg-point: segment endpoint at 1e-13, point at 1e13 scale")
  {
    // A tiny segment at origin and a far-away point
    CHECK_FALSE(CollisionPredicates::collides_segment_point_3d(
      Point(0.0, 0.0, 0.0), Point(1e-13, 0.0, 0.0),
      Point(1e13, 0.0, 0.0)));
  }

  SECTION("tet not colliding at wildly different scales")
  {
    // Tiny tet near origin, point at 1e13
    CHECK_FALSE(CollisionPredicates::collides_tetrahedron_point_3d(
      Point(0.0, 0.0, 0.0), Point(1e-13, 0.0, 0.0),
      Point(0.0, 1e-13, 0.0), Point(0.0, 0.0, 1e-13),
      Point(1e13, 1e13, 1e13)));
  }
}
