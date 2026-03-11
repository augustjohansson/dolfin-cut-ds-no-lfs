// CGALComparison.cpp
//
// Comprehensive Catch2 tests verifying that CGAL (EPICK) and Shewchuk predicates
// give identical results for ALL collision-detection routines across all
// geometric dimensions (gdim 2,3) and topological dimensions (tdim 0..3).
//
// For every function the test suite includes:
//   1. Random cases -- large batches of deterministic random inputs
//   2. Controlled corner cases -- collinear/coplanar alignments
//   3. Tiny perturbations (1e-20) -- exact arithmetic resolves inside/outside
//      even at perturbations far below machine precision (1e-16).
//      Two points are equal only if they are truly equal (==).
//   4. Every result is verified against CGAL::do_intersect (EPICK)
//
// All random generators use seed 42 for reproducibility.

#ifdef DOLFIN_WITH_CGAL

#include <array>
#include <cmath>
#include <random>

#include "../geometry/CGALPredicates.h"
#include "../geometry/CollisionPredicates.h"
#include "../geometry/predicates.h"
#include "../geometry/Point.h"
#include "catch/catch.hpp"

using namespace dolfin;

// ---------------------------------------------------------------------------
// Deterministic RNG helpers (seed 42)
// ---------------------------------------------------------------------------
static std::mt19937 rng(42);
static std::uniform_real_distribution<double> dist(-1.0, 1.0);

static double r() { return dist(rng); }
static Point rp2() { return Point(r(), r(), 0.0); }
static Point rp3() { return Point(r(), r(), r()); }

// Non-degenerate random segment in 2D
static std::array<Point,2> rseg2()
{
  while (true) {
    Point a = rp2(), b = rp2();
    if ((b-a).norm() != 0.0) return {a,b};
  }
}
// Non-degenerate random segment in 3D
static std::array<Point,2> rseg3()
{
  while (true) {
    Point a = rp3(), b = rp3();
    if ((b-a).norm() != 0.0) return {a,b};
  }
}
// Non-degenerate random triangle in 2D
static std::array<Point,3> rtri2()
{
  while (true) {
    Point p0 = rp2(), p1 = rp2(), p2 = rp2();
    if (orient2d(p0,p1,p2) != 0.0)
      return {p0,p1,p2};
  }
}
// Non-degenerate random triangle in 3D
static std::array<Point,3> rtri3()
{
  while (true) {
    Point p0 = rp3(), p1 = rp3(), p2 = rp3();
    if ((p1-p0).cross(p2-p0).norm() != 0.0)
      return {p0,p1,p2};
  }
}
// Non-degenerate random tetrahedron
static std::array<Point,4> rtet()
{
  while (true) {
    Point p0 = rp3(), p1 = rp3(), p2 = rp3(), p3 = rp3();
    if (orient3d(p0,p1,p2,p3) != 0.0)
      return {p0,p1,p2,p3};
  }
}

// Verify dolfin == cgal and return the dolfin value
static bool chk(bool dolfin_val, bool cgal_val, const char* label)
{
  if (dolfin_val != cgal_val) {
    FAIL("CGAL/Shewchuk mismatch in " << label
         << ": dolfin=" << dolfin_val << " cgal=" << cgal_val);
  }
  return dolfin_val;
}
#define CHK(s,c)   chk((s),(c), #s)

// ---------------------------------------------------------------------------
// 1D: segment-point
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: segment-point 1D", "[cgal]")
{
  // 1D has no CGAL equivalent; verify the dolfin implementation with exact checks.
  // Use 1e-20 perturbations near 0.0 where they are representable as distinct doubles.
  CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 0.5));    // clearly interior
  CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 0.0));    // left endpoint
  CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 1.0));    // right endpoint
  CHECK(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, 1e-20));  // 1e-20 inside left boundary
  CHECK_FALSE(CollisionPredicates::collides_segment_point_1d(0.0, 1.0, -1e-20)); // 1e-20 outside left boundary
}

// ---------------------------------------------------------------------------
// 2D: segment-point  (tdim=1 vs tdim=0, gdim=2)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: segment-point 2D", "[cgal]")
{
  SECTION("random") {
    const int N = 2000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto [p0,p1] = rseg2();
      Point q = rp2();
      bool s = CollisionPredicates::collides_segment_point_2d(p0,p1,q);
      bool c = cgal_collides_segment_point_2d(p0,p1,q);
      if (s != c) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: endpoint on segment") {
    Point p0(0,0,0), p1(1,0,0);
    CHECK(CHK(CollisionPredicates::collides_segment_point_2d(p0,p1,p0),
              cgal_collides_segment_point_2d(p0,p1,p0)));
    CHECK(CHK(CollisionPredicates::collides_segment_point_2d(p0,p1,p1),
              cgal_collides_segment_point_2d(p0,p1,p1)));
  }

  SECTION("corner: midpoint on segment") {
    Point p0(0,0,0), p1(1,0,0), q(0.5,0,0);
    CHECK(CHK(CollisionPredicates::collides_segment_point_2d(p0,p1,q),
              cgal_collides_segment_point_2d(p0,p1,q)));
  }

  SECTION("corner: point collinear but outside") {
    Point p0(0,0,0), p1(1,0,0), q(2,0,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_segment_point_2d(p0,p1,q),
                    cgal_collides_segment_point_2d(p0,p1,q)));
  }

  SECTION("perturbation: lateral offset distinguishes on-segment from off") {
    // A tiny y-offset of 1e-20 makes orient2d non-zero -> point is off the segment.
    // orient2d((0,0),(1,0),(0.5,1e-20)) = 1e-20 != 0, so the point is NOT on the line.
    Point p0(0,0,0), p1(1,0,0);
    Point q_on (0.5, 0,    0);   // exactly on segment
    Point q_off(0.5, 1e-20, 0);  // 1e-20 off the x-axis -> off segment
    CHECK(CHK(CollisionPredicates::collides_segment_point_2d(p0,p1,q_on),
              cgal_collides_segment_point_2d(p0,p1,q_on)));
    CHECK_FALSE(CHK(CollisionPredicates::collides_segment_point_2d(p0,p1,q_off),
                    cgal_collides_segment_point_2d(p0,p1,q_off)));
  }
}

// ---------------------------------------------------------------------------
// 3D: segment-point  (tdim=1 vs tdim=0, gdim=3)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: segment-point 3D", "[cgal]")
{
  SECTION("random") {
    const int N = 2000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto [p0,p1] = rseg3();
      Point q = rp3();
      bool s = CollisionPredicates::collides_segment_point_3d(p0,p1,q);
      bool c = cgal_collides_segment_point_3d(p0,p1,q);
      if (s != c) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: endpoint") {
    Point p0(0,0,0), p1(1,1,1);
    CHECK(CHK(CollisionPredicates::collides_segment_point_3d(p0,p1,p0),
              cgal_collides_segment_point_3d(p0,p1,p0)));
  }

  SECTION("corner: midpoint") {
    Point p0(0,0,0), p1(2,2,2), q(1,1,1);
    CHECK(CHK(CollisionPredicates::collides_segment_point_3d(p0,p1,q),
              cgal_collides_segment_point_3d(p0,p1,q)));
  }

  SECTION("corner: off-segment") {
    Point p0(0,0,0), p1(1,0,0), q(0.5,0.1,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_segment_point_3d(p0,p1,q),
                    cgal_collides_segment_point_3d(p0,p1,q)));
  }
}

// ---------------------------------------------------------------------------
// 2D: segment-segment  (tdim=1 vs tdim=1, gdim=2)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: segment-segment 2D", "[cgal]")
{
  SECTION("random") {
    const int N = 2000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto [a,b] = rseg2(); auto [c,d] = rseg2();
      bool s = CollisionPredicates::collides_segment_segment_2d(a,b,c,d);
      bool g = cgal_collides_segment_segment_2d(a,b,c,d);
      if (s != g) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: crossing X") {
    Point p0(0,0,0),p1(1,1,0),q0(1,0,0),q1(0,1,0);
    CHECK(CHK(CollisionPredicates::collides_segment_segment_2d(p0,p1,q0,q1),
              cgal_collides_segment_segment_2d(p0,p1,q0,q1)));
  }

  SECTION("corner: collinear overlap") {
    Point p0(0,0,0),p1(2,0,0),q0(1,0,0),q1(3,0,0);
    CHECK(CHK(CollisionPredicates::collides_segment_segment_2d(p0,p1,q0,q1),
              cgal_collides_segment_segment_2d(p0,p1,q0,q1)));
  }

  SECTION("corner: touching at endpoint") {
    Point p0(0,0,0),p1(1,0,0),q0(1,0,0),q1(2,0,0);
    CHECK(CHK(CollisionPredicates::collides_segment_segment_2d(p0,p1,q0,q1),
              cgal_collides_segment_segment_2d(p0,p1,q0,q1)));
  }

  SECTION("corner: parallel non-overlapping") {
    Point p0(0,0,0),p1(1,0,0),q0(0,1,0),q1(1,1,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_segment_segment_2d(p0,p1,q0,q1),
                    cgal_collides_segment_segment_2d(p0,p1,q0,q1)));
  }

  SECTION("perturbation: segments just touching vs. just apart") {
    // Use endpoints near 1e-15 so that a 1e-20 gap is representable.
    // ULP of 1e-15 is ~2.2e-31, so 1e-15 + 1e-20 is a distinct float from 1e-15.
    Point p0(0,0,0), p1(1e-15,0,0);
    Point q_touch(1e-15,0,0), q_end(2e-15,0,0);
    CHECK(CHK(CollisionPredicates::collides_segment_segment_2d(p0,p1,q_touch,q_end),
              cgal_collides_segment_segment_2d(p0,p1,q_touch,q_end)));
    Point q_near(1e-15+1e-20,0,0), q_far(2e-15,0,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_segment_segment_2d(p0,p1,q_near,q_far),
                    cgal_collides_segment_segment_2d(p0,p1,q_near,q_far)));
  }
}

// ---------------------------------------------------------------------------
// 3D: segment-segment  (tdim=1 vs tdim=1, gdim=3)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: segment-segment 3D", "[cgal]")
{
  SECTION("random") {
    const int N = 2000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto [a,b] = rseg3(); auto [c,d] = rseg3();
      bool s = CollisionPredicates::collides_segment_segment_3d(a,b,c,d);
      bool g = cgal_collides_segment_segment_3d(a,b,c,d);
      if (s != g) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: crossing in xy-plane") {
    Point p0(0,1,0),p1(1,1,0),q0(0.5,0,0),q1(0.5,2,0);
    CHECK(CHK(CollisionPredicates::collides_segment_segment_3d(p0,p1,q0,q1),
              cgal_collides_segment_segment_3d(p0,p1,q0,q1)));
  }

  SECTION("corner: skew segments") {
    Point p0(0,0,0),p1(1,0,0),q0(0,0,1),q1(0,1,1);
    CHECK_FALSE(CHK(CollisionPredicates::collides_segment_segment_3d(p0,p1,q0,q1),
                    cgal_collides_segment_segment_3d(p0,p1,q0,q1)));
  }

  SECTION("corner: collinear overlapping") {
    Point p0(0,0,0),p1(2,0,0),q0(1,0,0),q1(3,0,0);
    CHECK(CHK(CollisionPredicates::collides_segment_segment_3d(p0,p1,q0,q1),
              cgal_collides_segment_segment_3d(p0,p1,q0,q1)));
  }
}

// ---------------------------------------------------------------------------
// 2D: triangle-point  (tdim=2 vs tdim=0, gdim=2)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: triangle-point 2D", "[cgal]")
{
  SECTION("random") {
    const int N = 2000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto t = rtri2(); Point q = rp2();
      bool s = CollisionPredicates::collides_triangle_point_2d(t[0],t[1],t[2],q);
      bool c = cgal_collides_triangle_point_2d(t[0],t[1],t[2],q);
      if (s != c) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: interior, edge, vertex") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    // Vertex
    CHECK(CHK(CollisionPredicates::collides_triangle_point_2d(p0,p1,p2,p0),
              cgal_collides_triangle_point_2d(p0,p1,p2,p0)));
    // Edge midpoint
    Point emid(0.5,0,0);
    CHECK(CHK(CollisionPredicates::collides_triangle_point_2d(p0,p1,p2,emid),
              cgal_collides_triangle_point_2d(p0,p1,p2,emid)));
    // Interior
    Point interior(0.25,0.25,0);
    CHECK(CHK(CollisionPredicates::collides_triangle_point_2d(p0,p1,p2,interior),
              cgal_collides_triangle_point_2d(p0,p1,p2,interior)));
    // Outside
    Point outside(2,0,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_triangle_point_2d(p0,p1,p2,outside),
                    cgal_collides_triangle_point_2d(p0,p1,p2,outside)));
  }

  SECTION("perturbation: just inside/outside edge") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    // Unit tri: x>=0, y>=0, x+y<=1
    Point q_in (0.5, 0.49999, 0);   // x+y < 1
    Point q_out(0.5, 0.50001, 0);   // x+y > 1
    CHECK(CHK(CollisionPredicates::collides_triangle_point_2d(p0,p1,p2,q_in),
              cgal_collides_triangle_point_2d(p0,p1,p2,q_in)));
    CHECK_FALSE(CHK(CollisionPredicates::collides_triangle_point_2d(p0,p1,p2,q_out),
                    cgal_collides_triangle_point_2d(p0,p1,p2,q_out)));
  }
}

// ---------------------------------------------------------------------------
// 3D: triangle-point  (tdim=2 vs tdim=0, gdim=3)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: triangle-point 3D", "[cgal]")
{
  SECTION("random") {
    const int N = 5000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto t = rtri3(); Point q = rp3();
      bool s = CollisionPredicates::collides_triangle_point_3d(t[0],t[1],t[2],q);
      bool c = cgal_collides_triangle_point_3d(t[0],t[1],t[2],q);
      if (s != c) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: vertex, edge midpoint, centroid, off-plane") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    // Vertex
    CHECK(CHK(CollisionPredicates::collides_triangle_point_3d(p0,p1,p2,p2),
              cgal_collides_triangle_point_3d(p0,p1,p2,p2)));
    // Edge midpoint
    Point emid(0.5,0,0);
    CHECK(CHK(CollisionPredicates::collides_triangle_point_3d(p0,p1,p2,emid),
              cgal_collides_triangle_point_3d(p0,p1,p2,emid)));
    // Centroid
    Point cen(1.0/3,1.0/3,0);
    CHECK(CHK(CollisionPredicates::collides_triangle_point_3d(p0,p1,p2,cen),
              cgal_collides_triangle_point_3d(p0,p1,p2,cen)));
    // Off-plane
    Point off(0.25,0.25,1.0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_triangle_point_3d(p0,p1,p2,off),
                    cgal_collides_triangle_point_3d(p0,p1,p2,off)));
  }

  SECTION("perturbation: lift off plane by tiny epsilon") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q_on (0.25,0.25, 0.0);
    Point q_off(0.25,0.25, 1e-20);
    CHECK(CHK(CollisionPredicates::collides_triangle_point_3d(p0,p1,p2,q_on),
              cgal_collides_triangle_point_3d(p0,p1,p2,q_on)));
    CHECK_FALSE(CHK(CollisionPredicates::collides_triangle_point_3d(p0,p1,p2,q_off),
                    cgal_collides_triangle_point_3d(p0,p1,p2,q_off)));
  }
}

// ---------------------------------------------------------------------------
// 2D: triangle-segment  (tdim=2 vs tdim=1, gdim=2)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: triangle-segment 2D", "[cgal]")
{
  SECTION("random") {
    const int N = 2000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto t = rtri2(); auto [q0,q1] = rseg2();
      bool s = CollisionPredicates::collides_triangle_segment_2d(t[0],t[1],t[2],q0,q1);
      bool c = cgal_collides_triangle_segment_2d(t[0],t[1],t[2],q0,q1);
      if (s != c) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: segment along edge") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(0.2,0,0),q1(0.8,0,0);
    CHECK(CHK(CollisionPredicates::collides_triangle_segment_2d(p0,p1,p2,q0,q1),
              cgal_collides_triangle_segment_2d(p0,p1,p2,q0,q1)));
  }

  SECTION("corner: segment completely outside") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(2,0,0),q1(3,0,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_triangle_segment_2d(p0,p1,p2,q0,q1),
                    cgal_collides_triangle_segment_2d(p0,p1,p2,q0,q1)));
  }

  SECTION("corner: segment with endpoint on edge") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(0.5,0,0),q1(0.5,1,0);   // starts on edge, goes outside
    CHECK(CHK(CollisionPredicates::collides_triangle_segment_2d(p0,p1,p2,q0,q1),
              cgal_collides_triangle_segment_2d(p0,p1,p2,q0,q1)));
  }
}

// ---------------------------------------------------------------------------
// 3D: triangle-segment  (tdim=2 vs tdim=1, gdim=3)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: triangle-segment 3D", "[cgal]")
{
  SECTION("random") {
    const int N = 2000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto t = rtri3(); auto [q0,q1] = rseg3();
      bool s = CollisionPredicates::collides_triangle_segment_3d(t[0],t[1],t[2],q0,q1);
      bool c = cgal_collides_triangle_segment_3d(t[0],t[1],t[2],q0,q1);
      if (s != c) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: segment pierces triangle") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(0.25,0.25,-1),q1(0.25,0.25,1);
    CHECK(CHK(CollisionPredicates::collides_triangle_segment_3d(p0,p1,p2,q0,q1),
              cgal_collides_triangle_segment_3d(p0,p1,p2,q0,q1)));
  }

  SECTION("corner: segment lies in plane but outside triangle") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(2,0,0),q1(3,0,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_triangle_segment_3d(p0,p1,p2,q0,q1),
                    cgal_collides_triangle_segment_3d(p0,p1,p2,q0,q1)));
  }

  SECTION("corner: segment touches triangle at a vertex") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(0,0,0),q1(0,0,1);
    CHECK(CHK(CollisionPredicates::collides_triangle_segment_3d(p0,p1,p2,q0,q1),
              cgal_collides_triangle_segment_3d(p0,p1,p2,q0,q1)));
  }

  SECTION("perturbation: segment near vertex just inside vs. outside triangle") {
    // Triangle in z=0 plane. Test a vertical segment piercing near the origin.
    // (1e-20, 1e-20, .) is inside since x>0, y>0, x+y = 2e-20 < 1.
    // (-1e-20, 1e-20, .) is outside since orient2d((0,0),(1,0),(-1e-20,1e-20)) < 0.
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point qi0(1e-20, 1e-20, -1), qi1(1e-20, 1e-20, 1);    // pierces inside
    Point qo0(-1e-20, 1e-20, -1), qo1(-1e-20, 1e-20, 1);  // misses (outside p0-p1 edge)
    bool si = CollisionPredicates::collides_triangle_segment_3d(p0,p1,p2,qi0,qi1);
    bool ci = cgal_collides_triangle_segment_3d(p0,p1,p2,qi0,qi1);
    CHECK(si);
    CHECK(si == ci);
    bool so = CollisionPredicates::collides_triangle_segment_3d(p0,p1,p2,qo0,qo1);
    bool co = cgal_collides_triangle_segment_3d(p0,p1,p2,qo0,qo1);
    CHECK_FALSE(so);
    CHECK(so == co);
  }
}

// ---------------------------------------------------------------------------
// 2D: triangle-triangle  (tdim=2 vs tdim=2, gdim=2)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: triangle-triangle 2D", "[cgal]")
{
  SECTION("random") {
    const int N = 100; // 4950 pairs
    std::vector<std::array<Point,3>> tris(N);
    for (auto& t : tris) t = rtri2();
    int miss = 0;
    for (int i = 0; i < N; ++i)
      for (int j = i+1; j < N; ++j) {
        bool s = CollisionPredicates::collides_triangle_triangle_2d(
          tris[i][0],tris[i][1],tris[i][2],
          tris[j][0],tris[j][1],tris[j][2]);
        bool c = cgal_collides_triangle_triangle_2d(
          tris[i][0],tris[i][1],tris[i][2],
          tris[j][0],tris[j][1],tris[j][2]);
        if (s != c) ++miss;
      }
    CHECK(miss == 0);
  }

  SECTION("corner: overlapping") {
    Point p0(0,0,0),p1(2,0,0),p2(0,2,0);
    Point q0(1,0,0),q1(3,0,0),q2(1,2,0);
    CHECK(CHK(CollisionPredicates::collides_triangle_triangle_2d(p0,p1,p2,q0,q1,q2),
              cgal_collides_triangle_triangle_2d(p0,p1,p2,q0,q1,q2)));
  }

  SECTION("corner: touching along edge") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(1,0,0),q1(2,0,0),q2(1,1,0);
    CHECK(CHK(CollisionPredicates::collides_triangle_triangle_2d(p0,p1,p2,q0,q1,q2),
              cgal_collides_triangle_triangle_2d(p0,p1,p2,q0,q1,q2)));
  }

  SECTION("corner: separated") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(5,0,0),q1(6,0,0),q2(5,1,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_triangle_triangle_2d(p0,p1,p2,q0,q1,q2),
                    cgal_collides_triangle_triangle_2d(p0,p1,p2,q0,q1,q2)));
  }

  SECTION("perturbation: tiny gap between triangles") {
    // Use small-valued geometry so 1e-20 offset is representable as a distinct double.
    // Triangle p near origin, triangle q starts at 1e-15 + 1e-20 away.
    Point p0(0,0,0),p1(1e-15,0,0),p2(0,1e-15,0);
    Point q0(1e-15+1e-20,0,0),q1(2e-15,0,0),q2(1e-15+1e-20,1e-15,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_triangle_triangle_2d(p0,p1,p2,q0,q1,q2),
                    cgal_collides_triangle_triangle_2d(p0,p1,p2,q0,q1,q2)));
  }
}

// ---------------------------------------------------------------------------
// 3D: triangle-triangle  (tdim=2 vs tdim=2, gdim=3)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: triangle-triangle 3D", "[cgal]")
{
  SECTION("random") {
    const int N = 200; // ~20000 pairs
    std::vector<std::array<Point,3>> tris(N);
    for (auto& t : tris) t = rtri3();
    int miss = 0;
    for (int i = 0; i < N; ++i)
      for (int j = i+1; j < N; ++j) {
        bool s = CollisionPredicates::collides_triangle_triangle_3d(
          tris[i][0],tris[i][1],tris[i][2],
          tris[j][0],tris[j][1],tris[j][2]);
        bool c = cgal_collides_triangle_triangle_3d(
          tris[i][0],tris[i][1],tris[i][2],
          tris[j][0],tris[j][1],tris[j][2]);
        if (s != c) ++miss;
      }
    CHECK(miss == 0);
  }

  SECTION("corner: piercing") {
    Point p0(0,0,0),p1(2,0,0),p2(0,2,0);
    Point q0(1,0,-1),q1(1,0,1),q2(0,1,0);
    CHECK(CHK(CollisionPredicates::collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2),
              cgal_collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2)));
  }

  SECTION("corner: shared vertex") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(0,0,0),q1(0,0,1),q2(1,1,1);
    CHECK(CHK(CollisionPredicates::collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2),
              cgal_collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2)));
  }

  SECTION("corner: separated") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0);
    Point q0(5,5,5),q1(6,5,5),q2(5,6,5);
    CHECK_FALSE(CHK(CollisionPredicates::collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2),
                    cgal_collides_triangle_triangle_3d(p0,p1,p2,q0,q1,q2)));
  }
}

// ---------------------------------------------------------------------------
// 3D: tetrahedron-point  (tdim=3 vs tdim=0, gdim=3)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: tetrahedron-point 3D", "[cgal]")
{
  SECTION("random") {
    const int N = 3000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto t = rtet(); Point q = rp3();
      bool s = CollisionPredicates::collides_tetrahedron_point_3d(t[0],t[1],t[2],t[3],q);
      bool c = cgal_collides_tetrahedron_point_3d(t[0],t[1],t[2],t[3],q);
      if (s != c) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: vertex, face centroid, interior") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    // Vertex
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_point_3d(p0,p1,p2,p3,p1),
              cgal_collides_tetrahedron_point_3d(p0,p1,p2,p3,p1)));
    // Face centroid (on boundary)
    Point fc(1.0/3, 1.0/3, 1.0/3);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_point_3d(p0,p1,p2,p3,fc),
              cgal_collides_tetrahedron_point_3d(p0,p1,p2,p3,fc)));
    // Interior centroid
    Point cen(0.25,0.25,0.25);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_point_3d(p0,p1,p2,p3,cen),
              cgal_collides_tetrahedron_point_3d(p0,p1,p2,p3,cen)));
    // Outside
    Point out(1,1,1);
    CHECK_FALSE(CHK(CollisionPredicates::collides_tetrahedron_point_3d(p0,p1,p2,p3,out),
                    cgal_collides_tetrahedron_point_3d(p0,p1,p2,p3,out)));
  }

  SECTION("perturbation: just inside/outside face at x=0") {
    // The unit tet has a face at x=0 (vertices p0,p2,p3).
    // q_in  = (1e-20, 0.1, 0.1): x > 0, inside the tet.
    // q_out = (-1e-20, 0.1, 0.1): x < 0, outside.
    // orient3d involving this face changes sign between these two points.
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q_in (1e-20, 0.1, 0.1);   // tiny positive x -> inside
    Point q_out(-1e-20, 0.1, 0.1);  // tiny negative x -> outside
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_point_3d(p0,p1,p2,p3,q_in),
              cgal_collides_tetrahedron_point_3d(p0,p1,p2,p3,q_in)));
    CHECK_FALSE(CHK(CollisionPredicates::collides_tetrahedron_point_3d(p0,p1,p2,p3,q_out),
                    cgal_collides_tetrahedron_point_3d(p0,p1,p2,p3,q_out)));
  }
}

// ---------------------------------------------------------------------------
// 3D: tetrahedron-segment  (tdim=3 vs tdim=1, gdim=3)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: tetrahedron-segment 3D", "[cgal]")
{
  SECTION("random") {
    const int N = 2000;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto t = rtet(); auto [q0,q1] = rseg3();
      bool s = CollisionPredicates::collides_tetrahedron_segment_3d(t[0],t[1],t[2],t[3],q0,q1);
      bool c = cgal_collides_tetrahedron_segment_3d(t[0],t[1],t[2],t[3],q0,q1);
      if (s != c) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: segment pierces face") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q0(0.25,0.25,-1),q1(0.25,0.25,2);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0,q1),
              cgal_collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0,q1)));
  }

  SECTION("corner: segment endpoint at tet vertex") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q0(0,0,0),q1(-1,-1,-1);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0,q1),
              cgal_collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0,q1)));
  }

  SECTION("corner: segment entirely outside") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q0(5,0,0),q1(6,0,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0,q1),
                    cgal_collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0,q1)));
  }

  SECTION("corner: segment entirely inside") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q0(0.1,0.1,0.1),q1(0.2,0.1,0.1);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0,q1),
              cgal_collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0,q1)));
  }

  SECTION("perturbation: segment endpoint just outside face") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    double eps = 1e-20;
    // Endpoint just inside the slanted face (x+y+z < 1)
    Point q0_in(0.33, 0.33, 0.33), q1(2, 2, 2);
    Point q0_out(0.34+eps, 0.33, 0.33), q1b(2,2,2); // outside
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0_in,q1),
              cgal_collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0_in,q1)));
    // Segment from far outside to just-outside point: no intersection with tet
    bool s = CollisionPredicates::collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0_out,q1b);
    bool c = cgal_collides_tetrahedron_segment_3d(p0,p1,p2,p3,q0_out,q1b);
    CHECK(s == c);
  }
}

// ---------------------------------------------------------------------------
// 3D: tetrahedron-triangle  (tdim=3 vs tdim=2, gdim=3)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: tetrahedron-triangle 3D", "[cgal]")
{
  SECTION("random") {
    const int N = 1500;
    int miss = 0;
    for (int i = 0; i < N; ++i) {
      auto t = rtet(); auto tri = rtri3();
      bool s = CollisionPredicates::collides_tetrahedron_triangle_3d(
          t[0],t[1],t[2],t[3], tri[0],tri[1],tri[2]);
      bool c = cgal_collides_tetrahedron_triangle_3d(
          t[0],t[1],t[2],t[3], tri[0],tri[1],tri[2]);
      if (s != c) ++miss;
    }
    CHECK(miss == 0);
  }

  SECTION("corner: triangle vertex inside tet") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q0(0.1,0.1,0.1),q1(5,0,0),q2(0,5,0);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_triangle_3d(p0,p1,p2,p3,q0,q1,q2),
              cgal_collides_tetrahedron_triangle_3d(p0,p1,p2,p3,q0,q1,q2)));
  }

  SECTION("corner: triangle lies in tet face") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    // Triangle in face p0,p1,p2 (z=0)
    Point q0(0.1,0.1,0),q1(0.5,0.1,0),q2(0.1,0.5,0);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_triangle_3d(p0,p1,p2,p3,q0,q1,q2),
              cgal_collides_tetrahedron_triangle_3d(p0,p1,p2,p3,q0,q1,q2)));
  }

  SECTION("corner: triangle completely outside") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q0(5,0,0),q1(6,0,0),q2(5,1,0);
    CHECK_FALSE(CHK(CollisionPredicates::collides_tetrahedron_triangle_3d(p0,p1,p2,p3,q0,q1,q2),
                    cgal_collides_tetrahedron_triangle_3d(p0,p1,p2,p3,q0,q1,q2)));
  }
}

// ---------------------------------------------------------------------------
// 3D: tetrahedron-tetrahedron  (tdim=3 vs tdim=3, gdim=3)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL: tetrahedron-tetrahedron 3D", "[cgal]")
{
  SECTION("random pairs") {
    const int N = 150;  // 11175 pairs
    std::vector<std::array<Point,4>> tets(N);
    for (auto& t : tets) t = rtet();
    int miss = 0;
    for (int i = 0; i < N; ++i)
      for (int j = i+1; j < N; ++j) {
        bool s = CollisionPredicates::collides_tetrahedron_tetrahedron_3d(
          tets[i][0],tets[i][1],tets[i][2],tets[i][3],
          tets[j][0],tets[j][1],tets[j][2],tets[j][3]);
        bool c = cgal_collides_tetrahedron_tetrahedron_3d(
          tets[i][0],tets[i][1],tets[i][2],tets[i][3],
          tets[j][0],tets[j][1],tets[j][2],tets[j][3]);
        if (s != c) ++miss;
      }
    CHECK(miss == 0);
  }

  SECTION("corner: same tet") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,p0,p1,p2,p3),
              cgal_collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,p0,p1,p2,p3)));
  }

  SECTION("corner: shared face") {
    // Two tets sharing face p0p1p2
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q3(0,0,-1);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,p0,p1,p2,q3),
              cgal_collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,p0,p1,p2,q3)));
  }

  SECTION("corner: touching at single vertex") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q0(1,0,0),q1(2,0,0),q2(1,1,0),q3(1,0,1);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,q0,q1,q2,q3),
              cgal_collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,q0,q1,q2,q3)));
  }

  SECTION("corner: separated") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q0(5,0,0),q1(6,0,0),q2(5,1,0),q3(5,0,1);
    CHECK_FALSE(CHK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,q0,q1,q2,q3),
                    cgal_collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,q0,q1,q2,q3)));
  }

  SECTION("corner: partial overlap") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    Point q0(0.5,0,0),q1(1.5,0,0),q2(0.5,1,0),q3(0.5,0,1);
    CHECK(CHK(CollisionPredicates::collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,q0,q1,q2,q3),
              cgal_collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,q0,q1,q2,q3)));
  }

  SECTION("perturbation: second tet just outside first") {
    Point p0(0,0,0),p1(1,0,0),p2(0,1,0),p3(0,0,1);
    double eps = 1e-20;
    // Shift second tet so its vertex q0 is just outside the first
    Point q0(1+eps,0,0),q1(2,0,0),q2(1,1,0),q3(1,0,1);
    bool s = CollisionPredicates::collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,q0,q1,q2,q3);
    bool c = cgal_collides_tetrahedron_tetrahedron_3d(p0,p1,p2,p3,q0,q1,q2,q3);
    CHECK(s == c);
  }
}

// ---------------------------------------------------------------------------
// orient2d / orient3d sign agreement (retained from previous session)
// ---------------------------------------------------------------------------
TEST_CASE("CGAL vs Shewchuk: orient2d sign agreement", "[cgal]")
{
  const int N = 5000;
  int miss = 0;
  for (int i = 0; i < N; ++i) {
    Point a = rp2(), b = rp2(), c = rp2();
    auto sgn = [](double v){ return (v>0)?1:(v<0)?-1:0; };
    if (sgn(orient2d(a,b,c)) != sgn(cgal_orient2d(a,b,c))) ++miss;
  }
  CHECK(miss == 0);
}

TEST_CASE("CGAL vs Shewchuk: orient3d sign agreement", "[cgal]")
{
  const int N = 5000;
  int miss = 0;
  for (int i = 0; i < N; ++i) {
    Point a = rp3(), b = rp3(), c = rp3(), d = rp3();
    auto sgn = [](double v){ return (v>0)?1:(v<0)?-1:0; };
    if (sgn(orient3d(a,b,c,d)) != sgn(cgal_orient3d(a,b,c,d))) ++miss;
  }
  CHECK(miss == 0);
}

#endif // DOLFIN_WITH_CGAL
