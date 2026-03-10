// Copyright (C) 2016-2017 Anders Logg, August Johansson and Benjamin Kehlet
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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <vector>
#include "CGALExactArithmetic.h"
#include "CollisionPredicates.h"
#include "ConvexTriangulation.h"
#include "GeometryPredicates.h"
#include "GeometryTools.h"
#include "IntersectionConstruction.h"
#include "predicates.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
namespace
{
  // Create a unique list of points in the sense that |p-q| > tol in each dimension
  std::vector<Point>
  unique_points(const std::vector<Point>& input_points,
		std::size_t gdim,
		double tol)
  {
    std::vector<Point> points;

    for (std::size_t i = 0; i < input_points.size(); ++i)
    {
      bool unique = true;
      for (std::size_t j = i+1; j < input_points.size(); ++j)
      {
	std::size_t cnt = 0;
	for (std::size_t d = 0; d < gdim; ++d)
	{
	  if (std::abs(input_points[i][d] - input_points[j][d]) > tol)
	  {
	    cnt++;
	  }
	}
	if (cnt == 0)
	{
	  unique = false;
	  break;
	}
      }

      if (unique)
	points.push_back(input_points[i]);
    }

    return points;
  }


struct ProjPt
{
  double u, v;
  std::size_t idx;
};

static inline int dominant_drop_axis_from_normal(const Point& n) noexcept
{
  const double ax = std::abs(n.x());
  const double ay = std::abs(n.y());
  const double az = std::abs(n.z());
  if (ax >= ay && ax >= az) return 0; // drop x -> use (y,z)
  if (ay >= ax && ay >= az) return 1; // drop y -> use (x,z)
  return 2;                           // drop z -> use (x,y)
}

static inline void project_drop_axis(const Point& p, const Point& origin,
                                     int drop, double& u, double& v) noexcept
{
  const double x = p.x() - origin.x();
  const double y = p.y() - origin.y();
  const double z = p.z() - origin.z();
  if (drop == 0)      { u = y; v = z; }
  else if (drop == 1) { u = x; v = z; }
  else                { u = x; v = y; }
}

// Find a “good” normal direction by selecting the max-area triangle.
// Returns false if all points are collinear/identical (no plane).
static bool find_stable_normal(const std::vector<Point>& P, Point& N_out,
                               std::size_t& anchor_out)
{
  const std::size_t n = P.size();
  double best = 0.0;
  Point bestN;
  std::size_t bestA = 0;

  for (std::size_t a = 0; a < n; ++a)
    for (std::size_t b = a + 1; b < n; ++b)
      for (std::size_t c = b + 1; c < n; ++c)
      {
        Point N = GeometryTools::cross_product(P[a], P[b], P[c]); // (b-a)x(c-a)
        const double s = N.squared_norm();
        if (s > best)
        {
          best = s;
          bestN = N;
          bestA = a;
        }
      }

  if (best == 0.0)
    return false;

  N_out = bestN;     // NOTE: not normalized
  anchor_out = bestA;
  return true;
}

// Hybrid turn: orient2d fast; if 0, use orient3d(a,b,a+N,c) fallback.
static inline double turn_sign_hybrid(const ProjPt& A, const ProjPt& B, const ProjPt& C,
                                      const std::vector<Point>& P,
                                      const Point& N)
{
  const double a2[2] = {A.u, A.v};
  const double b2[2] = {B.u, B.v};
  const double c2[2] = {C.u, C.v};
  const double o2 = _orient2d(a2, b2, c2);
  if (o2 != 0.0) return o2;

  // Fallback: robust in-plane orientation using 3D predicate
  const Point& a3 = P[A.idx];
  const Point& b3 = P[B.idx];
  const Point& c3 = P[C.idx];
  const Point r = a3 + N;
  return orient3d(a3, b3, r, c3);
}

static std::vector<std::size_t>
convex_hull_planar_indices_robust(const std::vector<Point>& points)
{
  const std::size_t n = points.size();
  if (n == 0) return {};
  if (n == 1) return {0};

  Point N;
  std::size_t anchor = 0;
  if (!find_stable_normal(points, N, anchor))
  {
    // 1D hull: choose endpoints along axis with largest spread
    double minx = points[0].x(), maxx = points[0].x();
    double miny = points[0].y(), maxy = points[0].y();
    double minz = points[0].z(), maxz = points[0].z();
    for (std::size_t i = 1; i < n; ++i)
    {
      minx = std::min(minx, points[i].x()); maxx = std::max(maxx, points[i].x());
      miny = std::min(miny, points[i].y()); maxy = std::max(maxy, points[i].y());
      minz = std::min(minz, points[i].z()); maxz = std::max(maxz, points[i].z());
    }
    const double sx = maxx - minx, sy = maxy - miny, sz = maxz - minz;
    int axis = (sx >= sy && sx >= sz) ? 0 : (sy >= sx && sy >= sz) ? 1 : 2;

    auto coord = [&](const Point& p) -> double {
      return (axis == 0) ? p.x() : (axis == 1) ? p.y() : p.z();
    };

    std::size_t lo = 0, hi = 0;
    for (std::size_t i = 1; i < n; ++i)
    {
      if (coord(points[i]) < coord(points[lo])) lo = i;
      if (coord(points[i]) > coord(points[hi])) hi = i;
    }
    if (lo == hi) return {lo};
    return {lo, hi};
  }

  const int drop = dominant_drop_axis_from_normal(N);
  const Point origin = points[anchor];

  // Project
  std::vector<ProjPt> P2;
  P2.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    double u, v;
    project_drop_axis(points[i], origin, drop, u, v);
    P2.push_back({u, v, i});
  }

  // Sort by (u,v)
  std::sort(P2.begin(), P2.end(), [](const ProjPt& a, const ProjPt& b) {
    if (a.u < b.u) return true;
    if (a.u > b.u) return false;
    if (a.v < b.v) return true;
    if (a.v > b.v) return false;
    return a.idx < b.idx;
  });

  // Remove exact duplicates in projection
  P2.erase(std::unique(P2.begin(), P2.end(), [](const ProjPt& a, const ProjPt& b) {
    return a.u == b.u && a.v == b.v;
  }), P2.end());

  if (P2.size() == 1) return {P2[0].idx};
  if (P2.size() == 2) return {P2[0].idx, P2[1].idx};

  // Andrew monotone chain, pop on <=0 to discard collinear boundary points
  std::vector<ProjPt> lower, upper;
  lower.reserve(P2.size());
  upper.reserve(P2.size());

  for (const auto& p : P2)
  {
    while (lower.size() >= 2)
    {
      const auto& A = lower[lower.size() - 2];
      const auto& B = lower[lower.size() - 1];
      const double o = turn_sign_hybrid(A, B, p, points, N);
      if (o <= 0.0) lower.pop_back();
      else break;
    }
    lower.push_back(p);
  }

  for (std::size_t k = P2.size(); k-- > 0; )
  {
    const auto& p = P2[k];
    while (upper.size() >= 2)
    {
      const auto& A = upper[upper.size() - 2];
      const auto& B = upper[upper.size() - 1];
      const double o = turn_sign_hybrid(A, B, p, points, N);
      if (o <= 0.0) upper.pop_back();
      else break;
    }
    upper.push_back(p);
  }

  lower.pop_back();
  upper.pop_back();

  std::vector<std::size_t> hull;
  hull.reserve(lower.size() + upper.size());
  for (const auto& p : lower) hull.push_back(p.idx);
  for (const auto& p : upper) hull.push_back(p.idx);

  return hull;
}

std::vector<std::pair<std::size_t, std::size_t>>
compute_convex_hull_planar(const std::vector<Point>& points)
{
  std::vector<std::pair<std::size_t, std::size_t>> edges;

  const auto hull = convex_hull_planar_indices_robust(points);
  if (hull.size() < 2) return edges;

  edges.reserve(hull.size());
  for (std::size_t i = 0; i < hull.size(); ++i)
  {
    const std::size_t a = hull[i];
    const std::size_t b = hull[(i + 1) % hull.size()];
    edges.emplace_back(a, b);
  }

  return edges;
}

} // anonymous namespace

//------------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::triangulate(const std::vector<Point>& p,
                                 std::size_t gdim,
                                 std::size_t tdim)
{
  if (p.empty())
    return std::vector<std::vector<Point>>();

  if (tdim == 1)
  {
    return triangulate_1d(p, gdim);
  }
  else if (tdim == 2 && gdim == 2)
  {
    return triangulate_graham_scan_2d(p);
  }
  else if (tdim == 3 && gdim == 3)
  {
    return triangulate_graham_scan_3d(p);
  }

  throw std::runtime_error("Triangulation of polyhedron of topological dimension %u and geometric dimension %u not implemented");

  return std::vector<std::vector<Point>>();
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_1d(const std::vector<Point>& p,
				     std::size_t gdim)
{
  // A convex polyhedron of topological dimension 1 can not have more
  // than two points. If more, they must be collinear (more or
  // less). This can happen due to tolerances in
  // IntersectionConstruction::intersection_segment_segment_2d.

  if (gdim != 2)
  {
    throw std::runtime_error("Function is only implemented for gdim = 2");
  }

  const std::vector<Point> unique_p = unique_points(p, gdim, 3.0e-16);

  if (unique_p.size() == 2)
  {
    std::vector<std::vector<Point>> t { unique_p };
    return t;
  }
  else if (unique_p.size() < 2)
  {
    return std::vector<std::vector<Point>>();
  }
  else
  {
    // Here unique_p.size() > 2: points must be approximately collinear
    bool collinear = true;
    for (std::size_t i = 2; i < unique_p.size(); ++i)
    {
      const double o = orient2d(unique_p[0], unique_p[1], unique_p[i]);
      if (std::abs(o) > 1e-14)
      {
	collinear = false;
	break;
      }
    }

    assert(collinear);

    // Return extremal points
    const Point v = unique_p[1] - unique_p[0];
    std::vector<std::pair<double, std::size_t>> order;
    order.emplace_back(0.0, 0);
    order.emplace_back(1.0, 1);
    for (std::size_t i = 2; i < unique_p.size(); ++i)
      order.emplace_back(v.dot(unique_p[i]-unique_p[0]), i);

    std::sort(order.begin(), order.end());

    return {{ unique_p[order.front().second],
	  unique_p[order.back().second] }};
  }

}

//------------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_graham_scan_2d(const std::vector<Point>& input_points)
{
  assert(GeometryPredicates::is_finite(input_points));

  constexpr std::size_t tdim = 2;
  constexpr std::size_t gdim = 2;

  std::vector<Point> points = unique_points(input_points, gdim, 3.0e-16);
  const std::size_t n = points.size();

  std::vector<std::vector<Point>> triangulation;
  if (n < 3) return triangulation;

  if (n == 3)
  {
    if (!GeometryPredicates::is_degenerate(points, tdim, gdim))
      triangulation.push_back(points);
    return triangulation;
  }

  // 1) Choose pivot: lowest y, break ties by lowest x
  std::size_t pivot = 0;
  for (std::size_t i = 1; i < n; ++i)
  {
    if (points[i].y() < points[pivot].y() ||
        (points[i].y() == points[pivot].y() && points[i].x() < points[pivot].x()))
      pivot = i;
  }
  std::swap(points[0], points[pivot]);
  const Point& p0 = points[0];

  // 2) Create index list for points[1..n-1]
  std::vector<std::size_t> idx;
  idx.reserve(n - 1);
  for (std::size_t i = 1; i < n; ++i) idx.push_back(i);

  auto dist2 = [&](std::size_t i) -> double
  {
    const double dx = points[i].x() - p0.x();
    const double dy = points[i].y() - p0.y();
    return dx*dx + dy*dy;
  };

  // 3) Polar sort by orientation around p0
  std::sort(idx.begin(), idx.end(),
            [&](std::size_t ia, std::size_t ib)
  {
    const Point& a = points[ia];
    const Point& b = points[ib];

    const double o = orient2d(p0, a, b);
    if (o > 0.0) return true;   // a before b (CCW)
    if (o < 0.0) return false;

    // Collinear with pivot: sort by increasing distance (near first)
    // This makes it easy to drop interior collinear points later.
    return dist2(ia) < dist2(ib);
  });

  // 4) Remove points that are collinear on the same ray from pivot
  // Keep only the farthest point on each ray to avoid degenerate fan triangles.
  std::vector<std::size_t> order;
  order.reserve(idx.size());
  for (std::size_t t = 0; t < idx.size(); )
  {
    std::size_t best = idx[t]; // farthest at the end of a collinear block (since dist ascending)
    std::size_t u = t + 1;
    while (u < idx.size() && orient2d(p0, points[idx[t]], points[idx[u]]) == 0.0)
    {
      best = idx[u];
      ++u;
    }
    order.push_back(best);
    t = u;
  }

  if (order.size() < 2)
    return triangulation; // all points collinear w.r.t pivot => no area

  // 5) Fan triangulation around p0
  triangulation.reserve(order.size() - 1);

  for (std::size_t i = 0; i + 1 < order.size(); ++i)
  {
    const Point& a = points[order[i]];
    const Point& b = points[order[i + 1]];

    // If you trust the filtering above, this orient2d check is usually enough:
    // if (orient2d(p0, a, b) == 0.0) continue;

    std::vector<Point> tri;
    tri.reserve(3);
    tri.push_back(p0);
    tri.push_back(a);
    tri.push_back(b);

    if (!GeometryPredicates::is_degenerate(tri, tdim, gdim))
      triangulation.push_back(std::move(tri));
  }

  return triangulation;
}
//-----------------------------------------------------------------------------
static inline int signum(double x) noexcept
{
  return (x > 0.0) - (x < 0.0);
}

// Pack an ordered triple (i<j<k) into a 64-bit key.
// Works as long as indices fit in 21 bits each (2 million). For n~100 you're fine.
static inline std::uint64_t pack3(std::uint32_t i, std::uint32_t j, std::uint32_t k) noexcept
{
  return (std::uint64_t(i) << 42) | (std::uint64_t(j) << 21) | std::uint64_t(k);
}

std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_graham_scan_3d(const std::vector<Point>& input_points)
{
  assert(GeometryPredicates::is_finite(input_points));

  constexpr std::size_t tdim = 3;
  constexpr std::size_t gdim = 3;

  std::vector<Point> points = unique_points(input_points, gdim, 3.0e-16);

  std::vector<std::vector<Point>> triangulation;
  triangulation.reserve(points.size()); // heuristic

  const std::size_t n = points.size();
  if (n < 4)
    return triangulation;

  if (n == 4)
  {
    if (!GeometryPredicates::is_degenerate(points, tdim, gdim))
      triangulation.push_back(points);
    return triangulation;
  }

  // Polyhedron center
  Point polyhedroncenter(0, 0, 0);
  for (const Point& p : points) polyhedroncenter += p;
  polyhedroncenter /= static_cast<double>(n);

  // Much faster than std::set<tuple<...>>
  std::unordered_set<std::uint64_t> checked;
  checked.reserve(n * n * 4); // heuristic, avoids rehash churn

  // Helper to push a tet candidate without repeated small allocations
  auto push_tet = [&](const Point& a, const Point& b, const Point& c, const Point& d)
  {
    std::vector<Point> cand;
    cand.reserve(4);
    cand.push_back(a);
    cand.push_back(b);
    cand.push_back(c);
    cand.push_back(d);

#ifdef DOLFIN_ENABLE_GEOMETRY_DEBUGGING
    if (cgal_tet_is_degenerate(cand))
      throw std::runtime_error("tet is degenerate");
#endif

    if (!GeometryPredicates::is_degenerate(cand, tdim, gdim))
      triangulation.push_back(std::move(cand));
  };

  for (std::size_t i = 0; i < n; ++i)
  {
    for (std::size_t j = i + 1; j < n; ++j)
    {
      for (std::size_t k = j + 1; k < n; ++k)
      {
        const std::uint64_t key = pack3((std::uint32_t)i, (std::uint32_t)j, (std::uint32_t)k);
        if (!checked.emplace(key).second)
          continue;

        // Skip collinear triples
        if ((points[j] - points[i]).cross(points[k] - points[i]).squared_norm() == 0.0)
          continue;

        bool on_convex_hull = true;
        int ref_sign = 0;

        // Collect coplanar indices only if we actually encounter any.
        // Start with i,j,k in the list.
        std::vector<std::size_t> coplanar;
        coplanar.reserve(16);
        coplanar.push_back(i);
        coplanar.push_back(j);
        coplanar.push_back(k);

        for (std::size_t m = 0; m < n; ++m)
        {
          if (m == i || m == j || m == k) continue;

          const double o = orient3d(points[i], points[j], points[k], points[m]);

          if (o == 0.0)
          {
            coplanar.push_back(m);
            continue;
          }

          const int s = signum(o);
          if (ref_sign == 0)
          {
            ref_sign = s;
          }
          else if (s != ref_sign)
          {
            on_convex_hull = false;
            break; // stop calling orient3d once we know it fails
          }
        }

        if (!on_convex_hull)
          continue;

        if (coplanar.size() == 3)
        {
          // Simple face
          push_tet(points[i], points[j], points[k], polyhedroncenter);
          continue;
        }

        // ---- coplanar.size() > 3 ----
        // Build coplanar points (local indexing)
        std::vector<Point> coplanar_points;
        coplanar_points.reserve(coplanar.size());
        for (std::size_t idx : coplanar)
          coplanar_points.push_back(points[idx]);

        // Compute planar convex hull edges among coplanar points (local indices)
        const auto coplanar_convex_hull = compute_convex_hull_planar(coplanar_points);

        // Center of the coplanar set
        Point coplanar_center(0, 0, 0);
        for (const Point& p : coplanar_points) coplanar_center += p;
        coplanar_center /= static_cast<double>(coplanar_points.size());

        // For each hull edge, create a tet fan
        for (const auto& edge : coplanar_convex_hull)
        {
          const Point& a = coplanar_points[edge.first];
          const Point& b = coplanar_points[edge.second];

          // order matches your original: {polyhedroncenter, coplanar_center, a, b}
          push_tet(polyhedroncenter, coplanar_center, a, b);

#ifdef DOLFIN_ENABLE_GEOMETRY_DEBUGGING
          if (!triangulation.empty() && cgal_triangulation_overlap(triangulation))
            throw std::runtime_error("now triangulation overlaps");
#endif
        }

        // Only mark coplanar triples once
        std::sort(coplanar.begin(), coplanar.end());
        const std::size_t csz = coplanar.size();
        for (std::size_t a = 0; a + 2 < csz; ++a)
        {
          const std::uint32_t ia = (std::uint32_t)coplanar[a];
          for (std::size_t b = a + 1; b + 1 < csz; ++b)
          {
            const std::uint32_t ib = (std::uint32_t)coplanar[b];
            for (std::size_t c = b + 1; c < csz; ++c)
            {
              const std::uint32_t ic = (std::uint32_t)coplanar[c];
              checked.emplace(pack3(ia, ib, ic));
            }
          }
        }
      }
    }
  }

  return triangulation;
}

//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::triangulate_graham_scan_3d(const std::vector<Point>& pm)
{
  std::vector<std::vector<Point>> triangulation =
    _triangulate_graham_scan_3d(pm);

#ifdef DOLFIN_ENABLE_GEOMETRY_DEBUGGING

  if (cgal_triangulation_has_degenerate(triangulation))
    throw std::runtime_error("triangulation contains degenerate tetrahedron");

  if (cgal_triangulation_overlap(triangulation))
  {
    throw std::runtime_error("tetrahedrons overlap");
  }

  double volume = .0;
  for (const std::vector<Point>& tet : triangulation)
  {
    assert(tet.size() == 4);
    const double tet_volume = std::abs(orient3d(tet[0], tet[1], tet[2], tet[3]))/6.0;
    volume += tet_volume;
  }

  const double reference_volume = cgal_polyhedron_volume(pm);

  if (std::abs(volume - reference_volume) > 1e-14)
    throw std::runtime_error("computed volume %f, but reference volume is %f (diff %e)");

#endif
  return triangulation;
}
//-----------------------------------------------------------------------------
bool ConvexTriangulation::selfintersects(const std::vector<std::vector<Point>>& p)
{
  for (std::size_t i = 0; i < p.size(); i++)
  {
    for (std::size_t j = i+1; j < p.size(); j++)
    {
      assert(p[i].size() == p[j].size());
      if (p[i].size() == 4)
      {
	if (CollisionPredicates::collides_tetrahedron_tetrahedron_3d(p[i][0],
								     p[i][1],
								     p[i][2],
								     p[i][3],
								     p[j][0],
								     p[j][1],
								     p[j][2],
								     p[j][3]))
	{
	  auto intersection =
	    IntersectionConstruction::intersection_tetrahedron_tetrahedron_3d(p[i][0],
									      p[i][1],
									      p[i][2],
									      p[i][3],
									      p[j][0],
									      p[j][1],
									      p[j][2],
									      p[j][3]);
	  if (intersection.size() > 3)
	  {
	    for (std::size_t k = 3; k < intersection.size(); k++)
	    {
	      if (orient3d(intersection[0],
			   intersection[1],
			   intersection[2],
			   intersection[k]) != 0)
	      {
		return true;
	      }
	    }
	  }
	}
      }
      else if (p[i].size() == 3)
      {
	throw std::runtime_error("Not implemented");
      }
    }
  }

  return false;
}
//-----------------------------------------------------------------------------
