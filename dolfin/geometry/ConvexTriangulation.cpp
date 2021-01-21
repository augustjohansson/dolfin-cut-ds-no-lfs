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
//
// First added:  2016-06-01
// Last changed: 2018-02-09

#include <algorithm>
#include <numeric>
#include <tuple>
#include <set>
#include "predicates.h"
#include "GeometryPredicates.h"
#include "GeometryTools.h"
#include "CollisionPredicates.h"
#include "IntersectionConstruction.h"
#include "ConvexTriangulation.h"
#include "CGALExactArithmetic.h"

//#ifdef DOLFIN_GEOMETRY_PRINT
#include "GeometryDebugging.h"
//#endif

using namespace dolfin;

//-----------------------------------------------------------------------------
namespace
{
  //-----------------------------------------------------------------------------
  // Create a unique list of points in the sense that |p-q| > tol in each dimension
  std::vector<Point>
  unique_points(const std::vector<Point>& input_points,
		std::size_t gdim,
		double tol)
  {
    std::vector<Point> points;
#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
  GeometryDebugging::print(input_points);
#endif

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
#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
  GeometryDebugging::print(points);
#endif

    return points;
  }
  //------------------------------------------------------------------------------
  // Check if q lies between p0 and p1. p0, p1 and q are assumed to be colinear
  bool is_between(Point p0, Point p1, Point q)
  {
    const double sqnorm = (p1-p0).squared_norm();
    return (p0-q).squared_norm() < sqnorm && (p1-q).squared_norm() < sqnorm;
  }
  //------------------------------------------------------------------------------
  // Return the indices to the points that forms the polygon that is
  // the convex hull of the points. The points are assumed to be coplanar
  std::vector<std::pair<std::size_t, std::size_t>>
  compute_convex_hull_planar(const std::vector<Point>& points)
  {
#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
#endif
    // Try to compute normal
    Point normal;
    bool found = false;
    for (std::size_t i = 0; i < points.size() and !found; ++i)
    {
      for (std::size_t j = i+1; j < points.size() and !found; ++j)
      {
	for (std::size_t k = j+1; k < points.size() and !found; ++k)
	{
	  normal = GeometryTools::cross_product(points[0], points[1], points[2]);
	  const double length = normal.norm();
	  if (length > 0.0)
	  {
	    normal /= length;
	    found = true;
	  }
	}
      }
    }

    dolfin_assert(found);

    std::vector<std::pair<std::size_t, std::size_t>> edges;

    // Filter out points which are in the interior of the
    // convex hull of the planar points.
    for (std::size_t i = 0; i < points.size(); i++)
    {
      for (std::size_t j = i+1; j < points.size(); j++)
      {
        // Form at plane of i, j  and i + the normal of plane
        const Point r = points[i]+normal;

        // search for the first point which is not in the
        // i, j, p plane to determine sign of orientation
        double edge_orientation = 0;
        {
          std::size_t a = 0;
          while (edge_orientation == 0)
          {
            if (a != i && a != j)
            {
              edge_orientation = orient3d(points[i],
                                          points[j],
                                          r,
                                          points[a]);
            }
            a++;
          }
        }

        bool on_convex_hull = true;
	std::vector<std::size_t> colinear;
        for (std::size_t p = 0; p < points.size(); p++)
        {
          if (p != i && p != j)
          {
            const double orientation = orient3d(points[i],
                                                points[j],
                                                r,
                                                points[p]);

	    if (orientation == 0.0)
	    {
	      colinear.push_back(p);
	    }

            // Sign change: triangle is not on convex hull
            if (edge_orientation * orientation < 0)
            {
              on_convex_hull = false;
            }
          }
        }

        if (on_convex_hull)
        {
	  if (!colinear.empty())
	  {
	    // Several points are colinear. Only add if i and j are
	    // the 1d convex hull of the colinear points
	    bool is_linear_convex_hull = true;
	    for (std::size_t q : colinear)
	    {
	      if (!is_between(points[i], points[j], points[q]))
	      {
		is_linear_convex_hull = false;
		break;
	      }
	    }

	    if (is_linear_convex_hull)
	    {
	      edges.push_back(std::make_pair(i,  j));
	    }
	  }
	  else
	  {
	    edges.push_back(std::make_pair(i, j));
	  }
        }
      }
    }

    return edges;
  }

  //-----------------------------------------------------------------------------
  std::vector<std::pair<double, std::size_t>>
  graham_scan_2d(const std::vector<Point>& points)
  {
#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
#endif
    // Sometimes we can get an extra point on an edge: a-----c--b. This
    // point c may cause problems for the graham scan. To avoid this,
    // use an extra center point.  Use this center point and point no 0
    // as reference for the angle calculation
    Point pointscenter = points[0];
    for (std::size_t m = 1; m < points.size(); ++m)
      pointscenter += points[m];
    pointscenter /= points.size();

    // Reference
    const Point ref = points[0] - pointscenter;

    // Calculate and store angles
    std::vector<std::pair<double, std::size_t>> order;
    for (std::size_t m = 1; m < points.size(); ++m)
    {
      const double A = orient2d(pointscenter, points[0], points[m]);
      const Point s = points[m] - pointscenter;
      double alpha = std::atan2(A, s.dot(ref));
      if (alpha < 0)
	alpha += 2.0*DOLFIN_PI;
      order.emplace_back(alpha, m);
    }

    // Sort angles
    std::sort(order.begin(), order.end());
    return order;
  }
} // end empty namespace
//------------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::triangulate(const std::vector<Point>& p,
                                 std::size_t gdim,
                                 std::size_t tdim)
{
#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << ' '<<gdim<<' '<<tdim<<'\n';
#endif
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
  else if (tdim == 2 && gdim == 3)
  {
    return triangulate_graham_scan_2d_3d(p);
  }
  else if (tdim == 3 && gdim == 3)
  {
    return triangulate_graham_scan_3d(p);
  }

  dolfin_error("ConvexTriangulation.cpp",
               "triangulate convex polyhedron",
               "Triangulation of polyhedron of topological dimension %u and geometric dimension %u not implemented", tdim, gdim);

  return std::vector<std::vector<Point>>();
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::triangulate_1d(const std::vector<Point>& p,
				    std::size_t gdim)
{
  return _triangulate_1d(p, gdim);
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_1d(const std::vector<Point>& p,
				     std::size_t gdim)
{
#ifdef DOLFIN_GEOMETRY_PRINT
    std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
#endif
  // A convex polyhedron of topological dimension 1 can not have more
  // than two points. If more, they must be collinear (more or
  // less). This can happen due to tolerances in
  // IntersectionConstruction::intersection_segment_segment_2d.

  if (gdim != 2)
  {
    dolfin_error("ConvexTriangulation.cpp",
		 "triangulate topological 1d",
		 "Function is only implemented for gdim = 2");
  }

  const std::vector<Point> unique_p = unique_points(p, gdim, DOLFIN_EPS);

  if (unique_p.size() == 2)
  {
    // Return the point list. Since it is unique it is also
    // non-degenerate
    std::vector<std::vector<Point>> t { unique_p };
    return t;
  }
  else if (unique_p.size() < 2)
  {
    // Return empty if 0 or 1 point
    return std::vector<std::vector<Point>>();
  }
  else
  {
    // Here unique_p.size() > 2. Make sure the points are
    // approximately collinear
    bool collinear = true;
    for (std::size_t i = 2; i < unique_p.size(); ++i)
    {
      const double o = orient2d(unique_p[0], unique_p[1], unique_p[i]);
      if (std::abs(o) > DOLFIN_EPS)
      {
	collinear = false;
	break;
      }
    }

    if (!collinear)
    {
      dolfin_error("ConvexTriangulation.cpp",
		   "triangulate topological 1d",
		   "A convex polyhedron of topological dimension 1 can not have more than 2 points unless they are collinear.");
    }

    // Return extremal
    const Point v = unique_p[1] - unique_p[0];
    std::vector<std::pair<double, std::size_t>> order;
    order.emplace_back(0.0, 0);
    order.emplace_back(1.0, 1);
    for (std::size_t i = 2; i < unique_p.size(); ++i)
      order.emplace_back(v.dot(unique_p[i]-unique_p[0]), i);

    std::sort(order.begin(), order.end());

    // Return first and last
    return {{ unique_p[order.front().second],
	  unique_p[order.back().second] }};
  }
}
//------------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::triangulate_graham_scan_2d(const std::vector<Point>& input_points)
{
  const std::vector<std::vector<Point>> triangulation
    = _triangulate_graham_scan_2d_3d(input_points);

  return triangulation;
}
//------------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_graham_scan_2d(const std::vector<Point>& input_points)
{
#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
#endif
  dolfin_assert(GeometryPredicates::is_finite(input_points));
  const std::size_t tdim = 2;
  const std::size_t gdim = 2;

  // Make sure the input points are unique. We assume this has
  // negligble effect on the area
  std::vector<Point> points = unique_points(input_points, gdim, DOLFIN_EPS);

  // Too few points
  if (points.size() < 3)
    return std::vector<std::vector<Point>>();

  // Exact one triangle
  if (points.size() == 3)
  {
    std::vector<std::vector<Point>> triangulation;
    if (!GeometryPredicates::is_degenerate(points, tdim, gdim))
      triangulation.push_back(points);
    return triangulation;
  }

  // Graham scan
  const std::vector<std::pair<double, std::size_t>> order
    = graham_scan_2d(points);

  // Tessellate
  std::vector<std::vector<Point>> triangulation;
  for (std::size_t m = 0; m < order.size()-1; ++m)
  {
    const std::vector<Point> tri {{ points[0],
	  points[order[m].second],
	  points[order[m + 1].second] }};
    if (!GeometryPredicates::is_degenerate(tri, tdim, gdim))
      triangulation.push_back(tri);
  }

  return triangulation;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::triangulate_graham_scan_2d_3d(const std::vector<Point>& input_points)
  {
  const std::vector<std::vector<Point>> triangulation
    = _triangulate_graham_scan_2d_3d(input_points);

  return triangulation;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_graham_scan_2d_3d(const std::vector<Point>& input_points)
{
#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
  GeometryDebugging::print(input_points);
#endif
  dolfin_assert(GeometryPredicates::is_finite(input_points));
  const std::size_t tdim = 2;
  const std::size_t gdim = 3;

  // Make sure the input points are unique. We assume this has
  // negligble effect on the area
  std::vector<Point> points = unique_points(input_points, gdim, DOLFIN_EPS);

  // Too few points
  if (points.size() < 3)
    return std::vector<std::vector<Point>>();

  // Exactly one triangle
  if (points.size() == 3)
  {
    std::vector<std::vector<Point>> triangulation;
    if (!GeometryPredicates::is_degenerate(points, tdim, gdim))
      triangulation.push_back(points);
    return triangulation;
  }

  // Check that all points are (almost) coplanar
  for (std::size_t i = 3; i < points.size(); ++i)
  {
    const double o = orient3d(points[0], points[1], points[2],
			      points[i]);
    if (std::abs(o) > DOLFIN_EPS_LARGE)
    {
      // dolfin_error("ConvexTriangulation.cpp",
      //   	   "triangulate tdim 2, gdim 3",
      //   	   "points are not coplanar since o=%e", o);
      std::cout << __FUNCTION__<<' '<<__LINE__<<" warning, perhaps not coplanar: " << std::abs(o) << '\n';
    }
  }

  // Project to major axis and perform Graham scan.
  const Point n = GeometryTools::cross_product(points[0], points[1], points[2]);
  const std::size_t major_axis = GeometryTools::major_axis_3d(n);

  // Project
  std::vector<Point> proj(points.size());
  for (std::size_t i = 0; i < points.size(); ++i)
    proj[i] = GeometryTools::project_to_plane_3d(points[i], major_axis);

  // Graham scan in projected plane
  const std::vector<std::pair<double, std::size_t>> order
    = graham_scan_2d(proj);

  // Tessellate
  std::vector<std::vector<Point>> triangulation;
  for (std::size_t m = 0; m < order.size()-1; ++m)
  {
    // FIXME: We could consider only triangles with area > tolerance here.
    const std::vector<Point> tri {{ points[0],
	  points[order[m].second],
	  points[order[m + 1].second] }};
    if (!GeometryPredicates::is_degenerate(tri, tdim, gdim))
      triangulation.push_back(tri);
  }

#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout<<"   " << __FUNCTION__ << ' ' << __LINE__ << " done\n";
  GeometryDebugging::print(triangulation);
#endif
  
  return triangulation;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_graham_scan_3d(const std::vector<Point>& input_points)
{
 #ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
#endif
 dolfin_assert(GeometryPredicates::is_finite(input_points));
  const std::size_t tdim = 3;
  const std::size_t gdim = 3;

  // Make sure the input points are unique. We assume this has
  // negligble effect on volume
  std::vector<Point> points = unique_points(input_points, gdim, DOLFIN_EPS);

  std::vector<std::vector<Point>> triangulation;

  if (points.size() < 4)
  {
    // Empty
    return triangulation;
  }
  else if (points.size() == 4)
  {
    // Single tetrahedron
    if (!GeometryPredicates::is_degenerate(points, tdim, gdim))
      triangulation.push_back(points);
    return triangulation;
  }
  else
  {
    const Point polyhedroncenter = std::accumulate(points.begin()+1, points.end(),
						   points.front()) / points.size();

    // FIXME: Better data structure than set?
    std::set<std::tuple<std::size_t, std::size_t, std::size_t> > checked;

    // Loop over all triplets
    for (std::size_t i = 0; i < points.size()-2; ++i)
    {
      for (std::size_t j = i+1; j < points.size()-1; ++j)
      {
        for (std::size_t k = j+1; k < points.size(); ++k)
        {
	  if (checked.emplace(std::make_tuple(i, j, k)).second)
	  {
            // Test for the special case where i, j, k are collinear
            {
              const Point ij = points[j] - points[i];
              const Point ik = points[k] - points[i];
	      if ( -(std::abs((ij).dot(ik)) - ij.norm()*ik.norm()) < DOLFIN_EPS)
                continue;
            }

            // Check whether all other points are on one side of this
            // (i,j,k) facet, i.e. we're on the convex
            // hull. Initialize as true for the case of only three
            // coplanar points.
	    bool on_convex_hull = true;

	    // Use orient3d to determine if the plane (i,j,k) is on the
	    // convex hull.
	    std::vector<std::size_t> coplanar { i, j, k };
	    double previous_orientation;
	    bool first = true;

	    for (std::size_t m = 0; m < points.size(); ++m)
	    {
	      if (m != i and m != j and m != k)
	      {
		const double orientation = orient3d(points[i],
						    points[j],
						    points[k],
						    points[m]);
		// Save point index if we find coplanar points
		if (orientation == 0.0)
		  coplanar.push_back(m);
                else
                {
                  if (first)
                  {
                    previous_orientation = orientation;
                    first = false;
                  }
                  else
                  {
                    // Sign change: triangle is not on convex hull
		    if ((previous_orientation < 0.0 and orientation > 0.0) or
			(previous_orientation > 0.0 and orientation < 0.0))
                    {
                      on_convex_hull = false;
                    }
                  }
		}
	      }
	    }

	    if (on_convex_hull)
	    {
	      if (coplanar.size() == 3)
	      {
		// Form one tetrahedron (can be degenerate due to
		// calculation of polyhedroncenter)
		const std::vector<Point> cand { points[i], points[j], points[k],
					    polyhedroncenter };

		// FIXME: Here we could include if determinant is
		// sufficiently large
		if (!GeometryPredicates::is_degenerate(cand, tdim, gdim))
		{
		  triangulation.push_back(cand);
	        }
	      }
	      else // At least four coplanar points
	      {
		// Tessellate as in the triangle-triangle intersection
		// case: First sort points using a Graham scan, then
		// connect to form triangles. Finally form tetrahedra
		// using the center of the polyhedron.

		// Use the center of the coplanar points and point no 0
		// as reference for the angle calculation

		std::vector<Point> coplanar_points;
		for (std::size_t i : coplanar)
		  coplanar_points.push_back(points[i]);

		std::vector<std::pair<std::size_t, std::size_t>> coplanar_convex_hull =
		  compute_convex_hull_planar(coplanar_points);

		const Point coplanar_center
		  = std::accumulate(coplanar_points.begin()+1, coplanar_points.end(),
				    coplanar_points.front()) / coplanar_points.size();

		// Tessellate
		for (const std::pair<std::size_t, std::size_t>& edge : coplanar_convex_hull)
		{
		  const std::vector<Point> cand {{ polyhedroncenter,
			coplanar_center,
			coplanar_points[edge.first],
			coplanar_points[edge.second] }};
		  if (!GeometryPredicates::is_degenerate(cand, tdim, gdim))
		  {
		    triangulation.push_back(cand);
		  }

		  // Mark all combinations of the coplanar vertices as
		  // checked to avoid duplicating triangles
		  std::sort(coplanar.begin(), coplanar.end());

		  for (std::size_t i = 0; i < coplanar.size()-2; ++i)
		    for (std::size_t j = i+1; j < coplanar.size()-1; ++j)
		      for (std::size_t k = j+1; k < coplanar.size(); k++)
			checked.emplace( std::make_tuple(coplanar[i], coplanar[j], coplanar[k]) );
		}
	      } // end coplanar.size() > 3
	    } // end on_convexhull
	  }
	}
      }
    }

    return triangulation;
  }
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::triangulate_graham_scan_3d(const std::vector<Point>& pm)
{
#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
#endif
  std::vector<std::vector<Point>> triangulation =
    _triangulate_graham_scan_3d(pm);

#ifdef DOLFIN_ENABLE_GEOMETRY_DEBUGGING

  if (cgal_triangulation_has_degenerate(triangulation))
    dolfin_error("ConvexTriangulation.cpp",
		 "verify convex triangulation",
		 "triangulation contains degenerate tetrahedron");

  const std::size_t tdim = 3;
  const double reference_volume = cgal_polyhedron_volume(pm, tdim);
  // If volume is really small it's no need to check overlap
  if (std::abs(reference_volume) > DOLFIN_EPS and
      cgal_triangulation_overlap(triangulation))
  {
    std::cout << "\n\n cgal_triangulation_overlap error line " << __LINE__<<"\n";

    for (const Point& p : pm)
      std::cout << tools::plot(p);
    std::cout << '\n';
    double v=0;
    for (const std::vector<Point>& tet : triangulation)
    {
      v += std::abs(orient3d(tet[0], tet[1], tet[2], tet[3]))/6.0;
      std::cout << tools::drawtriangle(tet)<<"% "<<orient3d(tet[0], tet[1], tet[2], tet[3])/6.0<<'\n';
    }
    std::cout << '\n';
    std::cout << "approx volume " << v << '\n';

    std::cout << "reference volume " << reference_volume << '\n';

    dolfin_error("ConvexTriangulation.cpp",
		 "verify convex triangulation",
		 "tetrahedrons overlap");
  }


  double volume = 0.0;
  for (const std::vector<Point>& tet : triangulation)
  {
    dolfin_assert(tet.size() == 4);
    const double tet_volume = std::abs(orient3d(tet[0], tet[1], tet[2], tet[3]))/6.0;
    volume += tet_volume;
  }

  if (std::abs(volume - reference_volume) > DOLFIN_EPS_LARGE)
  {
    for (const Point& p : pm)
      std::cout << tools::plot(p);
    std::cout << '\n';

    const double error = std::abs(reference_volume-volume);
    std::cout << "error vol " << volume<<' '<<reference_volume<<' '<<error<<'\n';
    // dolfin_error("ConvexTriangulation.cpp",
    // 		 "verifying convex triangulation",
    // 		 "computed volume %e, but reference volume is %e (diff %e)",
    // 		 volume, reference_volume, error);
  }
#endif

  return triangulation;
}
//-----------------------------------------------------------------------------
bool ConvexTriangulation::selfintersects(const std::vector<std::vector<Point>>& p)
{
#ifdef DOLFIN_GEOMETRY_PRINT
  std::cout << __FUNCTION__ << ' ' << __LINE__ << '\n';
#endif
  for (std::size_t i = 0; i < p.size(); i++)
  {
    for (std::size_t j = i+1; j < p.size(); j++)
    {
      dolfin_assert(p[i].size() == p[j].size());
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
	      // FIXME: Note that this fails if the first three points are colinear!
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
	dolfin_not_implemented();
      }
    }
  }

  return false;
}
//-----------------------------------------------------------------------------
