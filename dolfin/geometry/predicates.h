// This is a DOLFIN header file for predicates.cpp which provides
//
//   Routines for Arbitrary Precision Floating-point Arithmetic
//   and Fast Robust Geometric Predicates
//
// by
//
//   Jonathan Richard Shewchuk
//
// Code is placed in the public domain.

#ifndef __PREDICATES_H
#define __PREDICATES_H

#include <map>
#include <tuple>
#include "Point.h"

namespace dolfin
{

  class Point;

  /// Initialize tolerances for exact arithmetic
  void exactinit();

  /// Compute relative orientation of point x wrt segment [a, b]
  double orient1d(double a, double b, double x);

  /// Compute relative orientation of points a, b, c. The orientation
  /// is such that orient2d(a, b, c) > 0 if a, b, c are ordered
  /// counter-clockwise.
  double _orient2d(const double* a, const double* b, const double* c);

  /// Convenience function using dolfin::Point
  double orient2d(const Point& a, const Point& b, const Point& c);

  /// Compute relative orientation of points a, b, c, d. The
  /// orientation is such that orient3d(a, b, c, d) > 0 if a, b, c, d
  /// are oriented according to the left hand rule.
  double _orient3d(const double* a, const double* b, const double* c, const double* d);

  /// Convenience function using dolfin::Point
  double orient3d(const Point& a, const Point& b, const Point& c, const Point& d);

  // // Memoized version
  // struct point_strictly_less
  // {
  //   bool operator()(const dolfin::Point& p0,
  // 		    const dolfin::Point& p1) const
  //   {
  //     return p0[0] < p1[0] and p0[1] < p1[1] and p0[2] < p1[2];
  //   }
  // };

  // typedef std::tuple<Point, Point, Point, Point> Key;

  // struct point_strictly_less
  // {
  //   bool operator()(const Point& _x,
  // 		    const Point& a) const
  //   {
  //     if (_x[0] != a[0]) {
  // 	return _x[0] < a[0];
  //     }

  //     if (_x[1] != a[1]) {
  // 	return _x[1] < a[1];
  //     }

  //     return _x[2] < a[2];
  //   }
  // };

  // struct tuple_strictly_less
  // {
  //   bool operator(const Key& k1, const Key& k2)
  //   {
  //     const auto p0 = std::get<0>(k1);
  //     const auto p1 = std::get<1>(k1);
  //     const auto p2 = std::get<2>(k1);
  //     const auto p3 = std::get<3>(k1);
  //     const auto q0 = std::get<0>(k2);
  //     const auto q1 = std::get<1>(k2);
  //     const auto q2 = std::get<2>(k2);
  //     const auto q3 = std::get<3>(k2);

  //   }
  // };

  //static std::map<Key, double> hash_o3d;
  static std::map<std::tuple<Point, Point, Point, Point>, double> hash_o3d;
  double memoized_orient3d(const Point& a, const Point& b, const Point& c, const Point& d);

  /// Class used for automatic initialization of tolerances at startup.
  /// A global instance is defined inside predicates.cpp to ensure that
  /// the constructor and thus exactinit() is called.

  class PredicateInitialization
  {
  public:

    PredicateInitialization()
    {
      exactinit();
    }

  };

}

#endif
