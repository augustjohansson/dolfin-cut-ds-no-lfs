// Copyright (C) 2006-2014 Anders Logg
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
// Modified by Garth N. Wells 2006
// Modified by Andre Massing 2009
//
// First added:  2006-06-12
// Last changed: 2017-09-28

#ifndef __POINT_H
#define __POINT_H

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <string>
#include <ostream>
#include <sstream>

// Optional: aggressive inlining hint for hot builds
#if defined(_MSC_VER)
#define SIMPEX_FORCEINLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define SIMPEX_FORCEINLINE inline __attribute__((always_inline))
#else
#define SIMPEX_FORCEINLINE inline
#endif

namespace dolfin
{

  /// A Point represents a point in R^3 with coordinates x, y, z, or
  /// alternatively, a vector in R^3, supporting standard operations.
  /// This is the optimised implementation from the simpex geometry library.
  class Point
  {
  public:
    // Construction / assignment: keep trivial where possible
    constexpr explicit Point(double x = 0.0, double y = 0.0, double z = 0.0) noexcept
      : _x{{x, y, z}} {}

    // Copy the first dim entries (dim expected <= 3).
    explicit Point(std::size_t dim, const double* x) noexcept : _x{{0.0, 0.0, 0.0}}
    {
      assert(dim <= 3);
      assert(x != nullptr || dim == 0);
      for (std::size_t i = 0; i < dim; ++i) _x[i] = x[i];
    }

    constexpr Point(const Point&) noexcept = default;
    constexpr Point(Point&&) noexcept = default;
    constexpr Point& operator=(const Point&) noexcept = default;
    constexpr Point& operator=(Point&&) noexcept = default;
    ~Point() = default;

    // Element access
    SIMPEX_FORCEINLINE double& operator[](std::size_t i) noexcept
    { assert(i < 3); return _x[i]; }

    SIMPEX_FORCEINLINE const double& operator[](std::size_t i) const noexcept
    { assert(i < 3); return _x[i]; }

    // Named coordinates
    SIMPEX_FORCEINLINE constexpr double x() const noexcept { return _x[0]; }
    SIMPEX_FORCEINLINE constexpr double y() const noexcept { return _x[1]; }
    SIMPEX_FORCEINLINE constexpr double z() const noexcept { return _x[2]; }

    // Raw pointer access
    SIMPEX_FORCEINLINE double* coordinates() noexcept { return _x.data(); }
    SIMPEX_FORCEINLINE const double* coordinates() const noexcept { return _x.data(); }

    SIMPEX_FORCEINLINE constexpr std::array<double, 3> array() const noexcept { return _x; }

    // Arithmetic
    SIMPEX_FORCEINLINE constexpr Point operator+(const Point& p) const noexcept
    { return Point(_x[0] + p._x[0], _x[1] + p._x[1], _x[2] + p._x[2]); }

    SIMPEX_FORCEINLINE constexpr Point operator-(const Point& p) const noexcept
    { return Point(_x[0] - p._x[0], _x[1] - p._x[1], _x[2] - p._x[2]); }

    SIMPEX_FORCEINLINE Point& operator+=(const Point& p) noexcept
    { _x[0] += p._x[0]; _x[1] += p._x[1]; _x[2] += p._x[2]; return *this; }

    SIMPEX_FORCEINLINE Point& operator-=(const Point& p) noexcept
    { _x[0] -= p._x[0]; _x[1] -= p._x[1]; _x[2] -= p._x[2]; return *this; }

    SIMPEX_FORCEINLINE constexpr Point operator-() const noexcept
    { return Point(-_x[0], -_x[1], -_x[2]); }

    SIMPEX_FORCEINLINE constexpr Point operator*(double a) const noexcept
    { return Point(a*_x[0], a*_x[1], a*_x[2]); }

    SIMPEX_FORCEINLINE Point& operator*=(double a) noexcept
    { _x[0] *= a; _x[1] *= a; _x[2] *= a; return *this; }

    SIMPEX_FORCEINLINE Point operator/(double a) const noexcept
    { return Point(_x[0]/a, _x[1]/a, _x[2]/a); }

    SIMPEX_FORCEINLINE Point& operator/=(double a) noexcept
    { _x[0] /= a; _x[1] /= a; _x[2] /= a; return *this; }

    SIMPEX_FORCEINLINE constexpr bool operator==(const Point& p) const noexcept
    { return _x[0] == p._x[0] && _x[1] == p._x[1] && _x[2] == p._x[2]; }

    SIMPEX_FORCEINLINE constexpr bool operator!=(const Point& p) const noexcept
    { return !(*this == p); }

    // Lexicographic less-than for use in ordered containers.
    SIMPEX_FORCEINLINE constexpr bool operator<(const Point& a) const noexcept
    {
      if (_x[0] != a._x[0]) return _x[0] < a._x[0];
      if (_x[1] != a._x[1]) return _x[1] < a._x[1];
      return _x[2] < a._x[2];
    }

    SIMPEX_FORCEINLINE constexpr double squared_distance(const Point& p) const noexcept
    {
      const double dx = _x[0] - p._x[0];
      const double dy = _x[1] - p._x[1];
      const double dz = _x[2] - p._x[2];
      return dx*dx + dy*dy + dz*dz;
    }

    SIMPEX_FORCEINLINE double distance(const Point& p) const noexcept
    { return std::sqrt(squared_distance(p)); }

    SIMPEX_FORCEINLINE constexpr double squared_norm() const noexcept
    { return _x[0]*_x[0] + _x[1]*_x[1] + _x[2]*_x[2]; }

    SIMPEX_FORCEINLINE double norm() const noexcept
    { return std::sqrt(squared_norm()); }

    SIMPEX_FORCEINLINE constexpr Point cross(const Point& p) const noexcept
    {
      return Point(
        _x[1]*p._x[2] - _x[2]*p._x[1],
        _x[2]*p._x[0] - _x[0]*p._x[2],
        _x[0]*p._x[1] - _x[1]*p._x[0]
      );
    }

    SIMPEX_FORCEINLINE constexpr double dot(const Point& p) const noexcept
    { return _x[0]*p._x[0] + _x[1]*p._x[1] + _x[2]*p._x[2]; }

    Point rotate(const Point& a, double theta) const
    {
      assert(std::abs(a.norm() - 1.0) < 3.0e-16);

      const double s = std::sin(theta);
      const double c = std::cos(theta);
      const double one_c = 1.0 - c;

      const double ax = a[0], ay = a[1], az = a[2];
      const double xv = _x[0], yv = _x[1], zv = _x[2];

      const double d = ax*xv + ay*yv + az*zv;

      const double cx = ay*zv - az*yv;
      const double cy = az*xv - ax*zv;
      const double cz = ax*yv - ay*xv;

      return Point(xv*c + cx*s + ax*d*one_c,
                   yv*c + cy*s + ay*d*one_c,
                   zv*c + cz*s + az*d*one_c);
    }

    std::string str() const
    {
      std::stringstream ss;
      ss << x() << ' ' << y() << ' ' << z();
      return ss.str();
    }

    // Backwards-compatibility overload (verbose parameter is ignored).
    std::string str(bool /*verbose*/) const { return str(); }

  private:
    std::array<double, 3> _x;
  };

  /// Multiplication with scalar
  SIMPEX_FORCEINLINE constexpr Point operator*(double a, const Point& p) noexcept
  { return p * a; }

  /// Output of Point to stream
  inline std::ostream& operator<<(std::ostream& stream, const Point& point)
  { stream << point.str(false); return stream; }

} // namespace dolfin

// Provide simpex::Point as an alias for dolfin::Point so that the new
// optimised geometry routines (which use namespace simpex) work without
// any type-conversion overhead.
namespace simpex
{
  using Point = dolfin::Point;
}

#endif
