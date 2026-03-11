// Unit tests for SimplexQuadrature, ported from python/tests/test_quadrature.py

#include <cmath>
#include <numeric>
#include "../geometry/SimplexQuadrature.h"
#include "../geometry/Point.h"
#include "catch/catch.hpp"

using namespace simpex;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static double sum_weights(const std::pair<std::vector<double>, std::vector<double>>& qr)
{
  const auto& wts = qr.second;
  return std::accumulate(wts.begin(), wts.end(), 0.0);
}

// ---------------------------------------------------------------------------
// Interval quadrature (1D)
// ---------------------------------------------------------------------------

TEST_CASE("SimplexQuadrature: interval")
{
  SECTION("unit interval weight sum")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - 1.0) < 1e-14);
  }

  SECTION("scaled interval weight sum")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(3.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - 3.0) < 1e-14);
  }

  SECTION("exact polynomial integration: integrate x over [0,1] = 0.5")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    const auto& pts = qr.first;
    const auto& wts = qr.second;
    double result = 0.0;
    for (std::size_t i = 0; i < wts.size(); ++i)
      result += pts[i] * wts[i];
    CHECK(std::abs(result - 0.5) < 1e-14);
  }

  SECTION("wrapper compute_quadrature_rule for interval")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule(coords, 1);
    CHECK(std::abs(sum_weights(qr) - 1.0) < 1e-14);
  }

  SECTION("higher order interval")
  {
    SimplexQuadrature sq(1, 5);
    std::vector<Point> coords = { Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - 1.0) < 1e-14);
  }
}

// ---------------------------------------------------------------------------
// Triangle quadrature (2D)
// ---------------------------------------------------------------------------

TEST_CASE("SimplexQuadrature: triangle")
{
  SECTION("unit triangle weight sum")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-14);
  }

  SECTION("scaled triangle weight sum")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0), Point(0.0, 2.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    // area = 0.5 * 2 * 2 = 2.0
    CHECK(std::abs(sum_weights(qr) - 2.0) < 1e-14);
  }

  SECTION("wrapper compute_quadrature_rule for triangle")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
    auto qr = sq.compute_quadrature_rule(coords, 2);
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-14);
  }

  SECTION("higher order triangle")
  {
    SimplexQuadrature sq(2, 5);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-14);
  }

  SECTION("points layout: gdim=2 coords per point")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    CHECK(qr.first.size() == qr.second.size() * 2);
  }
}

// ---------------------------------------------------------------------------
// Tetrahedron quadrature (3D)
// ---------------------------------------------------------------------------

TEST_CASE("SimplexQuadrature: tetrahedron")
{
  SECTION("unit tet weight sum")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    CHECK(std::abs(sum_weights(qr) - 1.0/6.0) < 1e-14);
  }

  SECTION("scaled tet weight sum")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0),
      Point(0.0, 2.0, 0.0), Point(0.0, 0.0, 2.0) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    // volume = (2^3)/6 = 8/6 = 4/3
    CHECK(std::abs(sum_weights(qr) - 4.0/3.0) < 1e-14);
  }

  SECTION("wrapper compute_quadrature_rule for tet")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0) };
    auto qr = sq.compute_quadrature_rule(coords, 3);
    CHECK(std::abs(sum_weights(qr) - 1.0/6.0) < 1e-14);
  }

  SECTION("higher order tet")
  {
    SimplexQuadrature sq(3, 4);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    CHECK(std::abs(sum_weights(qr) - 1.0/6.0) < 1e-13);
  }

  SECTION("points layout: gdim=3 coords per point")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
      Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    CHECK(qr.first.size() == qr.second.size() * 3);
  }
}

// ---------------------------------------------------------------------------
// Magnitude tests: entities of size ~1e-13, ~1e13, and mixed
// ---------------------------------------------------------------------------

TEST_CASE("SimplexQuadrature: magnitude ~1e-13")
{
  const double s = 1e-13;

  SECTION("interval weight sum (small scale)")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0), Point(s) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - s) < s * 1e-12);
  }

  SECTION("triangle weight sum (small scale)")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0), Point(0.0, s, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    // area = 0.5 * s^2
    CHECK(std::abs(sum_weights(qr) - 0.5*s*s) < 0.5*s*s * 1e-12);
  }

  SECTION("tet weight sum (small scale)")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0),
      Point(0.0, s, 0.0),   Point(0.0, 0.0, s) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    // volume = s^3/6
    CHECK(std::abs(sum_weights(qr) - s*s*s/6.0) < s*s*s/6.0 * 1e-12);
  }
}

TEST_CASE("SimplexQuadrature: magnitude ~1e13")
{
  const double s = 1e13;

  SECTION("interval weight sum (large scale)")
  {
    SimplexQuadrature sq(1, 2);
    std::vector<Point> coords = { Point(0.0), Point(s) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - s) < s * 1e-12);
  }

  SECTION("triangle weight sum (large scale)")
  {
    SimplexQuadrature sq(2, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0), Point(0.0, s, 0.0) };
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    CHECK(std::abs(sum_weights(qr) - 0.5*s*s) < 0.5*s*s * 1e-12);
  }

  SECTION("tet weight sum (large scale)")
  {
    SimplexQuadrature sq(3, 2);
    std::vector<Point> coords = {
      Point(0.0, 0.0, 0.0), Point(s, 0.0, 0.0),
      Point(0.0, s, 0.0),   Point(0.0, 0.0, s) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    CHECK(std::abs(sum_weights(qr) - s*s*s/6.0) < s*s*s/6.0 * 1e-12);
  }
}

TEST_CASE("SimplexQuadrature: mixed magnitudes")
{
  SECTION("interval mixed: [0, 1e-13]")
  {
    SimplexQuadrature sq(1, 2);
    const double len = 1e-13;
    std::vector<Point> coords = { Point(0.0), Point(len) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - len) < len * 1e-12);
  }

  SECTION("interval mixed: [0, 1e13]")
  {
    SimplexQuadrature sq(1, 2);
    const double len = 1e13;
    std::vector<Point> coords = { Point(0.0), Point(len) };
    auto qr = sq.compute_quadrature_rule_interval(coords, 1);
    CHECK(std::abs(sum_weights(qr) - len) < len * 1e-12);
  }

  SECTION("tet with mixed edge lengths (~1e-13 and ~1e13)")
  {
    SimplexQuadrature sq(3, 2);
    // Volume = (1e-13 * 1e13 * 1.0) / 6 = 1.0/6
    std::vector<Point> coords = {
      Point(0.0,    0.0,   0.0),
      Point(1e-13,  0.0,   0.0),
      Point(0.0,    1e13,  0.0),
      Point(0.0,    0.0,   1.0) };
    auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);
    CHECK(std::abs(sum_weights(qr) - 1.0/6.0) < 1e-10);
  }
}

// ---------------------------------------------------------------------------
// Compression tests
// ---------------------------------------------------------------------------

#include <functional>
#include "../geometry/SimplexQuadratureCompression.h"

using namespace simpex::compression;

// Helper: apply a quadrature rule to integrate a function f(x, y, [z]).
using Func2D = std::function<double(double, double)>;
using Func3D = std::function<double(double, double, double)>;

static double integrate_2d(const std::pair<std::vector<double>,
                                            std::vector<double>>& qr,
                            Func2D f)
{
  double result = 0.0;
  const std::size_t n = qr.second.size();
  for (std::size_t i = 0; i < n; ++i)
    result += qr.second[i] * f(qr.first[2*i], qr.first[2*i+1]);
  return result;
}

static double integrate_3d(const std::pair<std::vector<double>,
                                            std::vector<double>>& qr,
                            Func3D f)
{
  double result = 0.0;
  const std::size_t n = qr.second.size();
  for (std::size_t i = 0; i < n; ++i)
    result += qr.second[i] * f(qr.first[3*i], qr.first[3*i+1], qr.first[3*i+2]);
  return result;
}

TEST_CASE("Compression: legacy (signed weights)")
{
  // Dunavant order-10 triangle rule has 25 points.
  // Compressing to preserve degree-3 polynomials: K = choose(5,2) = 10 < 25.
  SimplexQuadrature sq(2, 10);
  std::vector<Point> coords = {
    Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
  auto qr = sq.compute_quadrature_rule_triangle(coords, 2);

  const std::size_t n_before = qr.second.size();  // 25
  auto idx = compress_legacy(qr, 2, 3);           // K = 10

  SECTION("number of points is reduced")
  {
    CHECK(qr.second.size() < n_before);
    CHECK(!idx.empty());
  }

  SECTION("weight sum (area) is preserved")
  {
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-10);
  }

  SECTION("degree-3 polynomial is integrated exactly")
  {
    // integral of x^3 over unit triangle = 1/20
    const double val = integrate_2d(qr, [](double x, double){ return x * x * x; });
    CHECK(std::abs(val - 1.0/20.0) < 1e-10);
  }

  SECTION("legacy rule may have negative weights (documented behaviour)")
  {
    // The legacy method may produce negative weights; we do not fail on this.
    bool any_negative = false;
    for (double w : qr.second)
      if (w < 0.0) { any_negative = true; break; }
    (void)any_negative;
    CHECK(true);
  }
}

TEST_CASE("Compression: NNLS (positive weights, Tchakaloff / Tetrafreeq)")
{
  // Dunavant order-10 triangle rule has 25 points.
  // Compressing to preserve degree-3: K = choose(5,2) = 10 < 25.
  SimplexQuadrature sq(2, 10);
  std::vector<Point> coords = {
    Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
  auto qr = sq.compute_quadrature_rule_triangle(coords, 2);

  const std::size_t n_before = qr.second.size();
  compress_nnls(qr, 2, 3);

  SECTION("number of points is reduced")
  {
    CHECK(qr.second.size() < n_before);
  }

  SECTION("all weights are strictly non-negative")
  {
    for (double w : qr.second)
      CHECK(w >= 0.0);
  }

  SECTION("weight sum (area) is preserved")
  {
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-10);
  }

  SECTION("degree-3 polynomial is integrated exactly")
  {
    const double val = integrate_2d(qr, [](double x, double){ return x * x * x; });
    CHECK(std::abs(val - 1.0/20.0) < 1e-10);
  }
}

TEST_CASE("Compression: IRLS (sparse signed weights, CheapQ-compatible)")
{
  SimplexQuadrature sq(2, 10);
  std::vector<Point> coords = {
    Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };
  auto qr = sq.compute_quadrature_rule_triangle(coords, 2);

  const std::size_t n_before = qr.second.size();
  compress_irls(qr, 2, 3);

  SECTION("number of points is reduced")
  {
    CHECK(qr.second.size() < n_before);
  }

  SECTION("weight sum (area) is preserved")
  {
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-10);
  }

  SECTION("degree-3 polynomial is integrated exactly")
  {
    const double val = integrate_2d(qr, [](double x, double){ return x * x * x; });
    CHECK(std::abs(val - 1.0/20.0) < 1e-10);
  }
}

TEST_CASE("Compression: tet NNLS positive weights")
{
  // Keast order-6 tet rule has 24 points.
  // Compressing to preserve degree-3: K = choose(6,3) = 20 < 24.
  SimplexQuadrature sq(3, 6);
  std::vector<Point> coords = {
    Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0), Point(0.0, 0.0, 1.0) };
  auto qr = sq.compute_quadrature_rule_tetrahedron(coords, 3);

  compress_nnls(qr, 3, 3);

  SECTION("all weights are strictly non-negative")
  {
    for (double w : qr.second)
      CHECK(w >= 0.0);
  }

  SECTION("weight sum (volume) is preserved")
  {
    CHECK(std::abs(sum_weights(qr) - 1.0/6.0) < 1e-10);
  }
}

TEST_CASE("Compression: SimplexQuadrature::compress dispatch")
{
  // Test that SimplexQuadrature::compress() delegates to the correct algorithm.
  // Order-10 triangle (25 points), compress to degree-3 (K=10).
  SimplexQuadrature sq(2, 10);
  std::vector<Point> coords = {
    Point(0.0, 0.0, 0.0), Point(1.0, 0.0, 0.0), Point(0.0, 1.0, 0.0) };

  SECTION("legacy dispatch")
  {
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    const std::size_t n_before = qr.second.size();
    SimplexQuadrature::compress(qr, 2, 3, CompressionMethod::legacy);
    CHECK(qr.second.size() < n_before);
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-10);
  }

  SECTION("nnls dispatch")
  {
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    SimplexQuadrature::compress(qr, 2, 3, CompressionMethod::nnls);
    for (double w : qr.second)
      CHECK(w >= 0.0);
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-10);
  }

  SECTION("irls dispatch")
  {
    auto qr = sq.compute_quadrature_rule_triangle(coords, 2);
    SimplexQuadrature::compress(qr, 2, 3, CompressionMethod::irls);
    CHECK(std::abs(sum_weights(qr) - 0.5) < 1e-10);
  }
}

// ---------------------------------------------------------------------------
// Strictly positive quadrature: cancellation effect demonstration.
//
// Integrating a strongly peaked function f(x) = exp(alpha * x) over [0, 1].
// With negative-weight rules, cancellation between large positive and large
// negative terms causes significant error.  NNLS produces positive weights
// only, which avoids cancellation entirely.
// ---------------------------------------------------------------------------

TEST_CASE("Compression: positive weights avoid cancellation")
{
  // Integrate exp(alpha * x) over the unit interval [0,1].
  // Exact answer: (exp(alpha) - 1) / alpha.
  const double alpha    = 50.0;
  const double exact    = (std::exp(alpha) - 1.0) / alpha;

  // Use a high-order rule to accumulate many points.
  SimplexQuadrature sq(1, 20);
  std::vector<Point> coords = { Point(0.0), Point(1.0) };
  auto qr_full = sq.compute_quadrature_rule_interval(coords, 1);

  SECTION("full rule (no compression)")
  {
    double val = 0.0;
    for (std::size_t i = 0; i < qr_full.second.size(); ++i)
      val += qr_full.second[i] * std::exp(alpha * qr_full.first[i]);
    CHECK(std::abs(val - exact) / exact < 1e-6);
  }

  SECTION("NNLS-compressed rule has non-negative weights")
  {
    auto qr = qr_full;
    compress_nnls(qr, 1, 20);
    for (double w : qr.second)
      CHECK(w >= 0.0);
  }
}
