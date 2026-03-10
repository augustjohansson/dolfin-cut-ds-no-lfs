// Quadrature rule compression for simplices.
//
// Implements three algorithms:
//   1. Legacy (QR-based, Sommariva-Vianello 2015): signed weights.
//   2. NNLS (positive weights, Tchakaloff / Tetrafreeq): all weights >= 0.
//   3. IRLS (sparse, CheapQ-compatible): signed sparse weights.
//
// References:
//   A. Sommariva, M. Vianello,
//   "Compression of multivariate discrete measures and applications",
//   Numerical Functional Analysis and Optimization, 36(9), 2015.
//
//   A. Sommariva, M. Vianello, CheapQ (positive quadrature via compressed
//   Chebyshev Vandermonde systems), https://github.com/alvisesommariva/CheapQ.
//
//   G. Migliorati, T. Pruneri, F. Nobile, "Tetrahedral-free quadrature on
//   simplices via moment fitting (Tetrafreeq)",
//   Journal of Computational Physics 493 (2023).
//   https://www.sciencedirect.com/science/article/abs/pii/S016892742300185X

#include "SimplexQuadratureCompression.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

namespace simpex::compression
{

namespace
{

// ---------------------------------------------------------------------------
// Binomial coefficient n choose k.
// ---------------------------------------------------------------------------
static inline std::size_t choose(std::size_t n, std::size_t k)
{
  if (k > n) return 0;
  if (k == 0 || k == n) return 1;
  if (k > n - k) k = n - k;
  long double res = 1.0L;
  for (std::size_t i = 1; i <= k; ++i)
    res = res * static_cast<long double>(n - k + i)
              / static_cast<long double>(i);
  return static_cast<std::size_t>(std::llround(res));
}

// ---------------------------------------------------------------------------
// Chebyshev polynomials of the first kind T_0 ... T_N evaluated at x.
//
//   T_0(x) = 1
//   T_1(x) = x
//   T_k(x) = 2 x T_{k-1}(x) - T_{k-2}(x), k >= 2
// ---------------------------------------------------------------------------
static std::vector<Eigen::VectorXd>
chebyshev_polynomial(const Eigen::VectorXd& x, std::size_t N)
{
  std::vector<Eigen::VectorXd> T(N + 1, Eigen::VectorXd(x.size()));
  T[0] = Eigen::VectorXd::Ones(x.size());
  if (N >= 1)
  {
    T[1] = x;
    for (std::size_t k = 2; k <= N; ++k)
      T[k] = 2.0 * x.cwiseProduct(T[k - 1]) - T[k - 2];
  }
  return T;
}

// ---------------------------------------------------------------------------
// Graded lexicographic ordering of multi-indices up to total degree N.
// ---------------------------------------------------------------------------
static std::vector<std::vector<std::size_t>>
grlex(std::size_t gdim, std::size_t N)
{
  const std::size_t n_rows = choose(N + gdim, gdim);
  std::vector<std::vector<std::size_t>> P(n_rows,
                                          std::vector<std::size_t>(gdim));
  switch (gdim)
  {
  case 2:
    for (std::size_t s = 0, row = 0; s <= N; ++s)
      for (std::size_t xi = 0; xi <= N; ++xi)
        for (std::size_t yi = 0; yi <= N; ++yi)
          if (xi + yi == s)
          {
            P[row][0] = xi;
            P[row][1] = yi;
            ++row;
          }
    break;
  case 3:
    for (std::size_t s = 0, row = 0; s <= N; ++s)
      for (std::size_t xi = 0; xi <= N; ++xi)
        for (std::size_t yi = 0; yi <= N; ++yi)
          for (std::size_t zi = 0; zi <= N; ++zi)
            if (xi + yi + zi == s)
            {
              P[row][0] = xi;
              P[row][1] = yi;
              P[row][2] = zi;
              ++row;
            }
    break;
  default:
    throw std::runtime_error(
        "grlex ordering only implemented for gdim = 2 or gdim = 3.");
  }
  return P;
}

// ---------------------------------------------------------------------------
// Chebyshev-Vandermonde matrix V of size M x K.
//   V[m, k] = T_{alpha_k}(x_m)   (product over dimensions)
// where M = number of quadrature points, K = choose(N + gdim, gdim).
// ---------------------------------------------------------------------------
static Eigen::MatrixXd
chebyshev_vandermonde(
    const std::pair<std::vector<double>, std::vector<double>>& qr,
    std::size_t gdim,
    std::size_t N)
{
  const auto M = static_cast<Eigen::Index>(qr.second.size());
  std::vector<std::vector<Eigen::VectorXd>> T(gdim);
  Eigen::VectorXd x(M);

  for (std::size_t d = 0; d < gdim; ++d)
  {
    for (Eigen::Index i = 0; i < M; ++i)
      x(i) = qr.first[static_cast<std::size_t>(i) * gdim + d];

    const double xmin = x.minCoeff();
    const double xmax = x.maxCoeff();
    const double hx   = xmax - xmin;
    // Map points to [-1, 1].
    const Eigen::VectorXd xmap =
        (2.0 / hx) * x
        - ((xmin + xmax) / hx) * Eigen::VectorXd::Ones(M);

    T[d] = chebyshev_polynomial(xmap, N);
  }

  const auto P      = grlex(gdim, N);
  const auto K      = static_cast<Eigen::Index>(P.size());
  Eigen::MatrixXd V = Eigen::MatrixXd::Ones(M, K);

  for (Eigen::Index i = 1; i < K; ++i)
  {
    Eigen::VectorXd col = T[0][P[static_cast<std::size_t>(i)][0]];
    for (std::size_t d = 1; d < gdim; ++d)
      col = col.cwiseProduct(T[d][P[static_cast<std::size_t>(i)][d]]);
    V.col(i) = col;
  }

  return V; // M x K
}

// ---------------------------------------------------------------------------
// Write back the selected indices from w_new into qr.
// Returns the vector of retained indices.
// ---------------------------------------------------------------------------
static std::vector<std::size_t>
apply_compressed_weights(
    std::pair<std::vector<double>, std::vector<double>>&       qr,
    const std::pair<std::vector<double>, std::vector<double>>& qr_input,
    const Eigen::VectorXd&                                     w_new,
    std::size_t                                                gdim,
    std::size_t                                                max_pts,
    double                                                     tol)
{
  std::vector<std::size_t> idx;
  idx.reserve(w_new.size());
  for (Eigen::Index i = 0; i < w_new.size(); ++i)
    if (std::abs(w_new[i]) > tol)
      idx.push_back(static_cast<std::size_t>(i));

  // If more than max_pts nonzero entries remain, keep those with the largest
  // absolute weight.
  if (idx.size() > max_pts)
  {
    std::vector<std::pair<double, std::size_t>> mag;
    mag.reserve(idx.size());
    for (auto j : idx)
      mag.emplace_back(std::abs(w_new[static_cast<Eigen::Index>(j)]), j);
    std::nth_element(mag.begin(),
                     mag.begin() + static_cast<Eigen::Index>(max_pts),
                     mag.end(),
                     [](const auto& a, const auto& b) {
                       return a.first > b.first;
                     });
    mag.resize(max_pts);
    idx.clear();
    idx.reserve(max_pts);
    for (const auto& p : mag) idx.push_back(p.second);
    std::sort(idx.begin(), idx.end());
  }

  qr.first.resize(gdim * idx.size());
  qr.second.resize(idx.size());
  for (std::size_t i = 0; i < idx.size(); ++i)
  {
    for (std::size_t d = 0; d < gdim; ++d)
      qr.first[gdim * i + d] = qr_input.first[gdim * idx[i] + d];
    qr.second[i] = w_new[static_cast<Eigen::Index>(idx[i])];
  }
  return idx;
}

// ---------------------------------------------------------------------------
// Lawson-Hanson non-negative least-squares solver.
// Solves  min ||A x - b||_2  subject to  x >= 0.
// ---------------------------------------------------------------------------
static Eigen::VectorXd
nnls_lawson_hanson(const Eigen::MatrixXd& A,
                   const Eigen::VectorXd& b,
                   double                 tol       = 1e-12,
                   int                    max_outer = 5000,
                   int                    max_inner = 5000)
{
  const int m = static_cast<int>(A.rows());
  const int n = static_cast<int>(A.cols());

  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd w = A.transpose() * (b - A * x);

  std::vector<bool> P(n, false);
  std::vector<bool> Z(n, true);

  for (int outer = 0; outer < max_outer; ++outer)
  {
    // Find the free variable with the largest gradient.
    int    t    = -1;
    double maxw = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i)
      if (Z[i] && w[i] > maxw) { maxw = w[i]; t = i; }

    if (t < 0 || maxw <= tol) break;

    P[t] = true;
    Z[t] = false;

    for (int inner = 0; inner < max_inner; ++inner)
    {
      // Collect the active-set column indices.
      std::vector<int> pidx;
      pidx.reserve(n);
      for (int i = 0; i < n; ++i) if (P[i]) pidx.push_back(i);

      const int        p   = static_cast<int>(pidx.size());
      Eigen::MatrixXd  A_P(m, p);
      for (int j = 0; j < p; ++j) A_P.col(j) = A.col(pidx[j]);

      const Eigen::VectorXd z_P = A_P.colPivHouseholderQr().solve(b);

      Eigen::VectorXd s = Eigen::VectorXd::Zero(n);
      for (int j = 0; j < p; ++j) s[pidx[j]] = z_P[j];

      bool all_pos = true;
      for (int j = 0; j < p; ++j)
        if (s[pidx[j]] <= tol) { all_pos = false; break; }

      if (all_pos) { x = s; break; }

      double alpha = std::numeric_limits<double>::infinity();
      for (int j = 0; j < p; ++j)
      {
        const int idx = pidx[j];
        if (s[idx] <= tol)
        {
          const double denom = x[idx] - s[idx];
          if (denom > 0.0)
            alpha = std::min(alpha, x[idx] / denom);
        }
      }
      if (!std::isfinite(alpha)) alpha = 0.0;

      x = x + alpha * (s - x);

      for (int i = 0; i < n; ++i)
        if (P[i] && x[i] <= tol)
        {
          x[i] = 0.0;
          P[i] = false;
          Z[i] = true;
        }
    }

    w = A.transpose() * (b - A * x);
  }

  // Ensure non-negativity.
  for (int i = 0; i < x.size(); ++i) if (x[i] < 0.0) x[i] = 0.0;
  return x;
}

// ---------------------------------------------------------------------------
// Iteratively reweighted least squares (IRLS) for sparse recovery.
// Minimises  ||x||_1  subject to  A x = b  (Candes-Wakin-Boyd 2008).
// This formulation satisfies the constraints exactly by construction.
// ---------------------------------------------------------------------------
static Eigen::VectorXd
irls_l1(const Eigen::MatrixXd& A,   // K x M
        const Eigen::VectorXd& b,   // K x 1
        int                    iters = 40,
        double                 eps   = 1e-12)
{
  const Eigen::Index M = A.cols();
  Eigen::VectorXd x = Eigen::VectorXd::Zero(M);

  for (int k = 0; k < iters; ++k)
  {
    // Diagonal weight matrix D = diag(|x_i| + eps).
    const Eigen::VectorXd d = (x.array().abs() + eps).matrix();

    // x = D A^T (A D A^T)^{-1} b.
    const Eigen::MatrixXd ADA = A * d.asDiagonal() * A.transpose(); // K x K
    const Eigen::VectorXd y   = ADA.ldlt().solve(b);                // K x 1
    x = d.asDiagonal() * (A.transpose() * y);                       // M x 1
  }
  return x;
}

} // anonymous namespace

// ---------------------------------------------------------------------------

std::vector<std::size_t>
compress(std::pair<std::vector<double>, std::vector<double>>& qr,
         std::size_t gdim,
         std::size_t quadrature_order,
         CompressionMethod method)
{
  switch (method)
  {
  case CompressionMethod::legacy: return compress_legacy(qr, gdim, quadrature_order);
  case CompressionMethod::nnls:   return compress_nnls  (qr, gdim, quadrature_order);
  case CompressionMethod::irls:   return compress_irls  (qr, gdim, quadrature_order);
  }
  throw std::runtime_error("Unknown compression method.");
}

// ---------------------------------------------------------------------------

std::vector<std::size_t>
compress_legacy(std::pair<std::vector<double>, std::vector<double>>& qr,
                std::size_t gdim,
                std::size_t quadrature_order)
{
  const std::size_t N    = quadrature_order;
  const std::size_t K    = choose(N + gdim, gdim);
  const std::size_t M    = qr.second.size();

  if (M <= K) return {};

  const auto qr_input = qr;

  // V is M x K.
  const Eigen::MatrixXd V = chebyshev_vandermonde(qr_input, gdim, N);

  // Economy-size QR: Q is M x K.
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(V);
  Eigen::MatrixXd Q = QR.householderQ();
  Q *= Eigen::MatrixXd::Identity(V.rows(), std::min(V.rows(), V.cols()));
  Q.transposeInPlace(); // Q is now K x M

  Eigen::Map<const Eigen::VectorXd> w_base(
      qr_input.second.data(), static_cast<Eigen::Index>(M));

  const Eigen::VectorXd nu    = Q * w_base;               // K x 1
  const Eigen::VectorXd w_new = Q.colPivHouseholderQr().solve(nu); // M x 1

  assert(w_new.size() == static_cast<Eigen::Index>(M));

  std::vector<std::size_t> idx;
  idx.reserve(M);
  for (Eigen::Index i = 0; i < w_new.size(); ++i)
    if (std::abs(w_new[i]) > 0.0)
      idx.push_back(static_cast<std::size_t>(i));

  assert(idx.size() <= K);

  qr.first.resize(gdim * idx.size());
  qr.second.resize(idx.size());
  for (std::size_t i = 0; i < idx.size(); ++i)
  {
    for (std::size_t d = 0; d < gdim; ++d)
      qr.first[gdim * i + d] = qr_input.first[gdim * idx[i] + d];
    qr.second[i] = w_new[static_cast<Eigen::Index>(idx[i])];
  }
  return idx;
}

// ---------------------------------------------------------------------------

std::vector<std::size_t>
compress_nnls(std::pair<std::vector<double>, std::vector<double>>& qr,
              std::size_t gdim,
              std::size_t quadrature_order)
{
  const std::size_t N = quadrature_order;
  const std::size_t K = choose(N + gdim, gdim);
  const std::size_t M = qr.second.size();

  if (M <= K) return {};

  const auto qr_input = qr;

  // V is M x K; A = V^T is K x M.
  const Eigen::MatrixXd V = chebyshev_vandermonde(qr_input, gdim, N);
  const Eigen::MatrixXd A = V.transpose(); // K x M

  Eigen::Map<const Eigen::VectorXd> w_orig(
      qr_input.second.data(), static_cast<Eigen::Index>(M));

  // b = A * w_orig is the K-vector of moments to reproduce.
  const Eigen::VectorXd b = A * w_orig;

  // Solve A w = b with w >= 0.
  const Eigen::VectorXd w_new = nnls_lawson_hanson(A, b, /*tol=*/1e-12);

  const double maxabs   = w_new.cwiseAbs().maxCoeff();
  const double tol_keep = std::max(1e-15, 1e-12 * maxabs);

  return apply_compressed_weights(qr, qr_input, w_new, gdim, K, tol_keep);
}

// ---------------------------------------------------------------------------

std::vector<std::size_t>
compress_irls(std::pair<std::vector<double>, std::vector<double>>& qr,
              std::size_t gdim,
              std::size_t quadrature_order)
{
  const std::size_t N = quadrature_order;
  const std::size_t K = choose(N + gdim, gdim);
  const std::size_t M = qr.second.size();

  if (M <= K) return {};

  const auto qr_input = qr;

  // V is M x K; A = V^T is K x M.
  const Eigen::MatrixXd V = chebyshev_vandermonde(qr_input, gdim, N);
  const Eigen::MatrixXd A = V.transpose(); // K x M

  Eigen::Map<const Eigen::VectorXd> w_orig(
      qr_input.second.data(), static_cast<Eigen::Index>(M));

  const Eigen::VectorXd b = A * w_orig;

  // IRLS: minimise ||w||_1 subject to A w = b.
  const Eigen::VectorXd w_irls =
      irls_l1(A, b, /*iters=*/40, /*eps=*/1e-12);

  // Select the support from the IRLS solution.
  const double maxabs   = w_irls.cwiseAbs().maxCoeff();
  const double tol_keep = std::max(1e-15, 1e-8 * maxabs);

  std::vector<std::size_t> idx;
  idx.reserve(static_cast<std::size_t>(w_irls.size()));
  for (Eigen::Index i = 0; i < w_irls.size(); ++i)
    if (std::abs(w_irls[i]) > tol_keep)
      idx.push_back(static_cast<std::size_t>(i));

  // Restrict to at most K points (keep those with the largest |w|).
  if (idx.size() > K)
  {
    std::vector<std::pair<double, std::size_t>> mag;
    mag.reserve(idx.size());
    for (auto j : idx)
      mag.emplace_back(
          std::abs(w_irls[static_cast<Eigen::Index>(j)]), j);
    std::nth_element(mag.begin(),
                     mag.begin() + static_cast<Eigen::Index>(K),
                     mag.end(),
                     [](const auto& a, const auto& b) {
                       return a.first > b.first;
                     });
    mag.resize(K);
    idx.clear();
    idx.reserve(K);
    for (const auto& p : mag) idx.push_back(p.second);
    std::sort(idx.begin(), idx.end());
  }

  // Re-solve the moment system restricted to the selected support to obtain
  // exact weights (moment-matching).
  Eigen::MatrixXd A_sub(static_cast<Eigen::Index>(K),
                         static_cast<Eigen::Index>(idx.size()));
  for (std::size_t j = 0; j < idx.size(); ++j)
    A_sub.col(static_cast<Eigen::Index>(j)) =
        A.col(static_cast<Eigen::Index>(idx[j]));
  const Eigen::VectorXd w_exact =
      A_sub.colPivHouseholderQr().solve(b);

  qr.first.resize(gdim * idx.size());
  qr.second.resize(idx.size());
  for (std::size_t i = 0; i < idx.size(); ++i)
  {
    for (std::size_t d = 0; d < gdim; ++d)
      qr.first[gdim * i + d] = qr_input.first[gdim * idx[i] + d];
    qr.second[i] = w_exact[static_cast<Eigen::Index>(i)];
  }
  return idx;
}

} // namespace simpex::compression
