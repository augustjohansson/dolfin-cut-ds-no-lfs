#pragma once

#include <cstddef>
#include <utility>
#include <vector>

/// Quadrature rule compression algorithms.
///
/// References:
///   A. Sommariva, M. Vianello,
///   "Compression of multivariate discrete measures and applications",
///   Numerical Functional Analysis and Optimization, 36(9), 2015.
///
///   A. Sommariva, M. Vianello, "CheapQ: fast positive quadrature via compressed
///   Chebyshev Vandermonde systems", https://github.com/alvisesommariva/CheapQ.
///
///   G. Migliorati, T. Pruneri, F. Nobile,
///   "Tetrahedral-free quadrature on simplices via moment fitting (Tetrafreeq)",
///   Journal of Computational Physics 493 (2023).
///   https://www.sciencedirect.com/science/article/abs/pii/S016892742300185X

namespace dolfin::compression
{

/// Compression algorithm to use.
enum class CompressionMethod
{
  /// QR-based legacy method (signed weights, matches the original DOLFIN
  /// implementation). Analogous to the algorithm in Sommariva–Vianello 2015.
  legacy,
  /// NNLS-based positive-weight compression (Tchakaloff / Tetrafreeq style).
  /// All compressed weights are non-negative.
  nnls,
  /// IRLS-based sparse compression (CheapQ-compatible, signed weights).
  irls
};

/// Compress a quadrature rule using the selected algorithm.
///
/// \param[in,out] qr               Quadrature rule to compress.
/// \param[in]     gdim             Geometric dimension.
/// \param[in]     quadrature_order Polynomial degree integrated exactly.
/// \param[in]     method           Compression algorithm.
/// \returns Indices of retained points, or an empty vector if no compression
///          was performed.
std::vector<std::size_t>
compress(std::pair<std::vector<double>, std::vector<double>>& qr,
         std::size_t gdim,
         std::size_t quadrature_order,
         CompressionMethod method = CompressionMethod::legacy);

/// QR-based legacy compression (may produce negative weights).
std::vector<std::size_t>
compress_legacy(std::pair<std::vector<double>, std::vector<double>>& qr,
                std::size_t gdim,
                std::size_t quadrature_order);

/// NNLS-based compression yielding non-negative weights (Tchakaloff / Tetrafreeq).
std::vector<std::size_t>
compress_nnls(std::pair<std::vector<double>, std::vector<double>>& qr,
              std::size_t gdim,
              std::size_t quadrature_order);

/// IRLS-based sparse compression (signed weights, CheapQ-compatible).
std::vector<std::size_t>
compress_irls(std::pair<std::vector<double>, std::vector<double>>& qr,
              std::size_t gdim,
              std::size_t quadrature_order);

} // namespace dolfin::compression
