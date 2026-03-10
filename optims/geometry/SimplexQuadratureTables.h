#pragma once
#include <cstddef>
#include <vector>

namespace simpex::tables {

  struct RuleRef
  {
    std::vector<double> pts; // packed as (r, s, t) triples; unused entries are 0
    std::vector<double> w;
  };

  // Interval rules on [0,1]; points stored as (r,0,0), sum(w) = 1.
  RuleRef interval_rule_legendre(std::size_t order);

  // Triangle rules on the reference triangle (0,0),(1,0),(0,1); sum(w) = 0.5.
  // Supported orders: 1..20.
  RuleRef triangle_rule_dunavant(std::size_t order);

  // Tetrahedron rules on the reference tetrahedron; sum(w) = 1/6.
  // Supported orders: 1..6.
  RuleRef tetra_rule_keast(std::size_t order);

} // namespace simpex::tables
