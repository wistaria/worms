/*****************************************************************************
*
* BCL: Balance Condition Library
*
* Copyright (C) 2009-2012 by Hidemaro Suwa <suwamaro@looper.t.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef BCL_METROPOLIS_HPP
#define BCL_METROPOLIS_HPP

#include <vector>
#include <algorithm>
#include <cmath>

namespace bcl {

class metropolis {
public:
  // Generate transition matrix
  template<class VEC, class MAT>
  static void generate_transition_matrix(VEC const& weights, MAT& tm) {
    using std::abs;
    typedef typename VEC::value_type value_type;
    std::size_t n = weights.size();
    for (std::size_t i = 0; i < n; ++i) {
      tm[i][i] = 1;
      if (weights[i] > 0) {
        for (std::size_t j = 0; j < n; ++j) {
          if (i != j) {
            tm[i][j] = std::min(weights[i], weights[j]) / weights[i] / (n-1);
            tm[i][i] -= tm[i][j];
          }
        }
        tm[i][i] = abs(tm[i][i]);
      } else {
        for (std::size_t j = 0; j < n; ++j)
          if (i != j) tm[i][j] = value_type(0);
      }
    }
  }
  template<class VEC, class MAT>
  static void generate_transition_matrix_resize(VEC const& weights, MAT& tm) {
    std::size_t n = weights.size();
    tm.resize(n);
    for (std::size_t i = 0; i < n; ++i) tm[i].resize(n);
    generate_transition_matrix(weights, tm);
  }

  template<class VEC, class RNG>
  static std::size_t choose_next(VEC const& weights, std::size_t present, RNG& rng) {
    std::size_t proposal = weights.size() * rng();
    return (weights(present) * rng() < weights(proposal)) ? proposal : present;
  }
};

} // end namespace bcl

#endif // BCL_METROPOLIS_HPP
