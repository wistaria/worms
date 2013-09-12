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

#ifndef BCL_ST2010_H
#define BCL_ST2010_H

#include <vector>
#include <algorithm>

namespace bcl {

class st2010 {
public:
  // Generate transition matrix
  template<class VEC, class MAT>
  static void generate_transition_matrix(VEC const& weights, MAT& tm) {
    using std::abs;
    typedef typename VEC::value_type value_type;
    std::size_t n = weights.size();
    std::vector<double> accum(n+1, 0);
    double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    double shift = *(std::max_element(weights.begin(), weights.end())) / sum;
    accum[0] = 0;
    for (std::size_t i = 0; i < n; ++i) accum[i+1] = accum[i] + weights[i] / sum;

    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        tm[i][j] = (std::max(std::min(accum[i+1] + shift, accum[j+1]) -
                             std::max(accum[i] + shift, accum[j]), 0.0) +
                    std::max(std::min(accum[i+1] + shift, accum[j+1] + 1) -
                             std::max(accum[i] + shift, accum[j] + 1), 0.0)) / (weights[i] / sum);
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
  static std::size_t choose_next(VEC const& weights, std::size_t present, RNG& rng){
    double sum = *(std::max_element(weights.begin(), weights.end()));
    sum -= weights[present] * rng();
    for (int i = 0; i < weights.size(); ++i) {
      int j = (present + i + 1) % weights.size();
      if (sum <= weights[j]) return j;
      sum -= weights[j];
    }
    return present;
  }
};

} // end namespace bcl

#endif // BCL_ST2010_H
