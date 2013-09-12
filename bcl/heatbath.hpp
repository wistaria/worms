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

#ifndef BCL_HEATBATH_HPP
#define BCL_HEATBATH_HPP

#include <vector>
#include <algorithm>

namespace bcl {

class heatbath {
public:
  // Generate transition matrix of heat-bath update 
  template<class VEC, class MAT>
  static void generate_transition_matrix(VEC const& weights, MAT& tm) {
    typedef typename VEC::value_type value_type;
    value_type sum = std::accumulate(weights.begin(), weights.end(), value_type(0));
    int n = weights.size();
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j) tm[i][j] = weights[j] / sum;
  }
  template<class VEC, class MAT>
  static void generate_transition_matrix_resize(VEC const& weights, MAT& tm) {
    int n = weights.size();
    tm.resize(n);
    for (int i = 0; i < n; ++i) tm[i].resize(n);
    generate_transition_matrix(weights, tm);
  }
};
  
} // end namespace bcl

#endif // BCL_HEATBATH_HPP
