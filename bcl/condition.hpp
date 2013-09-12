/*****************************************************************************
*
* BCL: Balance Condition Library
*
* Copyright (C) 2009-2012 by Hidemaro Suwa <suwamaro@looper.t.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>

* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef BCL_CONDITION_HPP
#define BCL_CONDITION_HPP

#include <cmath>
#include <numeric>

namespace bcl {

template<class MAT>
bool check_probability_conservation(MAT const& transition_matrix, double tolerance = 1.0e-10) {
  std::size_t n = transition_matrix.size();
  for (std::size_t i = 0; i < n; ++i) {
    double p = 0;
    for (std::size_t j = 0; j < n; ++j) p += transition_matrix[i][j];
    if (std::abs(p - 1) > tolerance) return false;
  }
  return true;
}

template<class VEC, class MAT>
bool check_balance_condition(VEC const& weights, MAT const& transition_matrix,
  double tolerance = 1.0e-10) {
  typedef typename VEC::value_type value_type;
  std::size_t n = weights.size();
  for (std::size_t j = 0; j < n; ++j) {
    value_type w = 0;
    for (std::size_t i = 0; i < n; ++i) w += weights[i] * transition_matrix[i][j];
    if (std::abs(w - weights[j]) > tolerance * weights[j]) return false;
  }
  return true;
}

template<class VEC, class MAT>
bool check_detailed_balance(VEC const& weights, MAT const& transition_matrix,
  double tolerance = 1.0e-10) {
  typedef typename VEC::value_type value_type;
  value_type sum = std::accumulate(weights.begin(), weights.end(), value_type(0));
  std::size_t n = weights.size();
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      if (std::abs(weights[i] * transition_matrix[i][j] - weights[j] * transition_matrix[j][i]) >
          tolerance * sum) return false;
  return true;
}

template<class VEC, class MAT>
bool check(VEC const& weights, MAT const& transition_matrix, double tolerance = 1.0e-10) {
  return check_probability_conservagion(transition_matrix, tolerance) &&
    check_balance_condition(weights, transition_matrix, tolerance);
}

template<class VEC, class MAT>
typename VEC::value_type average_rejection(VEC const& weights, MAT const& transition_matrix) {
  typename VEC::value_type vs = 0;
  typename VEC::value_type ws = 0;
  std::size_t n = weights.size();
  for (std::size_t i = 0; i < n; ++i) {
    vs = weights[i] * transition_matrix[i][i];
    ws = weights[i];
  }
  return vs / ws;
}

} // end namespace bcl

#endif // BCL_CONDITION_HPP
