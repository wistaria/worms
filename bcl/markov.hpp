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

#ifndef BCL_MARKOV_HPP
#define BCL_MARKOV_HPP

#include "random_choice.hpp"

namespace bcl {

template<typename CutoffType>
class markov_impl {
private:
  typedef random_choice<CutoffType> rc_type;
public:
  typedef std::size_t state_type;
  markov_impl() {}
  markov_impl(std::vector<std::vector<double> > const& tm) {
    init(tm);
  }
  template<typename MC, typename WVEC>
  markov_impl(MC const&, WVEC const& weights) {
    std::vector<std::vector<double> > tm;
    MC::generate_transition_matrix_resize(weights, tm);
    init(tm);
  }
  void init(std::vector<std::vector<double> > const& tm) {
    rc_.clear();
    for(std::size_t i = 0; i < tm.size(); ++i) rc_.push_back(rc_type(tm[i]));
  }
  template<typename RNG>
  state_type operator()(state_type prev, RNG& rng) {
    return rc_[prev](rng);
  }
private:
  std::size_t dim;
  std::vector<rc_type> rc_;
};

template<typename RNG>
class markov : public markov_impl<typename RNG::result_type> {
private:
  typedef markov_impl<typename RNG::result_type> base_type;
public:
  markov() : base_type() {}
  markov(std::vector<std::vector<double> > const& tm) : base_type(tm) {}
  template<typename MC, typename WVEC>
  markov(MC const& mc, WVEC const& weights) : base_type(mc, weights) {}
};

} // end namespace bcl

#endif // BCL_MARKOV_HPP
