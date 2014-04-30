/*****************************************************************************
*
* worms: a simple worm code
*
* Copyright (C) 2013-2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef WEIGHT_HPP
#define WEIGHT_HPP

#include <algorithm>
#include <vector>
#include "outgoing_weight.hpp"

class weight {
public:
  template<typename OP>
  weight(OP const& op, double extra_offset = 0) {
    init_table(op, extra_offset);
  }
  template<typename OP>
  void init_table(OP const& op, double extra_offset = 0) {
    typedef typename OP::spin_state_t spin_state_t;
    offset_ = op.max_diagonal() + extra_offset;
    max_diagonal_weight_ = 0;
    weights_.resize(16);
    for (int p = 0; p < 16; ++p) {
      if (spin_state_t::is_diagonal(p)) {
        weights_[p] = -(op(p) - offset_);
        max_diagonal_weight_ = std::max(max_diagonal_weight_, weights_[p]);
      } else {
        weights_[p] = -op(p);
      }
    }
  }
  int size() const { return weights_.size(); }
  double offset() const { return offset_; }
  double max_diagonal_weight() const { return max_diagonal_weight_; }
  double operator[](int p) const { return weights_[p]; }
private:
  double offset_, max_diagonal_weight_;
  std::vector<double> weights_;
};

#endif // WEIGHT_HPP
