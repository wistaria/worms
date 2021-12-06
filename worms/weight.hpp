/*
   worms: a simple worm code

   Copyright (C) 2013-2021 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#pragma once

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
    weights_.resize(spin_state_t::num_configurations);
    for (std::size_t p = 0; p < weights_.size(); ++p) {
      if (spin_state_t::is_diagonal(p)) {
        weights_[p] = -(op(p) - offset_);
        max_diagonal_weight_ = std::max(max_diagonal_weight_, weights_[p]);
      } else {
        weights_[p] = -op(p);
      }
    }
  }
  int num_configurations() const { return weights_.size(); }
  double offset() const { return offset_; }
  double max_diagonal_weight() const { return max_diagonal_weight_; }
  double operator[](int p) const { return weights_[p]; }
private:
  double offset_, max_diagonal_weight_;
  std::vector<double> weights_;
};
