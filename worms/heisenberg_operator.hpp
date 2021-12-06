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
#include "spin_state.hpp"

class heisenberg_operator {
public:
  typedef spin_state<2, 2> spin_state_t;
  heisenberg_operator(double h, int coord_num) : elements_(16, 0) {
    elements_[spin_state_t::c2p(0, 0, 0, 0)] = 0.25 - h / coord_num;
    elements_[spin_state_t::c2p(0, 1, 0, 1)] = -0.25;
    elements_[spin_state_t::c2p(1, 0, 1, 0)] = -0.25;
    elements_[spin_state_t::c2p(1, 1, 1, 1)] = 0.25 + h / coord_num;
    elements_[spin_state_t::c2p(0, 1, 1, 0)] = -0.5;
    elements_[spin_state_t::c2p(1, 0, 0, 1)] = -0.5;
  }
  double operator()(int s) const { return elements_[s]; }
  double operator()(int sp0, int sp1, int sn0, int sn1) const {
    return elements_[spin_state_t::c2p(sp0, sp1, sn0, sn1)];
  }
  double max_diagonal() const {
    return std::max(std::max(elements_[spin_state_t::c2p(0, 0, 0, 0)],
                             elements_[spin_state_t::c2p(0, 1, 0, 1)]),
                    std::max(elements_[spin_state_t::c2p(1, 0, 1, 0)],
                             elements_[spin_state_t::c2p(1, 1, 1, 1)]));
  }
  double min_diagonal() const {
    return std::min(std::min(elements_[spin_state_t::c2p(0, 0, 0, 0)],
                             elements_[spin_state_t::c2p(0, 1, 0, 1)]),
                    std::min(elements_[spin_state_t::c2p(1, 0, 1, 0)],
                             elements_[spin_state_t::c2p(1, 1, 1, 1)]));
  }
private:
  std::vector<double> elements_;
};
