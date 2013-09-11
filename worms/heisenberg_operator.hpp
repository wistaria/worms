/*****************************************************************************
*
* worms: a simple worm code
*
* Copyright (C) 2013 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef HEISENBERG_OPERATOR_HPP
#define HEISENBERG_OPERATOR_HPP

#include <algorithm>
#include <vector>
#include "spin_state.hpp"

class heisenberg_operator {
public:
  heisenberg_operator(double h) : elements_(16, 0) {
    elements_[spin_state::c2p(0, 0, 0, 0)] = 0.25 - h;
    elements_[spin_state::c2p(0, 1, 0, 1)] = -0.25;
    elements_[spin_state::c2p(1, 0, 1, 0)] = -0.25;
    elements_[spin_state::c2p(1, 1, 1, 1)] = 0.25 + h;
    elements_[spin_state::c2p(0, 1, 1, 0)] = -0.5;
    elements_[spin_state::c2p(1, 0, 0, 1)] = -0.5;
  }
  double operator()(int s) const { return elements_[s]; }
  double operator()(int sp0, int sp1, int sn0, int sn1) const {
    return elements_[spin_state::c2p(sp0, sp1, sn0, sn1)];
  }
  double max_diagonal() const {
    return std::max(std::max(elements_[spin_state::c2p(0, 0, 0, 0)],
                             elements_[spin_state::c2p(0, 1, 0, 1)]),
                    std::max(elements_[spin_state::c2p(1, 0, 1, 0)],
                             elements_[spin_state::c2p(1, 1, 1, 1)]));
  }
  double min_diagonal() const {
    return std::min(std::min(elements_[spin_state::c2p(0, 0, 0, 0)],
                             elements_[spin_state::c2p(0, 1, 0, 1)]),
                    std::min(elements_[spin_state::c2p(1, 0, 1, 0)],
                             elements_[spin_state::c2p(1, 1, 1, 1)]));
  }
private:
  std::vector<double> elements_;
};

#endif // HEISENBERG_OPERATOR_HPP
