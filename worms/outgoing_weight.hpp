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

#ifndef OUTGOING_WEIGHT_HPP
#define OUTGOING_WEIGHT_HPP

#include <algorithm>
#include <vector>
#include "weight.hpp"

class outgoing_weight {
public:
  outgoing_weight(weight const& w) { init_table(w); }
  void init_table(weight const& w) {
    weights_.clear();
    weights_.resize(16);
    for (int s = 0; s < 16; ++s) {
      weights_[s].resize(4);
      for (int g = 0; g < 4; ++g) weights_[s][g] = w[s ^ (1<<g)];
    }
  }
  std::vector<double> const& operator[](int s) { return weights_[s]; }
private:
  std::vector<std::vector<double> > weights_;
};

#endif // OUTGOING_WEIGHT_HPP
