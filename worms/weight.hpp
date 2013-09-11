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

#ifndef WEIGHT_HPP
#define WEIGHT_HPP

#include <algorithm>
#include <vector>

class weight {
public:
  template<typename OP>
  weight(OP const& op, double offset = 0) { init_table(op, offset); }
  template<typename OP>
  void init_table(OP const& op, double offset = 0) {
    double c = op.max_diagonal() + offset;
    weights_.clear();
    weights_.resize(16);
    for (int s = 0; s < 16; ++s) {
      weights_[s].resize(4);
      for (int g = 0; g < 4; ++g) {
        weights_[s][g] = -(op(s ^ (1<<g)) - c);
        if (weights_[s][g] < 0)
          std::cerr << "Warning: negative weights (" << weights_[s][g] << ") for (s,g)=(" << s
                    << "," << g << ")\n";
      }
    }
  }
  std::vector<double> const& operator[](int s) { return weights_[s]; }
private:
  std::vector<std::vector<double> > weights_;
};

#endif // WEIGHT_HPP
