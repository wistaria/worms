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

template <int Spin>
class outgoing_weight_raise {
public:
  outgoing_weight_raise(weight<Spin> const& w) { init_table(w); }
  void init_table(weight<Spin> const& w) {
    weights_.clear();
    weights_.resize(operatorsize<Spin>::val);
    for (int s = 0; s < operatorsize<Spin>::val; ++s) {
      weights_[s].resize(4);
      for (int g = 0; g < 4; ++g){
	int state = s;
	if (spin_state::raise<Spin>(g, &state))
	  weights_[s][g] = w[state];
	else
	  weights_[s][g] = 0;

	}
    }
  }
  std::vector<double> const& operator[](int s) { return weights_[s]; }
private:
  std::vector<std::vector<double> > weights_;
};

template <int Spin>
class outgoing_weight_lower {
public:
  outgoing_weight_lower(weight<Spin> const& w) { init_table(w); }
  void init_table(weight<Spin> const& w) {
    weights_.clear();
    weights_.resize(operatorsize<Spin>::val);
    for (int s = 0; s < operatorsize<Spin>::val; ++s) {
      weights_[s].resize(4);
      for (int g = 0; g < 4; ++g){
	  int state = s;
	  if (spin_state::lower<Spin>(g, &state))
	    weights_[s][g] = w[state];
	  else
	    weights_[s][g] = 0;
	}
    }
  }
  std::vector<double> const& operator[](int s) { return weights_[s]; }
private:
  std::vector<std::vector<double> > weights_;
};


class outgoing_weight {
public:
  outgoing_weight(weight<2> const& w) { init_table(w); }
  void init_table(weight<2> const& w) {
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
