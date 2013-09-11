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

#ifndef BOND_OPERATOR_HPP
#define BOND_OPERATOR_HPP

#include "spin_state.hpp"

class bond_operator {
public:
  bond_operator() {}
  // generate diagonal operator
  bond_operator(int s0, int s1, int stp0, int stp1, int state, double t) :
    s0_(s0), s1_(s1), stp0_(stp0), stp1_(stp1), state_(state), time_(t) {}
  int site0() const { return s0_; }
  int site1() const { return s1_; }
  int stp0() const { return stp0_; }
  int stp1() const { return stp1_; }
  double time() const { return time_; }
  int state() const { return state_; }
  bool is_diagonal() const { return (spin_state::p2u(state_, 0) == spin_state::p2u(state_, 1)); }
  bool is_offdiagonal() const { return !is_diagonal(); }
  void flip_state(int leg) { state_ ^= spin_state::maskp(leg); }
  void print(std::ostream& os) const {
    os << s0_ << ' ' << s1_ << ' ' << stp0_ << ' ' << stp1_ << ' ' << state_ << ' ' << time_;
  }
private:
  int s0_, s1_; // site index of source and target
  int stp0_, stp1_; // stp = spacetime point
  int state_;
  double time_;
};

std::ostream& operator<<(std::ostream& os, bond_operator const& op) {
  op.print(os);
}
  
#endif // BOND_OPERATOR_HPP
