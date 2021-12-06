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

#include "spin_state.hpp"

class bond_operator {
public:
  typedef spin_state<2, 2> spin_state_t;
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
  bool is_diagonal() const { return (spin_state_t::p2u(state_, 0) == spin_state_t::p2u(state_, 1)); }
  bool is_offdiagonal() const { return !is_diagonal(); }
  void flip_state(int leg) { state_ ^= spin_state_t::maskp(leg); }
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
  return os;
}
