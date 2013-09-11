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

#ifndef SPACETIME_POINT_HPP
#define SPACETIME_POINT_HPP

class spacetime_point {
public:
  spacetime_point() {}
  spacetime_point(int p, int n, int bop, int l) :
    prev_(p), next_(n), bop_(bop), leg_(l) {}
  int prev() const { return prev_; }
  int next() const { return next_; }
  int leg() const { return leg_; }
  int bond_operator() const { return bop_; }
  bool at_operator() const { return bop_ >= 0; }
  bool at_origin() const { return bop_ == -1; }
  bool at_staring() const { return bop_ == -2; }
  void set_prev(int p) { prev_ = p; }
  void set_next(int n) { next_ = n; }
  static spacetime_point origin(int s) { return spacetime_point(s, s, -1, 0); }
  static spacetime_point starting(int p, int n) { return spacetime_point(p, n, -2, 0); }
private:
  int prev_, next_;
  int bop_; // index of bond operator, -1 for t=0, -2 for worm starting point
  int leg_;
};

#endif // SPACETIME_POINT_HPP
