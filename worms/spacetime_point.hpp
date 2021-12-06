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
