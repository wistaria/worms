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

#include <vector>

class square_lattice {
public:
  square_lattice(unsigned int L) : length_x_(L), length_y_(L) { init(); }
  square_lattice(unsigned int Lx, unsigned int Ly) : length_x_(Lx), length_y_(Ly) { init(); }
  void init() {
    source_.resize(2 * length_x_ * length_y_);
    target_.resize(2 * length_x_ * length_y_);
    site_phase_.resize(length_x_ * length_y_);
    unsigned int n = length_x_ * length_y_ * 2;
    for (unsigned int s = 0; s < length_x_ * length_y_; ++s) {
      int ix = s % length_x_;
      int iy = s / length_y_;
      site_phase_[s] = 2 * ((ix + iy) % 2) - 1;
    }
    for (unsigned int b = 0; b < n; ++b) {
      unsigned int s = b / 2;
      unsigned int t;
      if (b % 2 == 0) {
        t = (s + 1) % length_x_ + (s / length_x_) * length_x_; // target right
      } else {
        t = (s + length_x_) % (length_x_ * length_y_); // target below
      }
      source_[b] = s;
      target_[b] = t;
    }
  }
  unsigned int coordination_num() const { return 4; }
  unsigned int num_sites() const { return length_x_ * length_y_; }
  unsigned int num_bonds() const { return 2 * num_sites(); }
  unsigned int source(unsigned int b) const { return source_[b]; }
  unsigned int target(unsigned int b) const { return target_[b]; }
  double phase(unsigned int s) const { return site_phase_[s]; }
private:
  unsigned int length_x_, length_y_;
  std::vector<unsigned int> source_;
  std::vector<unsigned int> target_;
  std::vector<double> site_phase_;
};
