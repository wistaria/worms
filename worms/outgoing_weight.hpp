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

class outgoing_weight {
public:
  template<typename WEIGHT>
  outgoing_weight(WEIGHT const& w) { init_table(w); }
  template<typename WEIGHT>
  void init_table(WEIGHT const& w) {
    weights_.clear();
    weights_.resize(w.num_configurations());
    for (int s = 0; s < w.num_configurations(); ++s) {
      weights_[s].resize(4);
      for (int g = 0; g < 4; ++g) weights_[s][g] = w[s ^ (1<<g)];
    }
  }
  std::vector<double> const& operator[](int s) { return weights_[s]; }
private:
  std::vector<std::vector<double> > weights_;
};
