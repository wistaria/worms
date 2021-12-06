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

#include <iostream>
#include <gtest/gtest.h>
#include "worms/heisenberg_operator.hpp"
#include "worms/outgoing_weight.hpp"
#include "worms/weight.hpp"

TEST(outgoing_weight, test1) {
  for (double h = 0; h < 3; h += 0.1) {
    heisenberg_operator op(h, 2);
    weight w(op);
    outgoing_weight ow(w);
    std::cerr << "magnetic field = " << h << std::endl
              << "energy offset = " << w.offset() << std::endl;
    for (int p = 0; p < 16; ++p) {
      for (int g = 0; g < 4; ++g) {
        std::cerr << "  " << p << ' ' << g << ' ' << ow[p][g] << std::endl;
        EXPECT_TRUE((ow[p][g] >= 0));
      }
    }
  }
}
