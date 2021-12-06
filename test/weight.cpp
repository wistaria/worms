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
#include "worms/weight.hpp"

TEST(weight, test1) {
  for (double h = 0; h < 3; h += 0.1) {
    heisenberg_operator op(h, 2);
    weight w(op);
    std::cerr << "magnetic field = " << h << std::endl
              << "energy offset = " << w.offset() << std::endl;
    for (int p = 0; p < 16; ++p) {
      std::cerr << "  " << p << ' ' << w[p] << std::endl;
      EXPECT_TRUE(w[p] >= 0);
    }
  }
}
TEST(weight, test2) {
  for (double h = 0; h < 3; h += 0.1) {
    heisenberg_operator op(h, 2);
    weight w(op, 0.5);
    std::cerr << "magnetic field = " << h << std::endl
              << "energy offset = " << w.offset() << std::endl;
    for (int p = 0; p < 16; ++p) {
      std::cerr << "  " << p << ' ' << w[p] << std::endl;
      EXPECT_TRUE(w[p] >= 0);
    }
  }
}
