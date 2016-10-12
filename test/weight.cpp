/*****************************************************************************
*
* worms: a simple worm code
*
* Copyright (C) 2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#define BOOST_TEST_MODULE test_weight
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "worms/heisenberg_operator.hpp"
#include "worms/weight.hpp"

#define SPIN 2

BOOST_AUTO_TEST_CASE(test_weight) {
  for (double h = 0; h < 3; h += 0.1) {
    heisenberg_operator<SPIN> op(h, 2);
    weight<SPIN> w(op);
    std::cerr << "magnetic field = " << h << std::endl
              << "energy offset = " << w.offset() << std::endl;
    for (int p = 0; p < operatorsize<SPIN>::val; ++p) {
      std::cerr << "  " << p << ' ' << w[p] << std::endl;
      BOOST_CHECK(w[p] >= 0);
    }
  }
  for (double h = 0; h < 3; h += 0.1) {
    heisenberg_operator<SPIN> op(h, 2);
    weight<SPIN> w(op, 0.5);
    std::cerr << "magnetic field = " << h << std::endl
              << "energy offset = " << w.offset() << std::endl;
    for (int p = 0; p <  operatorsize<SPIN>::val; ++p) {
      std::cerr << "  " << p << ' ' << w[p] << std::endl;
      BOOST_CHECK(w[p] >= 0);
    }
  }
}
