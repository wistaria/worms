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

#define BOOST_TEST_MODULE test_outgoing_weight
#ifndef BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#else
#include <boost/test/unit_test.hpp>
#endif

#include "worms/heisenberg_operator.hpp"
#include "worms/outgoing_weight.hpp"

#define SPIN 2

BOOST_AUTO_TEST_CASE(test_outgoing_weight) {
  for (double h = 0; h < 3; h += 0.1) {
    heisenberg_operator<SPIN> op(h, 2);
    weight<SPIN> w(op);
    outgoing_weight ow(w);
    std::cerr << "magnetic field = " << h << std::endl
              << "energy offset = " << w.offset() << std::endl;
    for (int p = 0; p < operatorsize<SPIN>::val; ++p) {
      for (int g = 0; g < 4; ++g) {
        std::cerr << "  " << p << ' ' << g << ' ' << ow[p][g] << std::endl;
        BOOST_CHECK((ow[p][g] >= 0));
      }
    }
  }
}
