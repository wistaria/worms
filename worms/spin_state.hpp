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

#ifndef SPIN_STATE_HPP
#define SPIN_STATE_HPP

// p: spin state of operator (4 bits)
// u: spin state of site pair (2 bits)
// c: spin state of a site (1 bit)

// [u <-> c]
//   0-1

// [p <-> u]
//   1
//   .
//   0

// [p <-> c]
//   2 3
//   .-.
//   0 1

namespace spin_state {
  inline int p2c(int p, int l) { return (p >> l) & 1; }
  inline int p2u(int p, int d) { return (p >> (2 * d)) & 3; }
  inline int c2u(int c0, int c1) { return (c0 | (c1 << 1)); }
  inline int c2p(int c0, int c1, int c2, int c3) {
    return (c0 | (c1 << 1) | (c2 << 2) | (c3 << 3));
  }
  inline int u2p(int u0, int u1) { return (u0 | (u1 << 2)); }
  inline int maskp(int l) { return (1 << l); }
};

#endif // SPIN_STATE_HPP
