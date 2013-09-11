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

#ifndef CHAIN_LATTICE_HPP
#define CHAIN_LATTICE_HPP

class chain_lattice {
public:
  chain_lattice(unsigned int L) : length_(L) {}
  unsigned int num_sites() const { return length_; }
  unsigned int num_bonds() const { return num_sites(); }
  unsigned int source(unsigned int b) const { return b; }
  unsigned int target(unsigned int b) const { return (b == length_-1) ? 0 : b+1; }
  double phase(unsigned int s) const { return (s & 1) ? 1.0 : -1.0; }
private:
  unsigned int length_;
};

#endif // CHAIN_LATTICE_HPP
