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

#include <worms/spin_state.hpp>
#include <vector>
#include <boost/foreach.hpp>

inline void check(bool test, std::string const& message) {
#ifdef CHECK_OPERATORS
  if (!test) {
    std::cerr << "Failed: " << message << std::endl;
    std::exit(127);
  }
#endif
}
  
template<typename LATTICE, typename BOND_OPERATOR, typename SPACETIME_POINT>
void check_operators(LATTICE const& lattice, std::vector<int> const& spins,
                     std::vector<BOND_OPERATOR> const& operators,
                     std::vector<SPACETIME_POINT> const& stpoints) {
#ifdef CHECK_OPERATORS
  check(lattice.num_sites() == spins.size(), "lattice.num_sites() == spins.size()");
  int n = lattice.num_sites();
  std::vector<int> current(spins);
  for (int b = 0; b < operators.size(); ++b) {
    BOND_OPERATOR bop = operators[b];
    int s0 = bop.site0();
    int s1 = bop.site1();
    check(s0 < n, "s0 < n");
    check(s1 < n, "s1 < n");
    check(current[s0] == spin_state::p2c(bop.state(), 0), "current[s0] == spin_state::p2c(bop.state(), 0)");
    check(current[s1] == spin_state::p2c(bop.state(), 1), "current[s1] == spin_state::p2c(bop.state(), 1)");
    check(bop.state() >= 0, "bop.state() >= 0");
    check(bop.state() < 16, "bop.state() < 16");
    current[s0] = spin_state::p2c(bop.state(), 2);
    current[s1] = spin_state::p2c(bop.state(), 3);
  }
  for (int s = 0; s < n; ++s)
    check(spins[s] == current[s], "spins[s] == current[s]");

  for (int b = 0; b < operators.size(); ++b) {
    BOND_OPERATOR bop = operators[b];
    int s0 = bop.site0();
    int s1 = bop.site1();
    int stp0 = bop.stp0();
    int stp1 = bop.stp1();
    check(stp0 < stpoints.size(), "stp0 < stpoints.size()");
    check(stp1 < stpoints.size(), "stp1 < stpoints.size()");
    check(stpoints[stp0].bond_operator() == b, "stpoints[stp0].bond_operator() == b");
    check(stpoints[stp1].bond_operator() == b, "stpoints[stp1].bond_operator() == b");
    check(stpoints[stp0].leg() == 0, "stpoints[stp0].leg() == 0");
    check(stpoints[stp1].leg() == 1, "stpoints[stp1].leg() == 1");
  }
  for (int s = 0; s < lattice.num_sites(); ++s) {
    check(stpoints[s].at_origin(), "stpoints[s].at_origin()");
  }
  for (int p = 0; p < stpoints.size(); ++p) {
    //// std::cerr << p << ' ' << stpoints[p].next() << ' ' << stpoints[stpoints[p].next()].prev() << std::endl;
    //// std::cerr << p << ' ' << stpoints[p].prev() << ' ' << stpoints[stpoints[p].prev()].next() << std::endl;
    check(stpoints[stpoints[p].next()].prev() == p, "stpoints[stpoints[p].next()].prev() == p");
    check(stpoints[stpoints[p].prev()].next() == p, "stpoints[stpoints[p].prev()].next() == p");
    if (stpoints[p].at_operator()) {
      int b = stpoints[p].bond_operator();
      BOND_OPERATOR bop = operators[b];
      int s0 = bop.site0();
      int s1 = bop.site1();
      int stp0 = bop.stp0();
      int stp1 = bop.stp1();
      check(stp0 < stpoints.size(), "stp0 < stpoints.size()");
      check(stp1 < stpoints.size(), "stp0 < stpoints.size()");
      check(stpoints[stp0].bond_operator() == b, "stpoints[stp0].bond_operator() == b");
      check(stpoints[stp1].bond_operator() == b, "stpoints[stp1].bond_operator() == b");
      check(stpoints[stp0].leg() == 0, "stpoints[stp0].leg() == 0");
      check(stpoints[stp1].leg() == 1, "stpoints[s1][stp1].leg() == 1");
    }
  }
#endif // CHECK_OPERATORS
}
