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
  if (!test) {
    std::cerr << "Failed: " << message << std::endl;
    std::exit(127);
  }
}
  
template<typename LATTICE, typename BOND_OPERATOR>
void check_operators(LATTICE const& lattice, std::vector<int> const& spins,
                 std::vector<BOND_OPERATOR> const& operators) {
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
#endif // CHECK_OPERATORS
}

template<typename LATTICE, typename BOND_OPERATOR, typename SPACETIME_POINT>
void check_spacetime(LATTICE const& lattice, std::vector<BOND_OPERATOR> const& operators,
                     std::vector<std::vector<SPACETIME_POINT> > const& stpoints) {
#ifdef CHECK_SPACETIME
  for (int b = 0; b < operators.size(); ++b) {
    BOND_OPERATOR bop = operators[b];
    int s0 = bop.site0();
    int s1 = bop.site1();
    int stp0 = bop.stp0();
    int stp1 = bop.stp1();
    check(stp0 < stpoints[s0].size(), "stp0 < stpoints[s0].size()");
    check(stp1 < stpoints[s1].size(), "stp0 < stpoints[s0].size()");
    check(stpoints[s0][stp0].bond_operator() == b, "stpoints[s0][stp0].bond_operator() == b");
    check(stpoints[s1][stp1].bond_operator() == b, "stpoints[s1][stp1].bond_operator() == b");
    check(stpoints[s0][stp0].leg() == 0, "stpoints[s0][stp0].leg() == 0");
    check(stpoints[s1][stp1].leg() == 1, "stpoints[s1][stp1].leg() == 1");
  }
  for (int s = 0; s < lattice.num_sites(); ++s) {
    check(stpoints[s][0].at_origin(), "stpoints[s][0].at_origin()");
    for (int p = 0; p < stpoints[s].size(); ++p) {
      check(stpoints[s][stpoints[s][p].next()].prev() == p, "stpoints[s][stpoints[s][p].next()].prev() == p");
      check(stpoints[s][stpoints[s][p].prev()].next() == p, "stpoints[s][stpoints[s][p].prev()].next() == p");
      if (stpoints[s][0].at_operator()) {
        int b = stpoints[s][0].bond_operator();
        BOND_OPERATOR bop = operators[b];
        int s0 = bop.site0();
        int s1 = bop.site1();
        int stp0 = bop.stp0();
        int stp1 = bop.stp1();
        check(stp0 < stpoints[s0].size(), "stp0 < stpoints[s0].size()");
        check(stp1 < stpoints[s1].size(), "stp0 < stpoints[s0].size()");
        check(stpoints[s0][stp0].bond_operator() == b, "stpoints[s0][stp0].bond_operator() == b");
        check(stpoints[s1][stp1].bond_operator() == b, "stpoints[s1][stp1].bond_operator() == b");
        check(stpoints[s0][stp0].leg() == 0, "stpoints[s0][stp0].leg() == 0");
        check(stpoints[s1][stp1].leg() == 1, "stpoints[s1][stp1].leg() == 1");
      }
    }
  }
#endif // CHECK_SPACETIME
}
