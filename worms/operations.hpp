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

#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include "bond_operator.hpp"
#include "spacetime_point.hpp"
#include <boost/tuple/tuple.hpp>

template <int Spin>
inline void append_operator(int s0, int s1, int p, double t, std::vector<bond_operator<Spin> >& operators,
                     std::vector<spacetime_point>& stpoints) {
  int bindex = operators.size();
  int s0index = stpoints.size();
  int s1index = s0index + 1;
  operators.push_back(bond_operator<Spin>(s0, s1, s0index, s1index, p, t));
  stpoints.push_back(spacetime_point(stpoints[s0].prev(), s0, bindex, 0));
  stpoints[stpoints[s0].prev()].set_next(s0index);
  stpoints[s0].set_prev(s0index);
  stpoints.push_back(spacetime_point(stpoints[s1].prev(), s1, bindex, 1));
  stpoints[stpoints[s1].prev()].set_next(s1index);
  stpoints[s1].set_prev(s1index);
}

template <int Spin>
inline void append_operator(bond_operator<Spin> const& bop, std::vector<bond_operator<Spin> >& operators,
			    std::vector<spacetime_point>& stpoints) {
  append_operator(bop.site0(), bop.site1(), bop.state(), bop.time(), operators, stpoints);
}

inline void append_wstart(int s, double t, std::vector<spacetime_point>& stpoints,
                   std::vector<boost::tuple<int, int, double> >& wstart) {
  int sindex = stpoints.size();
  stpoints.push_back(spacetime_point::starting(stpoints[s].prev(), s));
  stpoints[stpoints[s].prev()].set_next(sindex);
  stpoints[s].set_prev(sindex);
  wstart.push_back(boost::make_tuple(s, sindex, t));
}

#ifdef CHECK_OPERATORS
inline void check(bool test, std::string const& message) {
  if (!test) {
    std::cerr << "Failed: " << message << std::endl;
    std::exit(127);
  }
}
#endif
  
template<typename LATTICE, int Spin>
void check_operators(LATTICE const& lattice, std::vector<int> const& spins,
                     std::vector<bond_operator<Spin> > const& operators,
                     std::vector<spacetime_point> const& stpoints) {
#ifdef CHECK_OPERATORS
  check(lattice.num_sites() == spins.size(), "lattice.num_sites() == spins.size()");
  int n = lattice.num_sites();
  std::vector<int> current(spins);
  for (int b = 0; b < operators.size(); ++b) {
    bond_operator<Spin> bop = operators[b];
    int s0 = bop.site0();
    int s1 = bop.site1();
    check(s0 < n, "s0 < n");
    check(s1 < n, "s1 < n");
    check(current[s0] == spin_state::p2c<Spin>(bop.state(), 0), "current[s0] == spin_state::p2c(bop.state(), 0)");
    check(current[s1] == spin_state::p2c<Spin>(bop.state(), 1), "current[s1] == spin_state::p2c(bop.state(), 1)");
    check(bop.state() >= 0, "bop.state() >= 0");
    check(bop.state() < operatorsize<Spin>::val, "bop.state() < 16");
    current[s0] = spin_state::p2c<Spin>(bop.state(), 2);
    current[s1] = spin_state::p2c<Spin>(bop.state(), 3);
  }
  for (int s = 0; s < n; ++s)
    check(spins[s] == current[s], "spins[s] == current[s]");

  for (int b = 0; b < operators.size(); ++b) {
    bond_operator<Spin> bop = operators[b];
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
      bond_operator<Spin> bop = operators[b];
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

#endif // OPERATIONS_HPP
