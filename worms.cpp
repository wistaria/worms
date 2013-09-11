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

#include "options.hpp"
#include "check.hpp"
#include "worms/bond_operator.hpp"
#include "worms/chain_lattice.hpp"
#include "worms/heisenberg_operator.hpp"
#include "worms/observable.hpp"
#include "worms/spacetime_point.hpp"
#include "worms/weight.hpp"

#include <bcl/markov.hpp>
#include <bcl/st2010.hpp>
#include <bcl/random_shuffle.hpp>

#include <vector>
#include <boost/random.hpp>
#include <boost/tuple/tuple.hpp>

template<typename BOND_OPERATOR, typename SPACETIME_POINT>
void insert_operator(int s0, int s1, int p, double t, std::vector<BOND_OPERATOR>& operators,
                     std::vector<std::vector<SPACETIME_POINT> >& stpoints) {
  int bindex = operators.size();
  int s0index = stpoints[s0].size();
  int s1index = stpoints[s1].size();
  operators.push_back(BOND_OPERATOR(s0, s1, s0index, s1index, p, t));
  stpoints[s0].push_back(SPACETIME_POINT(s0index-1, 0, bindex, 0));
  stpoints[s0][0].set_prev(s0index);
  stpoints[s0][s0index-1].set_next(s0index);
  stpoints[s1].push_back(SPACETIME_POINT(s1index-1, 0, bindex, 1));
  stpoints[s1][0].set_prev(s1index);
  stpoints[s1][s1index-1].set_next(s1index);
}

template<typename BOND_OPERATOR, typename SPACETIME_POINT>
void insert_operator(BOND_OPERATOR const& bop, std::vector<BOND_OPERATOR>& operators,
                     std::vector<std::vector<SPACETIME_POINT> >& stpoints) {
  insert_operator(bop.site0(), bop.site1(), bop.state(), bop.time(), operators, stpoints);
}

template<typename SPACETIME_POINT>
void insert_wstart(int s, double t, std::vector<std::vector<SPACETIME_POINT> >& stpoints,
                   std::vector<boost::tuple<int, int, double> >& wstart) {
  int sindex = stpoints[s].size();
  stpoints[s].push_back(SPACETIME_POINT::starting(sindex-1, 0));
  stpoints[s][0].set_prev(sindex);
  stpoints[s][sindex-1].set_next(sindex);
  wstart.push_back(boost::make_tuple(s, sindex, t));
}

int main(int argc, char* argv[]) {
  std::cout << "worms: a simple worm code\n";
  options opt(argc, argv, 16, 1.0);
  if (!opt.valid) std::exit(-1);
  double beta = 1 / opt.T;

  // lattice
  chain_lattice lattice(opt.L);
  
  // random number generator
  typedef boost::mt19937 engine_t;
  typedef boost::uniform_01<engine_t&> random01_t;
  typedef boost::exponential_distribution<> expdist_t;
  engine_t engine(29833u);
  random01_t random01(engine);

  // observables
  observable umag; // uniform magnetization
  observable umag2; // uniform magnetization^2
  observable smag2; // staggered magnetizetion^2

  // configuration
  std::vector<int> spins(opt.L, 0 /* all up */);
  std::vector<bond_operator> operators, operators_p;
  std::vector<std::vector<spacetime_point> > stpoints(opt.L);

  // Hamiltonian operator
  heisenberg_operator op(opt.H);
  double offset = 0; // >= 0

  // table for diagonal update
  double lambda = (op.max_diagonal() - op.min_diagonal() + offset);
  std::vector<double> accept(4);
  for (int c0 = 0; c0 < 2; ++c0)
    for (int c1 = 0; c1 < 2; ++c1)
      accept[spin_state::c2u(c0, c1)] = (op.max_diagonal() - op(c0, c1, c0, c1) + offset) / lambda;

  // table for worm update
  weight wt(op, offset);
  typedef bcl::markov<random01_t> markov_t;
  std::vector<markov_t> markov;
  for (int c = 0; c < 16; ++c) markov.push_back(markov_t(bcl::st2010(), wt[c]));

  // weight for worm insertion
  double wdensity = opt.L;
  std::vector<boost::tuple<int, int, double> > wstart;

  // worm statistics
  int wcount = 0;
  double wlength = 0;
  
  // temporaries
  std::vector<double> times;
  std::vector<int> current(opt.L);

  for (int mcs = 0; mcs < (opt.therm + opt.sweeps); ++mcs) {
    // diagonal update
    double pstart = wdensity / (opt.L * lambda + wdensity);
    expdist_t expdist(beta * (opt.L * lambda + wdensity));
    times.resize(0);
    double t = 0;
    while (t < 1) {
      t += expdist(random01);
      times.push_back(t);
    } // a sentinel (t >= 1) will be appended
    std::swap(operators, operators_p);
    operators.resize(0);
    operators_p.push_back(bond_operator(0, 0, 0, 0, 0, 1)); // sentinel
    std::copy(spins.begin(), spins.end(), current.begin());
    for (int s = 0; s < opt.L; ++s) {
      stpoints[s].resize(0);
      stpoints[s].push_back(spacetime_point::origin(s));
    }
    wstart.resize(0);
    check_operators(lattice, spins, operators);
    check_spacetime(lattice, operators, stpoints);
    std::vector<double>::iterator tmi = times.begin();
    for (std::vector<bond_operator>::iterator opi = operators_p.begin();
         opi != operators_p.end();) {
      if (*tmi < opi->time()) {
        if (random01() < pstart) {
          // insert worm starting point
          int s = static_cast<int>(opt.L * random01());
          insert_wstart(s, *tmi, stpoints, wstart);
        } else {
          // insert diagonal operator
          int b = static_cast<int>(opt.L * random01());
          int s0 = lattice.source(b);
          int s1 = lattice.target(b);
          int u = spin_state::c2u(current[s0], current[s1]);
          if (random01() < accept[u])
            insert_operator(s0, s1, spin_state::u2p(u, u), *tmi, operators, stpoints);
        }
        ++tmi;
      } else {
        if (opi->is_offdiagonal()) {
          // keep offdiagonal operator
          insert_operator(*opi, operators, stpoints);
          current[opi->site0()] = spin_state::p2c(opi->state(), 2);
          current[opi->site1()] = spin_state::p2c(opi->state(), 3);
        }
        ++opi;
      }
    }
    check_operators(lattice, spins, operators);
    check_spacetime(lattice, operators, stpoints);

    // shuffle starting point of worms
    bcl::random_shuffle(wstart.begin(), wstart.end(), random01);

    // worm update
    for (std::vector<boost::tuple<int, int, double> >::iterator wsi = wstart.begin();
         wsi != wstart.end(); ++wsi) {
      int direc = 2 * random01();
      int site, stp;
      double time;
      boost::tie(site, stp, time) = *wsi;
      int site_start = site;
      int stp_start = stp;
      wcount++;
      wlength += (direc == 0) ? time : -time;
      while (true) {
        stp = (direc == 0) ? stpoints[site][stp].prev() : stpoints[site][stp].next();
        if (stpoints[site][stp].at_operator()) {
          int bop = stpoints[site][stp].bond_operator();
          time = operators[bop].time();
          wlength += (direc == 0) ? -time : time;
          int ent = ((direc ^ 1) << 1) | stpoints[site][stp].leg();
          operators[bop].flip_state(ent);
          int ext = markov[operators[bop].state()](ent, random01);
          int leg = ext & 1;
          operators[bop].flip_state(ext);
          direc = (ext >> 1);
          site = (leg == 0) ? operators[bop].site0() : operators[bop].site1();
          stp = (leg == 0) ? operators[bop].stp0() : operators[bop].stp1();
          wlength += (direc == 0) ? time : -time;
        } else if (stpoints[site][stp].at_origin()) {
          spins[site] ^= 1;
          wlength += 1;
          time  = (direc == 0) ? 1 : 0;
        } else { // spoints[site][stp].at_starting()
          if (site == site_start && stp == stp_start) {
            wlength += (direc == 0) ? -wsi->get<2>() : wsi->get<2>();
            break;
          }
        }
      }
    }
    check_operators(lattice, spins, operators);
    check_spacetime(lattice, operators, stpoints);

    // measurement of physical quantities
    if (mcs >= opt.therm) {
      double mu = 0;
      double ms = 0;
      for (int s = 0; s < opt.L; ++s) {
        mu += 0.5 - spins[s];
        ms += lattice.phase(s) * (0.5 - spins[s]);
      }
      mu /= opt.L;
      ms /= opt.L;
      umag << mu;
      umag2 << mu * mu;
      smag2 << ms * ms;
    }
    
    if ((wcount > 0) && (mcs < opt.therm) && (mcs % (opt.therm / 4) == 0)) {
      std::cout << "Info: worm average number is reset from " << wdensity;
      wlength = wlength / wcount;
      wdensity = (beta * opt.L / wlength) + 1;
      wcount = 0;
      wlength = 0;
      std::cout << " to " << wdensity << std::endl;
    }
  }

  std::cout << "Uniform Magnetization     = "
            << umag.mean() << " +- " << umag.error() << std::endl
            << "Uniform Magnetization^2   = "
            << umag2.mean() << " +- " << umag2.error() << std::endl
            << "Staggered Magnetization^2 = "
            << smag2.mean() << " +- " << smag2.error() << std::endl;
}
