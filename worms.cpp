/*****************************************************************************
*
* worms: a simple worm code
*
* Copyright (C) 2013-2014 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef NDEBUG
# define CHECK_OPERATORS 1
#endif

#include "options.hpp"
#include "worms/version.hpp"
#include "worms/chain_lattice.hpp"
#include "worms/heisenberg_operator.hpp"
#include "worms/operations.hpp"
#include "worms/weight.hpp"
#include "worms/outgoing_weight.hpp"

#include <bcl.hpp>

#include <vector>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/tuple/tuple.hpp>

int main(int argc, char* argv[]) {
  std::cerr << "worms: a simple worm code (release " WORMS_VERSION ")\n"
            << "  Copyright (C) 2013-2014 by Synge Todo <wistaria@comp-phys.org>\n"
            << "  " WORMS_URL "\n\n";
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
  bcl::observable ene; // energy density
  bcl::observable umag; // uniform magnetization
  bcl::observable umag2; // uniform magnetization^2
  bcl::observable smag2; // staggered magnetizetion^2

  // configuration
  std::vector<int> spins(opt.L, 0 /* all up */);
  std::vector<bond_operator> operators, operators_p;
  std::vector<spacetime_point> stpoints;

  // Hamiltonian operator
  heisenberg_operator op(opt.H, /* coordination number = */ 2);
  typedef heisenberg_operator::spin_state_t spin_state_t;
  double offset = 0; // >= 0
  weight wt(op, offset);

  // table for diagonal update
  double lambda = wt.max_diagonal_weight();
  std::vector<double> accept(4);
  for (int c0 = 0; c0 < 2; ++c0)
    for (int c1 = 0; c1 < 2; ++c1)
      accept[spin_state_t::c2u(c0, c1)] = wt[spin_state_t::c2p(c0, c1, c0, c1)] / lambda;

  // table for worm update
  outgoing_weight ogwt(wt);
  typedef bcl::markov<random01_t> markov_t;
  std::vector<markov_t> markov;
  for (int c = 0; c < 16; ++c) markov.push_back(markov_t(bcl::st2010(), ogwt[c]));

  // weight for worm insertion
  double wdensity = opt.L;
  std::vector<boost::tuple<int, int, double> > wstart;

  // worm statistics
  double wcount = 0;
  double wlength = 0;
  
  // temporaries
  std::vector<double> times;
  std::vector<int> current(opt.L);

  boost::timer tm;
  for (int mcs = 0; mcs < (opt.therm + opt.sweeps); ++mcs) {
    // diagonal update
    double pstart = wdensity / (beta * opt.L * lambda + wdensity);
    expdist_t expdist(beta * opt.L * lambda + wdensity);
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
    stpoints.resize(0);
    for (int s = 0; s < opt.L; ++s) stpoints.push_back(spacetime_point::origin(s));
    wstart.resize(0);
    check_operators(lattice, spins, operators, stpoints);
    std::vector<double>::iterator tmi = times.begin();
    for (std::vector<bond_operator>::iterator opi = operators_p.begin();
         opi != operators_p.end();) {
      if (*tmi < opi->time()) {
        if (random01() < pstart) {
          // insert worm starting point
          int s = static_cast<int>(opt.L * random01());
          append_wstart(s, *tmi, stpoints, wstart);
        } else {
          // insert diagonal operator
          int b = static_cast<int>(opt.L * random01());
          int s0 = lattice.source(b);
          int s1 = lattice.target(b);
          int u = spin_state_t::c2u(current[s0], current[s1]);
          if (random01() < accept[u])
            append_operator(s0, s1, spin_state_t::u2p(u, u), *tmi, operators, stpoints);
        }
        ++tmi;
      } else {
        if (opi->is_offdiagonal()) {
          // keep offdiagonal operator
          append_operator(*opi, operators, stpoints);
          current[opi->site0()] = spin_state_t::p2c(opi->state(), 2);
          current[opi->site1()] = spin_state_t::p2c(opi->state(), 3);
        }
        ++opi;
      }
    }
    check_operators(lattice, spins, operators, stpoints);

    // shuffle starting point of worms
    bcl::random_shuffle(wstart.begin(), wstart.end(), random01);

    // worm update
    for (std::vector<boost::tuple<int, int, double> >::iterator wsi = wstart.begin();
         wsi != wstart.end(); ++wsi) {
      int direc = 2 * random01();
      int site, stp;
      double time;
      boost::tie(site, stp, time) = *wsi;
      int stp_start = stp;
      double time_start = time;
      wcount += 1;
      wlength += (direc == 0) ? time : -time;
      while (true) {
        stp = (direc == 0) ? stpoints[stp].prev() : stpoints[stp].next();
        if (stpoints[stp].at_operator()) {
          int bop = stpoints[stp].bond_operator();
          time = operators[bop].time();
          wlength += (direc == 0) ? -time : time;
          int ent = ((direc ^ 1) << 1) | stpoints[stp].leg();
          operators[bop].flip_state(ent);
          int ext = markov[operators[bop].state()](ent, random01);
          int leg = ext & 1;
          operators[bop].flip_state(ext);
          direc = (ext >> 1);
          site = (leg == 0) ? operators[bop].site0() : operators[bop].site1();
          stp = (leg == 0) ? operators[bop].stp0() : operators[bop].stp1();
          wlength += (direc == 0) ? time : -time;
        } else if (stpoints[stp].at_origin()) {
          spins[site] ^= 1;
          wlength += 1;
          time  = (direc == 0) ? 1 : 0;
        } else {
          if (stp == stp_start) {
            wlength += (direc == 0) ? -time_start : time_start;
            break;
          }
        }
      }
    }
    check_operators(lattice, spins, operators, stpoints);

    // measurement of physical quantities
    if (mcs >= opt.therm) {
      ene << (lattice.num_bonds() * wt.offset() - operators.size() / beta) / lattice.num_sites();
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
    
    if (mcs <= opt.therm / 2) {
      wdensity = opt.L / (wlength / wcount);
      if (mcs % (opt.therm / 8) == 0) {
        wcount /= 2;
        wlength /= 2;
      }
    }
    if (mcs == opt.therm / 2)
      std::cout << "Info: average number worms per MCS is reset from " << opt.L << " to "
                << wdensity << "\n\n";
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << (opt.therm + opt.sweeps) / elapsed << " MCS/sec\n";
  std::cout << "Energy Density            = "
            << ene.mean() << " +- " << ene.error() << std::endl
            << "Uniform Magnetization     = "
            << umag.mean() << " +- " << umag.error() << std::endl
            << "Uniform Magnetization^2   = "
            << umag2.mean() << " +- " << umag2.error() << std::endl
            << "Staggered Magnetization^2 = "
            << smag2.mean() << " +- " << smag2.error() << std::endl;
}
