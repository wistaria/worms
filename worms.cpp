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



// #define DEBUG 1

// put 2S+1 here
#define SPIN 5

#ifdef DEBUG
# define CHECK_OPERATORS 1
#endif

#include "options.hpp"
#include "worms/version.hpp"
#include "worms/chain_lattice.hpp"
#include "worms/heisenberg_operator.hpp"
#include "worms/operations.hpp"
#include "worms/weight.hpp"
#include "worms/outgoing_weight.hpp"

#include "bcl.hpp"

#include <vector>
#include <fstream>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <boost/tuple/tuple.hpp>

// function to decide whether worm raises operator ahead
inline int raise_ahead(int raise_lower, int direc){
  return (~(raise_lower^direc))&1;
}

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
  std::vector<int> spins(opt.L, 0 /* all down */);
  std::vector<bond_operator<SPIN> > operators, operators_p;
  std::vector<spacetime_point> stpoints;

  // Hamiltonian operator
  heisenberg_operator<SPIN> op(opt.H, /* coordination number = */ 2);
  double offset = 0.3; // >= 0
  weight<SPIN> wt(op, offset);
  
  // table for diagonal update
  double lambda = wt.max_diagonal_weight();
  std::vector<double> accept(1 << (2*n_bits<SPIN>::val));
  for (int c0 = 0; c0 < SPIN; ++c0)
    for (int c1 = 0; c1 < SPIN; ++c1)
      accept[spin_state::c2u<SPIN>(c0, c1)] =
	wt[spin_state::c2p<SPIN>(c0, c1, c0, c1)] / lambda;

  
  // prepare tables for transition probabilities in worm update
  typedef bcl::markov<random01_t> markov_t;
  std::vector<markov_t> markov_raise_top;
  std::vector<markov_t> markov_raise_bottom;

  std::vector<std::vector<double> > raise_top_weights;
  std::vector<std::vector<double> > raise_bottom_weights;

  markov_raise_top.resize(operatorsize<SPIN>::val);
  markov_raise_bottom.resize(operatorsize<SPIN>::val);
  raise_top_weights.resize(operatorsize<SPIN>::val);
  raise_bottom_weights.resize(operatorsize<SPIN>::val);
  for (int c = 0; c < operatorsize<SPIN>::val; ++c){
    raise_top_weights[c].resize(4);
    raise_bottom_weights[c].resize(4);

    // Compute weights for raising worm
    for (int ext_leg = 0; ext_leg < 4; ++ext_leg){
      int direc = (ext_leg >> 1);
      int state = c;

      // Compute weights for raise top worms
      if (direc==1){
	if (spin_state::raise<SPIN>(ext_leg, &state))
	  raise_top_weights[c][ext_leg] = wt[state];
	else
	  raise_top_weights[c][ext_leg] = 0;
      }
      else {
	if (spin_state::lower<SPIN>(ext_leg, &state))
	  raise_top_weights[c][ext_leg] = wt[state];
	else
	  raise_top_weights[c][ext_leg] = 0;
      }
      markov_raise_top[c] = markov_t(bcl::st2010(), raise_top_weights[c]);  
      

      // Compute weights for raise bottom worm
      state = c;
      if (direc==1){
	if (spin_state::lower<SPIN>(ext_leg, &state))
	  raise_bottom_weights[c][ext_leg] = wt[state];
	else
	  raise_bottom_weights[c][ext_leg] = 0;
      }
      else {
	if (spin_state::raise<SPIN>(ext_leg, &state))
	  raise_bottom_weights[c][ext_leg] = wt[state];
	else
	  raise_bottom_weights[c][ext_leg] = 0;
      }
      markov_raise_bottom[c] = markov_t(bcl::st2010(), raise_bottom_weights[c]);  
    }

  }


  // weight for worm insertion
  double wdensity = opt.L;
  std::vector<boost::tuple<int, int, double> > wstart;

  // worm statistics
  double wcount = 0;
  double wlength = 0;
  
  // temporaries
  std::vector<double> times;         
  std::vector<int> current(opt.L);   // stores current spin configuration 

  boost::timer tm;

  std::ofstream debugfile;

  for (int mcs = 0; mcs < (opt.therm + opt.sweeps); ++mcs) { 
    if (mcs % 10000 == 0)
      std::cout << "mcs " << mcs << std::endl;
    
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
    operators_p.push_back(bond_operator<SPIN>(0, 0, 0, 0, 0, 1));  // sentinel
    std::copy(spins.begin(), spins.end(), current.begin());  
    stpoints.resize(0);                                      
    for (int s = 0; s < opt.L; ++s)
      stpoints.push_back(spacetime_point::origin(s)); 
    wstart.resize(0);   
    check_operators(lattice, spins, operators, stpoints);
    std::vector<double>::iterator tmi = times.begin();
    for (std::vector<bond_operator<SPIN> >::iterator opi = operators_p.begin();
	 //insert either worm or diagonal operator at every time
         opi != operators_p.end();) {
      if (*tmi < opi->time()) {
        if (random01() < pstart) {
          // insert worm starting point
          int s = static_cast<int>(opt.L * random01());
          append_wstart(s, *tmi, stpoints, wstart);
        } else {
          // insert diagonal operator
          int b = static_cast<int>(opt.L * random01());  // choose random position of operator 
          int s0 = lattice.source(b);
          int s1 = lattice.target(b);
          int u = spin_state::c2u<SPIN>(current[s0], current[s1]);

	  double random = random01();

          if (random < accept[u])
            append_operator(s0, s1, spin_state::u2p<SPIN>(u, u), *tmi, operators, stpoints);
        }
        ++tmi;
      } else {
        if (opi->is_offdiagonal()) {
          // keep offdiagonal operator
          append_operator(*opi, operators, stpoints);
          current[opi->site0()] = spin_state::p2c<SPIN>(opi->state(), 2);
          current[opi->site1()] = spin_state::p2c<SPIN>(opi->state(), 3); 
       }
        ++opi;
      }
    }
    check_operators(lattice, spins, operators, stpoints);

    // shuffle starting point of worms
    bcl::random_shuffle(wstart.begin(), wstart.end(), random01);
    
    // worm update
    int worm_ctr = 0;
    for (std::vector<boost::tuple<int, int, double> >::iterator wsi = wstart.begin();
         wsi != wstart.end(); ++wsi, ++worm_ctr) {
      int direc = 2 * random01();  // run either forward or backward
      int site, stp;               // these three variables define position of wormhead
      double time;
      boost::tie(site, stp, time) = *wsi;   // set site, stp, and time to worm starting point
      int stp_start = stp;
      double time_start = time;
      wcount += 1;
      wlength += (direc == 0) ? time : -time;
      
      check_operators(lattice, spins, operators, stpoints);

      // Find current Spin in worldline
      int stp_next = stp_start;
      int next_op = -3;
      while (next_op < -1){
	stp_next = (direc == 0) ? stpoints[stp_next].prev() : stpoints[stp_next].next();
	next_op = stpoints[stp_next].bond_operator();
      }
      int next_ent = ((direc ^ 1) << 1) | stpoints[stp_next].leg();
      int next_spinstate;

      if (next_op == -1) // at t=0
	next_spinstate = spins[site];
      else
      	next_spinstate = spin_state::p2c<SPIN>(operators[next_op].state(), next_ent);

      // Decide on type of worm (raising/lowering)
      int raise_lower;
      raise_lower = 2 * random01(); // random worm

      // Reject worms that don't work
      if (((next_spinstate == 0) && (raise_ahead(raise_lower, direc) == 0))||
	  ((next_spinstate == SPIN-1) && (raise_ahead(raise_lower, direc) == 1)))
	continue;
      
      check_operators(lattice, spins, operators, stpoints);
      
      while (true) {
	
	stp = (direc == 0) ? stpoints[stp].prev() : stpoints[stp].next();
        if (stpoints[stp].at_operator()) {
	  int bop = stpoints[stp].bond_operator();   
          time = operators[bop].time();              
          site = (stpoints[stp].leg() == 0) ?
	    operators[bop].site0() : operators[bop].site1();

	  wlength += (direc == 0) ? -time : time;
          int ent = ((direc ^ 1) << 1) | stpoints[stp].leg();  	  
	  
	  // Choose exit
	  operators[bop].flip_state(ent, raise_ahead(raise_lower, direc));
	  int ext;
	  
	  if (raise_lower==1)
	    ext = markov_raise_top[operators[bop].state()](ent, random01); 
	  else
	    ext = markov_raise_bottom[operators[bop].state()](ent, random01);

	  int leg = ext & 1;
	  direc = (ext >> 1);
	  // assert(ext!=ent); // Assertion holds if code is bouncefree

	  // flip state of operator at exit point
          operators[bop].flip_state(ext, raise_ahead(raise_lower, direc));         

          site = (leg == 0) ? operators[bop].site0() : operators[bop].site1(); 
          stp = (leg == 0) ? operators[bop].stp0() : operators[bop].stp1();
          wlength += (direc == 0) ? time : -time;

        } else if (stpoints[stp].at_origin()) {

	  if (raise_ahead(raise_lower, direc))
	    spins[site] += 1;
	  else
	    spins[site] -= 1;
          wlength += 1;
          time  = (direc == 0) ? 1 : 0;

        } else {
	  
	  if (stp == stp_start) {  
            wlength += (direc == 0) ? -time_start : time_start;
            break; // break the infinite loop 
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
