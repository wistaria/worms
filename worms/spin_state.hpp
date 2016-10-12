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

template<int Spin> struct n_bits{ static const int val = 1;};
template<> struct n_bits<2> { static const int val = 1;};
template<> struct n_bits<3> { static const int val = 2;};
template<> struct n_bits<4> { static const int val = 2;};
template<> struct n_bits<5> { static const int val = 3;};
      
template<int Spin> struct sitemask
{ static const int val = (1 << n_bits<Spin>::val) - 1;};
template<int Spin> struct pairmask
{ static const int val = (1 << (n_bits<Spin>::val*2)) - 1;};

template<int Spin, int Site> struct siteshift
{ static const int val = n_bits<Spin>::val*Site;};
template<int Spin, int Site> struct pairshift
{ static const int val = 2*n_bits<Spin>::val*Site;};

template<int Spin> struct operatorsize
{ static const int val = 1 << (4*n_bits<Spin>::val); };

namespace spin_state {
  // inline int p2c(int p, int l) { return (p >> l) & 1; }        // return state at leg l \in {0,1,2,3}
  // inline int p2u(int p, int d) { return (p >> (2 * d)) & 3; }  // return state at site pair u \in {0,1}
  // inline int c2u(int c0, int c1) { return (c0 | (c1 << 1)); }
  // inline int c2p(int c0, int c1, int c2, int c3) {
  //   return (c0 | (c1 << 1) | (c2 << 2) | (c3 << 3));
  // }
  // inline int u2p(int u0, int u1) { return (u0 | (u1 << 2)); }
  // inline int maskp(int l) { return (1 << l); }

  template<int Spin> inline int p2c(const int& p, const int& l)
  { return (p >> (l * n_bits<Spin>::val)) & sitemask<Spin>::val; }
   
  template<int Spin> inline int p2u(const int& p, const int& d)
  { return (p >> (2 * d * n_bits<Spin>::val)) & pairmask<Spin>::val; }  // return state at site pair u \in {0,1}

  template<int Spin> inline int c2u(const int& c0, const int& c1)
  { return (c0 | (c1 << siteshift<Spin, 1>::val)); }
  template<int Spin> inline int c2p(const int& c0, const int& c1,
				    const int& c2, const int& c3) {
    return (c0 | (c1 << siteshift<Spin, 1>::val) |
	    (c2 << siteshift<Spin, 2>::val) | (c3 << siteshift<Spin, 3>::val));
  }
  template<int Spin> inline int u2p(const int& u0, const int& u1)
  { return (u0 | (u1 << pairshift<Spin, 1>::val)); }
  template<int Spin> inline int u2c(const int& u, const int& l)
  { return (u >> (l * n_bits<Spin>::val)) & sitemask<Spin>::val; }



  template<int Spin> inline int maskp(const int& l)
  { return (sitemask<Spin>::val << (l*n_bits<Spin>::val)); }

  template<int Spin> inline bool raise(const int& site, int *state){
    const int val = p2c<Spin>(*state, site);
    if (val < Spin-1){
      *state &= ~maskp<Spin>(site);
      *state |= (val + 1) << (site * n_bits<Spin>::val);
      return true;
    }
    return false;
  }

  template<int Spin> inline bool lower(const int& site, int *state){
    const int val = p2c<Spin>(*state, site);
    if (val > 0){
      *state &= ~maskp<Spin>(site);
      *state |= (val - 1) << (site * n_bits<Spin>::val);
      return true;
    }
    return false;
  }

  template<int Spin> inline bool valid_c(const int& c)
  { return ((0 <= c) && (c < Spin)); }

  template<int Spin> inline bool valid_u(const int& u)
  { return (valid_c<Spin>(u2c<Spin>(u,0)) && valid_c<Spin>(u2c<Spin>(u,1))); }

  template<int Spin> inline bool valid_p(const int& p)
  { return (valid_u<Spin>(p2u<Spin>(p,0)) && valid_u<Spin>(p2u<Spin>(p,1))); }
    
};

#endif // SPIN_STATE_HPP
