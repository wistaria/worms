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

#ifndef HEISENBERG_OPERATOR_HPP
#define HEISENBERG_OPERATOR_HPP

#include <algorithm>
#include <vector>
#include <limits>

#include "spin_state.hpp"


template <int Spin>
class heisenberg_operator
{
public:
  heisenberg_operator(double h, int coord_num)
    : elements_(operatorsize<Spin>::val, 0),
      S((Spin-1.0)/2.)
  {
    // Diagonal Part
    for (int t1 = 0; t1 < Spin; ++t1)
      for (int t2 = 0; t2 < Spin; ++t2)
	{
	  const double m1 = double(t1) - S;
	  const double m2 = double(t2) - S;
	  elements_[spin_state::c2p<Spin>(t1, t2, t1, t2)] = m1*m2;
	}

    // Offdiagonal Part
    for (int t1 = 0; t1 < Spin; ++t1)
      for (int t2 = 0; t2 < Spin; ++t2)
	{
	  const double m1 = double(t1) - S;
	  const double m2 = double(t2) - S;
	  if ((t1 < (Spin - 1)) && (t2 > 0))
	    {
	      const double coeff = -std::sqrt( (S*(S+1) - m1*(m1+1)) * (S*(S+1) - m2*(m2-1)) ) / 2.;
	      // std::cout << t1 << " " << t2 << " " << t1+1 << " " << t2-1  << " " << coeff << std::endl;
	      elements_[spin_state::c2p<Spin>(t1, t2, t1+1, t2-1)] = coeff;
	    }
	  if ((t2 < (Spin - 1)) && (t1 > 0))
	    {
	      const double coeff = -std::sqrt( (S*(S+1) - m1*(m1-1)) * (S*(S+1) - m2*(m2+1)) ) / 2.;
	      // std::cout << t1 << " " << t2 << " " << t1-1 << " " << t2+1  << " " << coeff << std::endl;
	      elements_[spin_state::c2p<Spin>(t1, t2, t1-1, t2+1)] = coeff; 
	    }
	}
  }
  
  double operator()(int s) const { return elements_[s]; }
  double operator()(int sp0, int sp1, int sn0, int sn1) const
  { return elements_[spin_state::c2p<Spin>(sp0, sp1, sn0, sn1)]; }

  double max_diagonal() const
  {
    double max = std::numeric_limits<double>::min();
    for (int t1 = 0; t1 < Spin; ++t1)
      for (int t2 = 0; t2 < Spin; ++t2)
	if (elements_[spin_state::c2p<Spin>(t1, t2, t1, t2)] > max)
	  max = elements_[spin_state::c2p<Spin>(t1, t2, t1, t2)];
    return max;
  }
  
  double min_diagonal() const
  {
    double min = std::numeric_limits<double>::max();
    for (int t1 = 0; t1 < Spin; ++t1)
      for (int t2 = 0; t2 < Spin; ++t2)
	if (elements_[spin_state::c2p<Spin>(t1, t2, t1, t2)] < min)
	  min = elements_[spin_state::c2p<Spin>(t1, t2, t1, t2)];
    return min;
  }
    
private:
  const double S;
  std::vector<double> elements_;

};

// template <>
// class heisenberg_operator<2> {
// public:
//   heisenberg_operator(double h, int coord_num) : elements_(16, 0) {
//     elements_[spin_state::c2p<2>(0, 0, 0, 0)] = 0.25 - h / coord_num;
//     elements_[spin_state::c2p<2>(0, 1, 0, 1)] = -0.25;
//     elements_[spin_state::c2p<2>(1, 0, 1, 0)] = -0.25;
//     elements_[spin_state::c2p<2>(1, 1, 1, 1)] = 0.25 + h / coord_num;
//     elements_[spin_state::c2p<2>(0, 1, 1, 0)] = -0.5;
//     elements_[spin_state::c2p<2>(1, 0, 0, 1)] = -0.5;
//   }
//   double operator()(int s) const { return elements_[s]; }
//   double operator()(int sp0, int sp1, int sn0, int sn1) const {
//     return elements_[spin_state::c2p<2>(sp0, sp1, sn0, sn1)];
//   }
//   double max_diagonal() const {
//     return std::max(std::max(elements_[spin_state::c2p<2>(0, 0, 0, 0)],
//                              elements_[spin_state::c2p<2>(0, 1, 0, 1)]),
//                     std::max(elements_[spin_state::c2p<2>(1, 0, 1, 0)],
//                              elements_[spin_state::c2p<2>(1, 1, 1, 1)]));
//   }
//   double min_diagonal() const {
//     return std::min(std::min(elements_[spin_state::c2p<2>(0, 0, 0, 0)],
//                              elements_[spin_state::c2p<2>(0, 1, 0, 1)]),
//                     std::min(elements_[spin_state::c2p<2>(1, 0, 1, 0)],
//                              elements_[spin_state::c2p<2>(1, 1, 1, 1)]));
//   }
// private:
//   std::vector<double> elements_;
// };

#endif // HEISENBERG_OPERATOR_HPP
