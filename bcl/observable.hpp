/*****************************************************************************
*
* BCL: Balance Condition Library
*
* Copyright (C) 2009-2013 by Hidemaro Suwa <suwamaro@looper.t.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef BCL_OBSERVABLE_HPP
#define BCL_OBSERVABLE_HPP

#include <cmath> // for std::sqrt

namespace bcl {

class observable {
public:
  observable() : count_(0), sum_(0), esq_(0) {}
  void operator<<(double x) { sum_ += x; esq_ += x * x; ++count_; }
  double mean() const { return (count_ > 0) ? (sum_ / count_) : 0.; }
  double error() const {
    return (count_ > 1) ?
      std::sqrt((esq_ / count_ - mean() * mean()) / (count_ - 1)) : 0.;
  }
private:
  unsigned int count_;
  double sum_, esq_;
};

} // end namespace bcl

#endif // BCL_OBSERVABLE_HPP
