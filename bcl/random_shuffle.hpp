/*****************************************************************************
*
* BCL: Balance Condition Library
*
* Copyright (C) 2006-2013 by Synge Todo <wistaria@comp-phys.org>
*
* Distributed under the Boost Software License, Version 1.0. (See accompanying
* file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
*
*****************************************************************************/

#ifndef BCL_RANDOM_SHUFFLE_HPP
#define BCL_RANDOM_SHUFFLE_HPP

#include <algorithm>
#include <iterator>

namespace bcl {

template<class RandomAccessIter, class RandomNumberGenerator>
void random_shuffle(RandomAccessIter first, RandomAccessIter last,
                    RandomNumberGenerator& rng) {
  using std::iter_swap;
  for (typename std::iterator_traits<RandomAccessIter>::difference_type
         n = last - first; n > 1; ++first, --n)
    iter_swap(first, first + (int)(n * rng()));
}

} // end namespace bcl

#endif // BCL_RANDOM_SHUFFLE_HPP
