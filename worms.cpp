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
#include <vector>

int main(int argc, char* argv[]) {
  std::cout << "worms: a simple worm code\n";
  options opt(argc, argv, 16, 1.0);
  if (!opt.valid) std::exit(-1);
}
