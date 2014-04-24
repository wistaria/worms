##############################################################################
#
# worms: a simple worm code
#
# Copyright (C) 2013-2014 by Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
##############################################################################

set(TESTS outgoing_weight weight)
foreach(name ${TESTS})
  add_executable(${name} ${name}.cpp)
  add_test(${name} ${name})
endforeach(name)
unset(TESTS)
