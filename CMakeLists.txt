##############################################################################
#
# worms: a simple worm code
#
# Copyright (C) 2013-2018 by Synge Todo <wistaria@comp-phys.org>
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
##############################################################################

cmake_minimum_required(VERSION 3.1)
project(worms)

if(NOT CMAKE_BUILD_TYPE)
    message("Using default build type: Release")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

find_package(Boost)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR})

enable_testing()
include(CTest)
include(add_alps_test)
add_subdirectory(test)

FOREACH(name worms)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} ${Boost_LIBRARIES})
  add_alps_test(${name})
ENDFOREACH(name)
