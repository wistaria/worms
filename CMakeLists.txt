cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(worms CXX)
include(cmake/postfix.cmake)

include_directories(${worms_SOURCE_DIR})

include(eigen3)
include(googletest)
include(standards)
include(bcl)
include(lattice)
FetchContent_MakeAvailable(googletest standards bcl)
include_directories(${FetchContent_includes})

add_subdirectory(test)

FOREACH(name worms)
  add_executable(${name} ${name}.cpp)
  add_test(${name} ${name})
ENDFOREACH(name)
