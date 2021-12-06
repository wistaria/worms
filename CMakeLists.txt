cmake_minimum_required(VERSION 3.14)
include(cmake/prefix.cmake)
project(worms CXX)
include(cmake/postfix.cmake)

include_directories(${worms_SOURCE_DIR})

include(googletest)
include(standards)
include(bcl)
FetchContent_MakeAvailable(googletest standards bcl)
include_directories(${standards_SOURCE_DIR} ${bcl_SOURCE_DIR})

add_subdirectory(test)

FOREACH(name worms)
  add_executable(${name} ${name}.cpp)
  add_test(${name} ${name})
ENDFOREACH(name)
