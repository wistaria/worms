cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)
project(worms)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/bcl/config)

find_package(BoostSrc)

enable_testing()
include(CTest)

include_directories(${PROJECT_SOURCE_DIR}/bcl)
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${Boost_ROOT_DIR})

FOREACH(name worms)
  add_executable(${name} ${name}.cpp)
ENDFOREACH(name)
