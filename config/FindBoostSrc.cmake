# - FindBoostSrc.cmake

set(Boost_DEBUG 1)

if (NOT Boost_ROOT_DIR)
  if(BOOST_ROOT)
    find_path(Boost_ROOT_DIR "libs/program_options/src/cmdline.cpp" ${BOOST_ROOT})
  else(BOOST_ROOT)
    set(DIR0
      ${PROJECT_SOURCE_DIR}/..
      $ENV{HOME} $ENV{HOME}/src $ENV{HOME}/ALPS/src /usr/local /usr/local/src /opt /opt/local /opt/alps /opt/local/include
      "$ENV{HOMEDRIVE}/Program Files"
      "$ENV{HOMEDRIVE}$ENV{HOMEPATH}"
      "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/src"
      "$ENV{HOMEDRIVE}$ENV{HOMEPATH}/ALPS/src")
    # FIXIT: | TODO: | TBD: THIS IS UGLY!
    set(DIR1 boost boostsrc boost_1_52_0 boost_1_51_0 boost_1_50_0 boost_1_49_0 boost_1_48_0 boost_1_47_0)
    set(_boost_SEARCH_PATH "")
    foreach(D0 ${DIR0})
      foreach(D1 ${DIR1})
        set(_boost_SEARCH_PATH ${_boost_SEARCH_PATH} ${D0}/${D1})
      endforeach(D1)
    endforeach(D0)
    find_path(Boost_ROOT_DIR "libs/program_options/src/cmdline.cpp" ${_boost_SEARCH_PATH})
  endif(BOOST_ROOT)
endif (NOT Boost_ROOT_DIR)

if(Boost_ROOT_DIR)
  set(Boost_INCLUDE_DIR ${Boost_ROOT_DIR} CACHE PATH "Boost Include Directory")
  if (NOT Boost_INCLUDE_DIR)
    set(Boost_INCLUDE_DIR ${Boost_ROOT_DIR} CACHE PATH "Boost Include Directory" FORCE)
  endif (NOT Boost_INCLUDE_DIR)
endif(Boost_ROOT_DIR)

# check Boost version (from Modules/FindBoost.cmake)
if(Boost_INCLUDE_DIR)
  set(_boost_VERSION 0)
  set(_boost_LIB_VERSION "")
  file(READ "${Boost_INCLUDE_DIR}/boost/version.hpp" _boost_VERSION_HPP_CONTENTS)
  string(REGEX REPLACE ".*#define BOOST_VERSION ([0-9]+).*" "\\1" _boost_VERSION "${_boost_VERSION_HPP_CONTENTS}")
  string(REGEX REPLACE ".*#define BOOST_LIB_VERSION \"([0-9_]+)\".*" "\\1" _boost_LIB_VERSION "${_boost_VERSION_HPP_CONTENTS}")
  set(Boost_LIB_VERSION ${_boost_LIB_VERSION} CACHE INTERNAL "The library version string for boost libraries")
  set(Boost_VERSION ${_boost_VERSION} CACHE INTERNAL "The version number for boost libraries")
  MATH(EXPR Boost_MAJOR_VERSION "${Boost_VERSION} / 100000")
  MATH(EXPR Boost_MINOR_VERSION "${Boost_VERSION} / 100 % 1000")
  MATH(EXPR Boost_SUBMINOR_VERSION "${Boost_VERSION} % 100")
endif(Boost_INCLUDE_DIR)

if(Boost_ROOT_DIR)
  message(STATUS "Found Boost Source: ${Boost_ROOT_DIR}")
  message(STATUS "Boost Version: ${Boost_MAJOR_VERSION}_${Boost_MINOR_VERSION}_${Boost_SUBMINOR_VERSION}")
  if(Boost_MAJOR_VERSION LESS 1 OR Boost_MINOR_VERSION LESS 47)
    message(FATAL_ERROR "Boost library version is too old, boost chrono requires boost >= 1.47.0")
  endif(Boost_MAJOR_VERSION LESS 1 OR Boost_MINOR_VERSION LESS 47)
else(Boost_ROOT_DIR)
  message(FATAL_ERROR "Boost Source not Found")
endif(Boost_ROOT_DIR)

# Avoid auto link of Boost library
add_definitions(-DBOOST_ALL_NO_LIB=1)
if(BUILD_SHARED_LIBS)
  add_definitions(-DBOOST_ALL_DYN_LINK=1)
else(BUILD_SHARED_LIBS)
  add_definitions(-DBOOST_ALL_STATIC_LINK=1)
endif(BUILD_SHARED_LIBS)

mark_as_advanced(Boost_INCLUDE_DIR)
