cmake_minimum_required(VERSION 2.8)
project(Transfinite)

cmake_policy(SET CMP0015 NEW)   # link directories are relative to source

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++1y")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(CMAKE_COMPILER_IS_GNUCXX)

set(USE_OWN_GEOMETRY_LIB OFF)

add_subdirectory(geom)
add_subdirectory(transfinite)
add_subdirectory(test)

if(NOT USE_OWN_GEOMETRY_LIB)
  add_subdirectory(geometry)
  add_subdirectory(pytransurf)
endif()
