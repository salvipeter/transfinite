cmake_minimum_required(VERSION 3.2)
project(Transfinite)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(LibTriangle QUIET)

cmake_policy(SET CMP0015 NEW)   # link directories are relative to source

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|(Apple)?Clang")
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++17 -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

add_subdirectory(transfinite)
add_subdirectory(utils)
add_subdirectory(test)
add_subdirectory(geom)
