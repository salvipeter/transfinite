cmake_minimum_required(VERSION 2.8)
project(Transfinite)

cmake_policy(SET CMP0015 NEW)   # link directories are relative to source

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "-Wall -pedantic -std=c++1y -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(CMAKE_COMPILER_IS_GNUCXX)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "-Weverything -pedantic -std=c++11 -fPIC -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-missing-prototypes")
  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

set(USE_OWN_GEOMETRY_LIB ON)

add_subdirectory(transfinite)
add_subdirectory(test)

if(USE_OWN_GEOMETRY_LIB)
  add_subdirectory(geom)
endif()

if(NOT USE_OWN_GEOMETRY_LIB)
  add_subdirectory(geometry)
  add_subdirectory(pytransurf)
endif()
