if(USE_OWN_GEOMETRY_LIB)
  include_directories(../geom)
else()
  include_directories(../geometry)
endif()

add_library(transfinite rmf.cc)
