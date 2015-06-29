if(USE_OWN_GEOMETRY_LIB)
  include_directories(../geom)
else()
  include_directories(../geometry)
endif()

add_library(transfinite
  rmf.cc
  domain.cc domain-regular.cc domain-circular.cc
  parameterization.cc parameterization-bilinear.cc parameterization-barycentric.cc
  ribbon.cc ribbon-compatible.cc
)
