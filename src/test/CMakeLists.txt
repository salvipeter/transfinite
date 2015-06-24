if(USE_OWN_GEOMETRY_LIB)
  include_directories(../geom)
  set(GEOM_LIB geom)
else()
  include_directories(../geometry)
  set(GEOM_LIB -L../geometry -lgeometry)
endif()

include_directories(../transfinite)

add_executable(test test.cc)

target_link_libraries(test ${GEOM_LIB} transfinite)
