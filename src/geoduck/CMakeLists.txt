set(GEOM_LIB -L../geometry -lgeometry)
include_directories(../geometry ../transfinite)

add_library(transfinite-kapi SHARED transfinite-kapi.cc)

target_link_libraries(transfinite-kapi ${GEOM_LIB} transfinite)
add_dependencies(transfinite-kapi geometry)
