include_directories(../geometry)

add_executable(test test.cc)

set(GEOM_LIB -L../geometry -lgeometry)

target_link_libraries(test ${GEOM_LIB})
