include_directories(../geom)
set(GEOM_LIB geom)

include_directories(../transfinite)
include_directories(../utils)

add_executable(transfinite-test test.cc)

add_dependencies(transfinite-test geom)

target_link_libraries(transfinite-test ${GEOM_LIB} transfinite transfinite-utils)
