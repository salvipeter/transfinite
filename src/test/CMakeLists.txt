include_directories(../geom)
set(GEOM_LIB geom)

include_directories(../transfinite)

add_executable(transfinite-test test.cc gb-fit.cc io.cc bezier.cc nelder-mead.cc)

add_dependencies(transfinite-test geom)

target_link_libraries(transfinite-test ${GEOM_LIB} transfinite)
