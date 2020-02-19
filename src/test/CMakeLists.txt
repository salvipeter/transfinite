include_directories(../geom)
set(GEOM_LIB geom)

include_directories(../transfinite)

cmake_policy(SET CMP0037 OLD)   # allow executables with reserved names (like "test")

add_executable(test test.cc gb-fit.cc io.cc bezier.cc nelder-mead.cc)

add_dependencies(test geom)

target_link_libraries(test ${GEOM_LIB} transfinite)
