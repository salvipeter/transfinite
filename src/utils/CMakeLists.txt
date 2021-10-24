include_directories(../geom)

include_directories(../transfinite)

add_library(transfinite-utils STATIC gb-fit.cc io.cc bezier.cc nelder-mead.cc)

target_link_libraries(transfinite-utils geom transfinite)
