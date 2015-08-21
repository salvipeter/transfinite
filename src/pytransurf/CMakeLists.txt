find_package(PythonLibs)

set(GEOM_LIB -L../geometry -lgeometry)
include_directories(../geometry ../transfinite ${PYTHON_INCLUDE_DIRS})

add_library(pytransurf SHARED pytransurf.cc)

set_target_properties(pytransurf PROPERTIES PREFIX "")
target_link_libraries(pytransurf ${GEOM_LIB} transfinite)
add_dependencies(pytransurf geometry)
