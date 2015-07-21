find_package(PythonLibs)

if(USE_OWN_GEOMETRY_LIB)
  include_directories(../geom)
  set(GEOM_LIB geom)
else()
  include_directories(../geometry)
  set(GEOM_LIB -L../geometry -lgeometry)
endif()

include_directories(../transfinite ${PYTHON_INCLUDE_DIRS})

add_library(pytransurf SHARED pytransurf.cc)

set_target_properties(pytransurf PROPERTIES PREFIX "")
target_link_libraries(pytransurf ${GEOM_LIB} transfinite)
