set(SKETCHSRC "../../../sketches/engine/src")

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(SKETCHLIB "../../../sketches/engine/debug")
else()
  set(SKETCHLIB "../../../sketches/engine/release")
endif()

add_definitions(-DNO_AFFIN_EXPR)
add_definitions(-DEXPLICIT_TEMPLATES)

include_directories(.
  ${SKETCHSRC}/Affin
  ${SKETCHSRC}/Approx/ApproxCrv
  ${SKETCHSRC}/Approx/ApproxCrv/BSCF_Base
  ${SKETCHSRC}/Approx/ApproxCrv/BSCF_FN
  ${SKETCHSRC}/Approx/ApproxCrv/BSCF_KI
  ${SKETCHSRC}/Approx/ApproxSurf
  ${SKETCHSRC}/Approx/ApproxSurf/BSSF_Base
  ${SKETCHSRC}/Approx/ApproxSurf/BSSF_FN
  ${SKETCHSRC}/Approx/ApproxSurf/BSSF_KI
  ${SKETCHSRC}/Approx/ApproxUtil
  ${SKETCHSRC}/BaseFunctions
  ${SKETCHSRC}/Container
  ${SKETCHSRC}/Curves
  ${SKETCHSRC}/Filter/IGES
  ${SKETCHSRC}/Flattener4.1/src
  ${SKETCHSRC}/Globals
  ${SKETCHSRC}/Mesh
  ${SKETCHSRC}/Numeric/Factorize
  ${SKETCHSRC}/Numeric/LinAlgTypes
  ${SKETCHSRC}/Numeric/Misc
  ${SKETCHSRC}/Surfaces
  )


link_directories(
  ${SKETCHLIB}/Approx/ApproxCrv
  ${SKETCHLIB}/Approx/ApproxCrv/BSCF_Base
  ${SKETCHLIB}/Approx/ApproxCrv/BSCF_FN
  ${SKETCHLIB}/Approx/ApproxCrv/BSCF_KI
  ${SKETCHLIB}/Approx/ApproxSurf
  ${SKETCHLIB}/Approx/ApproxSurf/BSSF_Base
  ${SKETCHLIB}/Approx/ApproxSurf/BSSF_FN
  ${SKETCHLIB}/Approx/ApproxSurf/BSSF_KI
  ${SKETCHLIB}/Approx/ApproxUtil
  ${SKETCHLIB}/BaseFunctions
  ${SKETCHLIB}/Curves
  ${SKETCHLIB}/Filter/IGES
  ${SKETCHLIB}/Flattener4.1/src
  ${SKETCHLIB}/Globals
  ${SKETCHLIB}/Mesh
  ${SKETCHLIB}/Numeric/Factorize
  ${SKETCHLIB}/Numeric/LinAlgTypes
  ${SKETCHLIB}/Numeric/Misc
  ${SKETCHLIB}/Surfaces
  )

add_library(geometry geometry.cc)

target_link_libraries(geometry
  IGES
  ApproxSurf ApproxCrv ApproxUtil BSSF_KI BSSF_FN BSSF_Base BSCF_KI BSCF_FN BSCF_Base
  Surfaces Curves BaseFunctions
  tmf4
  Mesh
  Misc Factorize LinAlgTypes
  Globals
  )
