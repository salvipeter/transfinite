set(SKETCHSRC "../../../sketches/engine/src")
set(SKETCHLIB "../../../sketches/engine/debug")

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
  ${SKETCHLIB}/Globals
  ${SKETCHLIB}/Mesh
  ${SKETCHLIB}/Numeric/Factorize
  ${SKETCHLIB}/Numeric/LinAlgTypes
  ${SKETCHLIB}/Numeric/Misc
  ${SKETCHLIB}/Surfaces
  )

add_library(geometry SHARED geometry.cc)

target_link_libraries(geometry
  BaseFunctions Curves Surfaces
  Mesh
  Factorize LinAlgTypes Misc
  Globals
  ApproxUtil
  ApproxCrv BSCF_Base BSCF_FN BSCF_KI
  ApproxSurf BSSF_Base BSSF_FN BSSF_KI
  )
