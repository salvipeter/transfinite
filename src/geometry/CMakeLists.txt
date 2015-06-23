set(SKETCHSRC "../../../sketches/engine/src")
set(SKETCHLIB "../../../sketches/engine/debug")

add_definitions(-DNO_AFFIN_EXPR)
add_definitions(-DEXPLICIT_TEMPLATES)

include_directories(. ${SKETCHSRC}/Globals ${SKETCHSRC}/Affin ${SKETCHSRC}/Curves
                    ${SKETCHSRC}/Mesh ${SKETCHSRC}/Numeric/Misc ${SKETCHSRC}/BaseFunctions
                    ${SKETCHSRC}/Container)


link_directories(${SKETCHLIB}/BaseFunctions ${SKETCHLIB}/Curves ${SKETCHLIB}/Globals
                 ${SKETCHLIB}/Numeric/LinAlgTypes ${SKETCHLIB}/Numeric/Misc)

add_library(geometry SHARED geometry.cc)

target_link_libraries(geometry BaseFunctions Curves LinAlgTypes Misc Globals)
