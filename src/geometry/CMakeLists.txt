set(SKETCHES "../../../sketches/engine/src")

add_definitions(-DNO_AFFIN_EXPR)
add_definitions(-DEXPLICIT_TEMPLATES)

include_directories(. ${SKETCHES}/Globals ${SKETCHES}/Affin ${SKETCHES}/Curves
                    ${SKETCHES}/Mesh ${SKETCHES}/Numeric/Misc ${SKETCHES}/BaseFunctions
                    ${SKETCHES}/Container)

add_library(geometry SHARED geometry.cc)
