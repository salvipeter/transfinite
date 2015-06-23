include_directories(../geometry)

add_executable(test test.cc)
target_link_libraries(test geometry)
