set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++14 -Wall")

add_executable(
  rectangular_mesh_generator
  ${PROJECT_SOURCE_DIR}/mesh_generators/rectangular_mesh_generator.cpp
)