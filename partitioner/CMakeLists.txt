set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall")

list (APPEND PARTITIONER_SOURCES
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_input/swe_inputs.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp

  ${PROJECT_SOURCE_DIR}/partitioner/partition.cpp
  ${PROJECT_SOURCE_DIR}/partitioner/write_distributed_interfaces.cpp

  ${PROJECT_SOURCE_DIR}/partitioner/problem/swe_partitioner_inputs.cpp
  ${PROJECT_SOURCE_DIR}/partitioner/problem/default_partitioner_inputs.cpp
)

add_executable(
  partitioner
  partitioner.cpp
  ${PARTITIONER_SOURCES}
)

target_include_directories(partitioner PRIVATE ${YAML_CPP_INCLUDE_DIR} ${METIS_INCLUDE_DIRS})
target_link_libraries(partitioner ${YAML_CPP_LIBRARIES} ${METIS_LIBRARIES})

install(TARGETS partitioner DESTINATION bin/utilities/partitioner)
