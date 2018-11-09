list (APPEND PARTITIONER_SOURCES
  ${PROJECT_SOURCE_DIR}/partitioner/partition.cpp
  ${PROJECT_SOURCE_DIR}/partitioner/write_distributed_interfaces.cpp

  ${PROJECT_SOURCE_DIR}/partitioner/problem/swe_partitioner_inputs.cpp
  ${PROJECT_SOURCE_DIR}/partitioner/problem/default_partitioner_inputs.cpp
)

add_executable(
  partitioner
  partitioner.cpp
  ${SOURCES}
  ${PARTITIONER_SOURCES}
)

target_compile_definitions(partitioner PRIVATE ${LINALG_DEFINITION})
target_include_directories(partitioner PRIVATE ${YAML_CPP_INCLUDE_DIR} ${METIS_INCLUDE_DIRS})
target_link_libraries(partitioner ${YAML_CPP_LIBRARIES} ${METIS_LIBRARIES})
install(TARGETS partitioner DESTINATION bin/utilities/partitioner)
