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
  ${SWE_SOURCES}
  ${PARTITIONER_SOURCES}
)

set_target_properties(partitioner PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
install(TARGETS partitioner DESTINATION bin/utilities/partitioner)
