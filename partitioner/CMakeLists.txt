set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall")

find_package(METIS REQUIRED)

list (APPEND PARTITIONER_SOURCES
  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp

  ${PROJECT_SOURCE_DIR}/partitioner/partition.cpp
)

include_directories(${METIS_INCLUDE_DIRS})

add_executable(
  partitioner
  partitioner.cpp
  ${PARTITIONER_SOURCES}
)

target_link_libraries(partitioner ${METIS_LIBRARIES})
