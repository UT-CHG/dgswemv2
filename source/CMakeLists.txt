set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall")

list (APPEND DG_HYPER_SWE_SOURCES
  ${PROJECT_SOURCE_DIR}/source/preprocessor/input_parameters.cpp

  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp

  ${PROJECT_SOURCE_DIR}/source/simulation/stepper.cpp

  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp

  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp

  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp

  ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
)

include_directories(${YAML_CPP_INCLUDE_DIR} ${MPI_C_INCLUDE_PATH})

add_executable(
  DG_HYPER_SWE
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/ompi_main_swe.cpp
  ${DG_HYPER_SWE_SOURCES}
)

target_link_libraries(DG_HYPER_SWE ${YAML_CPP_LIBRARIES} ${MPI_C_LIBRARIES})
