set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14 -Wall")

list (APPEND DG_HYPER_SWE_SOURCES
  ${PROJECT_SOURCE_DIR}/source/mesh_metadata.cpp
  ${PROJECT_SOURCE_DIR}/source/ADCIRC_reader/adcirc_format.cpp

  ${PROJECT_SOURCE_DIR}/source/stepper.cpp

  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp

  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp

  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
)

add_executable(
  DG_HYPER_SWE
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/main_swe.cpp
  ${DG_HYPER_SWE_SOURCES}
)