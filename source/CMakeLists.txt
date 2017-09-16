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
)

add_executable(
  DG_HYPER_SWE_SERIAL
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/main_swe.cpp
  ${DG_HYPER_SWE_SOURCES}
)

target_link_libraries(DG_HYPER_SWE_SERIAL ${YAML_CPP_LIBRARIES})

if(USE_OMPI)
  add_executable(
    DG_HYPER_SWE_OMPI
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/ompi_main_swe.cpp
    ${DG_HYPER_SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  target_link_libraries(DG_HYPER_SWE_OMPI ${YAML_CPP_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif()

if(USE_HPX)
  add_executable(
    DG_HYPER_SWE_HPX
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/hpx_main_swe.cpp
    ${DG_HYPER_SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/hpx_communicator.cpp
  )


  hpx_setup_target(
    DG_HYPER_SWE_HPX
    COMPONENT_DEPENDENCIES iostreams
  )

  target_link_libraries(DG_HYPER_SWE_HPX ${YAML_CPP_LIBRARIES} ${HPX_LIBRARIES})
endif()