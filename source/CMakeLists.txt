set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g")

list (APPEND DG_HYPER_SWE_SOURCES
  ${PROJECT_SOURCE_DIR}/source/rkdg_problem/SWE/problem_input/swe_inputs.cpp
  ${PROJECT_SOURCE_DIR}/source/rkdg_problem/SWE/problem_parser/swe_parser.cpp

  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp

  ${PROJECT_SOURCE_DIR}/source/simulation/rkdg_simulation/rkdg_stepper.cpp

  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp

  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp

  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
)

add_executable(
  DG_HYPER_SWE_SERIAL
  ${PROJECT_SOURCE_DIR}/source/rkdg_problem/SWE/main_files/main_swe.cpp
  ${DG_HYPER_SWE_SOURCES}
)

target_include_directories(DG_HYPER_SWE_SERIAL PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(DG_HYPER_SWE_SERIAL ${YAML_CPP_LIBRARIES})

install(TARGETS DG_HYPER_SWE_SERIAL DESTINATION bin)

if(USE_OMPI)
  add_executable(
    DG_HYPER_SWE_OMPI
    ${PROJECT_SOURCE_DIR}/source/rkdg_problem/SWE/main_files/ompi_main_swe.cpp
    ${DG_HYPER_SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  target_include_directories(DG_HYPER_SWE_OMPI PRIVATE ${YAML_CPP_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH})  
  target_link_libraries(DG_HYPER_SWE_OMPI ${YAML_CPP_LIBRARIES} ${MPI_CXX_LIBRARIES})

  install(TARGETS DG_HYPER_SWE_OMPI DESTINATION bin)
endif()

if(USE_HPX)
  add_executable(
    DG_HYPER_SWE_HPX
    ${PROJECT_SOURCE_DIR}/source/rkdg_problem/SWE/main_files/hpx_main_swe.cpp
    ${DG_HYPER_SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/hpx_communicator.cpp
  )

  hpx_setup_target(
    DG_HYPER_SWE_HPX
    COMPONENT_DEPENDENCIES iostreams
  )

  target_include_directories(DG_HYPER_SWE_HPX PRIVATE ${YAML_CPP_INCLUDE_DIR} ${HPX_INCLUDE_DIRS})  
  target_link_libraries(DG_HYPER_SWE_HPX ${YAML_CPP_LIBRARIES} ${HPX_LIBRARIES})

  install(TARGETS DG_HYPER_SWE_HPX DESTINATION bin)
endif()