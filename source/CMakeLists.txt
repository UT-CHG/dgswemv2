add_executable(
  dgswemv2-serial
  ${PROJECT_SOURCE_DIR}/source/dgswemv2-serial.cpp
  ${SOURCES}
  ${PROJECT_SOURCE_DIR}/source/simulation/serial/simulation_base.cpp
)

target_compile_definitions(dgswemv2-serial PRIVATE ${LINALG_DEFINITION} ${PROBLEM_DEFINITIONS})
target_include_directories(dgswemv2-serial PRIVATE ${YAML_CPP_INCLUDE_DIR})
target_link_libraries(dgswemv2-serial ${YAML_CPP_LIBRARIES})
install(TARGETS dgswemv2-serial DESTINATION bin)

if(USE_OMPI)
  add_executable(
    dgswemv2-ompi
    ${PROJECT_SOURCE_DIR}/source/dgswemv2-ompi.cpp
    ${SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
    ${PROJECT_SOURCE_DIR}/source/simulation/ompi/simulation_ompi_base.cpp
  )

  target_compile_definitions(dgswemv2-ompi PRIVATE ${LINALG_DEFINITION} ${PROBLEM_DEFINITIONS})
  target_compile_options(dgswemv2-ompi PRIVATE ${OpenMP_CXX_FLAGS})
  target_include_directories(dgswemv2-ompi PRIVATE ${YAML_CPP_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(dgswemv2-ompi ${YAML_CPP_LIBRARIES} ${MPI_CXX_LIBRARIES} ${OpenMP_CXX_FLAGS})

  if(PETSC_FOUND)
    target_compile_definitions(dgswemv2-ompi PRIVATE HAS_PETSC)
    target_include_directories(dgswemv2-ompi PRIVATE ${PETSC_INCLUDES})
    target_link_libraries(dgswemv2-ompi ${PETSC_LIBRARIES})
  endif()

  install(TARGETS dgswemv2-ompi DESTINATION bin)
endif()

if(USE_HPX)
  add_executable(
    dgswemv2-hpx
    ${PROJECT_SOURCE_DIR}/source/dgswemv2-hpx.cpp
    ${SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/hpx_communicator.cpp
    ${PROJECT_SOURCE_DIR}/source/simulation/hpx/sim_unit_hpx_base.cpp
  )

  target_compile_definitions(dgswemv2-hpx PRIVATE HAS_HPX ${LINALG_DEFINITION} ${PROBLEM_DEFINITIONS})
  hpx_setup_target(
    dgswemv2-hpx
    COMPONENT_DEPENDENCIES iostreams
  )
  target_include_directories(dgswemv2-hpx PRIVATE ${YAML_CPP_INCLUDE_DIR} ${HPX_INCLUDE_DIRS})  
  target_link_libraries(dgswemv2-hpx ${YAML_CPP_LIBRARIES} ${HPX_LIBRARIES})
  install(TARGETS dgswemv2-hpx DESTINATION bin)
endif()