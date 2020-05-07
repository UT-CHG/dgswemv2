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