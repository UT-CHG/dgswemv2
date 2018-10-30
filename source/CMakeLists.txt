add_executable(
  dgswemv2-serial
  ${PROJECT_SOURCE_DIR}/source/dgswemv2-serial.cpp
  ${SOURCES}
  ${PROJECT_SOURCE_DIR}/source/simulation/serial/simulation_base.cpp
)

set_target_properties(dgswemv2-serial PROPERTIES COMPILE_FLAGS "${LINALG_FLAG} ${PROBLEM_FLAGS}")
target_include_directories(dgswemv2-serial PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(dgswemv2-serial ${YAML_CPP_LIBRARIES})
install(TARGETS dgswemv2-serial DESTINATION bin)

#if(USE_OMPI)
#  add_executable(
#    RKDG_SWE_OMPI
#    ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_RKDG/main_files/rkdg_ompi_main_swe.cpp
#    ${SOURCES}
#    ${SWE_SOURCES}
#    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
#  )

#  set_target_properties(RKDG_SWE_OMPI PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${LINALG_FLAG}")
#  target_include_directories(RKDG_SWE_OMPI PRIVATE ${YAML_CPP_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH})
#  target_link_libraries(RKDG_SWE_OMPI ${YAML_CPP_LIBRARIES} ${MPI_CXX_LIBRARIES} ${OpenMP_CXX_FLAGS})
#  install(TARGETS RKDG_SWE_OMPI DESTINATION bin)
#endif()

if(USE_HPX)
  add_executable(
    dgswemv2-hpx
    ${PROJECT_SOURCE_DIR}/source/dgswemv2-hpx.cpp
    ${SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/hpx_communicator.cpp
    ${PROJECT_SOURCE_DIR}/source/simulation/hpx/sim_unit_hpx_base.cpp
  )

  set_target_properties(dgswemv2-hpx PROPERTIES COMPILE_FLAGS "-DHAS_HPX ${LINALG_FLAG} ${PROBLEM_FLAGS}")
  hpx_setup_target(
    dgswemv2-hpx
    COMPONENT_DEPENDENCIES iostreams
  )
  target_include_directories(dgswemv2-hpx PRIVATE ${YAML_CPP_INCLUDE_DIR} ${HPX_INCLUDE_DIRS})  
  target_link_libraries(dgswemv2-hpx ${YAML_CPP_LIBRARIES} ${HPX_LIBRARIES})
  install(TARGETS dgswemv2-hpx DESTINATION bin)
endif()