if (SWE)
  list(APPEND SOURCES
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_input/swe_inputs.cpp
  )
  list (APPEND PROBLEM_FLAGS "${PROBLEM_FLAGS} -DBUILD_SWE")
endif()

if (GN)
  list(APPEND SOURCES
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_input/swe_inputs.cpp
  )
  list (APPEND PROBLEM_FLAGS "${PROBLEM_FLAGS} -DBUILD_SWE")
endif()

#add_executable(
#  RKDG_SWE_SERIAL
#  ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_RKDG/main_files/rkdg_main_swe.cpp
#  ${SOURCES}
#  ${SWE_SOURCES}
#)

#set_target_properties(RKDG_SWE_SERIAL PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
#target_include_directories(RKDG_SWE_SERIAL PRIVATE ${YAML_CPP_INCLUDE_DIR})  
#target_link_libraries(RKDG_SWE_SERIAL ${YAML_CPP_LIBRARIES})
#install(TARGETS RKDG_SWE_SERIAL DESTINATION bin)

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