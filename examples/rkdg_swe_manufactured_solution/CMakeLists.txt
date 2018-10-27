add_executable(
  RKDG_MANUFACTURED_SOLUTION_SERIAL
  ${PROJECT_SOURCE_DIR}/examples/rkdg_swe_manufactured_solution/manufactured_main_swe.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

set_target_properties(RKDG_MANUFACTURED_SOLUTION_SERIAL PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_include_directories(RKDG_MANUFACTURED_SOLUTION_SERIAL PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(RKDG_MANUFACTURED_SOLUTION_SERIAL ${YAML_CPP_LIBRARIES})
install(TARGETS RKDG_MANUFACTURED_SOLUTION_SERIAL DESTINATION bin)

if(USE_OMPI)
  add_executable(
    RKDG_MANUFACTURED_SOLUTION_OMPI
    ${PROJECT_SOURCE_DIR}/examples/rkdg_swe_manufactured_solution/manufactured_ompi_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  set_target_properties(RKDG_MANUFACTURED_SOLUTION_OMPI PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${LINALG_FLAG}")
  target_include_directories(RKDG_MANUFACTURED_SOLUTION_OMPI PRIVATE ${YAML_CPP_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(RKDG_MANUFACTURED_SOLUTION_OMPI ${YAML_CPP_LIBRARIES} ${MPI_CXX_LIBRARIES} ${OpenMP_CXX_FLAGS})
  install(TARGETS RKDG_MANUFACTURED_SOLUTION_OMPI DESTINATION bin)
endif()

if(USE_HPX)
  add_executable(
    RKDG_MANUFACTURED_SOLUTION_HPX
    ${PROJECT_SOURCE_DIR}/examples/rkdg_swe_manufactured_solution/manufactured_hpx_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}  
    ${PROJECT_SOURCE_DIR}/source/communication/hpx_communicator.cpp
  )

  set_target_properties(RKDG_MANUFACTURED_SOLUTION_HPX PROPERTIES COMPILE_FLAGS "-DHAS_HPX ${LINALG_FLAG}")
  hpx_setup_target(
    RKDG_MANUFACTURED_SOLUTION_HPX
    COMPONENT_DEPENDENCIES iostreams
  )
  target_include_directories(RKDG_MANUFACTURED_SOLUTION_HPX PRIVATE ${YAML_CPP_INCLUDE_DIR} ${HPX_INCLUDE_DIRS})
  target_link_libraries(RKDG_MANUFACTURED_SOLUTION_HPX ${YAML_CPP_LIBRARIES} ${HPX_LIBRARIES})
  install(TARGETS RKDG_MANUFACTURED_SOLUTION_HPX DESTINATION bin)
endif()