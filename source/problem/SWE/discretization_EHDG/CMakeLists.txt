add_executable(
  EHDG_SWE_SERIAL
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_EHDG/main_files/ehdg_main_swe.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

set_target_properties(EHDG_SWE_SERIAL PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_include_directories(EHDG_SWE_SERIAL PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(EHDG_SWE_SERIAL ${YAML_CPP_LIBRARIES})
install(TARGETS EHDG_SWE_SERIAL DESTINATION bin)

if(USE_OMPI)
  add_executable(
    EHDG_SWE_OMPI
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_EHDG/main_files/ehdg_ompi_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  set_target_properties(EHDG_SWE_OMPI PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS} ${LINALG_FLAG}")
  target_include_directories(EHDG_SWE_OMPI PRIVATE ${YAML_CPP_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(EHDG_SWE_OMPI ${YAML_CPP_LIBRARIES} ${MPI_CXX_LIBRARIES} ${OpenMP_CXX_FLAGS})
  install(TARGETS EHDG_SWE_OMPI DESTINATION bin)
endif()

if(USE_HPX)
  add_executable(
    EHDG_SWE_HPX
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_EHDG/main_files/ehdg_hpx_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/hpx_communicator.cpp
  )

  set_target_properties(EHDG_SWE_HPX PROPERTIES COMPILE_FLAGS "-DHAS_HPX ${LINALG_FLAG}")
  hpx_setup_target(
    EHDG_SWE_HPX
    COMPONENT_DEPENDENCIES iostreams
  )
  target_include_directories(EHDG_SWE_HPX PRIVATE ${YAML_CPP_INCLUDE_DIR} ${HPX_INCLUDE_DIRS})  
  target_link_libraries(EHDG_SWE_HPX ${YAML_CPP_LIBRARIES} ${HPX_LIBRARIES})
  install(TARGETS EHDG_SWE_HPX DESTINATION bin)
endif()