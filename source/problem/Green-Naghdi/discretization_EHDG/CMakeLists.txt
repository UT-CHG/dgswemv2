set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g")

add_executable(
  EHDG_GN_SERIAL
  ${PROJECT_SOURCE_DIR}/source/problem/Green-Naghdi/discretization_EHDG/main_files/ehdg_main_gn.cpp
  ${SOURCES}
  ${GN_SOURCES}
)

target_include_directories(EHDG_GN_SERIAL PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(EHDG_GN_SERIAL ${YAML_CPP_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS EHDG_GN_SERIAL DESTINATION bin)

#[[if(USE_OMPI)
  add_executable(
    EHDG_GN_OMPI
    ${PROJECT_SOURCE_DIR}/source/problem/GN/discretization_EHDG/main_files/ehdg_ompi_main_gn.cpp
    ${SOURCES}
    ${GN_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  set_target_properties(EHDG_GN_OMPI PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")

  target_include_directories(EHDG_GN_OMPI PRIVATE ${YAML_CPP_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(EHDG_GN_OMPI ${YAML_CPP_LIBRARIES} ${MPI_CXX_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${OpenMP_CXX_FLAGS})

  install(TARGETS EHDG_GN_OMPI DESTINATION bin)
endif()

if(USE_HPX)
  add_executable(
    EHDG_GN_HPX
    ${PROJECT_SOURCE_DIR}/source/problem/GN/discretization_EHDG/main_files/ehdg_hpx_main_gn.cpp
    ${SOURCES}
    ${GN_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/hpx_communicator.cpp
  )

  hpx_setup_target(
    EHDG_GN_HPX
    COMPONENT_DEPENDENCIES iostreams
  )

  target_include_directories(EHDG_GN_HPX PRIVATE ${YAML_CPP_INCLUDE_DIR} ${HPX_INCLUDE_DIRS})  
  target_link_libraries(EHDG_GN_HPX ${YAML_CPP_LIBRARIES} ${HPX_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

  install(TARGETS EHDG_GN_HPX DESTINATION bin)
endif()]]