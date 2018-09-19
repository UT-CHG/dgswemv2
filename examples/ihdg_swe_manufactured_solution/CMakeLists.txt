set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -g")

add_executable(
  IHDG_MANUFACTURED_SOLUTION_SERIAL
  ${PROJECT_SOURCE_DIR}/examples/ihdg_swe_manufactured_solution/manufactured_main_swe.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

target_include_directories(IHDG_MANUFACTURED_SOLUTION_SERIAL PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(IHDG_MANUFACTURED_SOLUTION_SERIAL ${YAML_CPP_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

if(USE_OMPI)
  add_executable(
    IHDG_MANUFACTURED_SOLUTION_OMPI
    ${PROJECT_SOURCE_DIR}/examples/ihdg_swe_manufactured_solution/manufactured_ompi_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  set_target_properties(IHDG_MANUFACTURED_SOLUTION_OMPI PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")

  target_include_directories(IHDG_MANUFACTURED_SOLUTION_OMPI PRIVATE ${YAML_CPP_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH} ${PETSC_INCLUDES})  
  target_link_libraries(IHDG_MANUFACTURED_SOLUTION_OMPI ${YAML_CPP_LIBRARIES} ${MPI_CXX_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${OpenMP_CXX_FLAGS} ${PETSC_LIBRARIES})

  install(TARGETS IHDG_MANUFACTURED_SOLUTION_OMPI DESTINATION bin)
endif()

#[[if(USE_HPX)
  add_executable(
    IHDG_MANUFACTURED_SOLUTION_HPX
    ${PROJECT_SOURCE_DIR}/examples/ihdg_swe_manufactured_solution/manufactured_hpx_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/hpx_communicator.cpp
  )

  hpx_setup_target(
    IHDG_MANUFACTURED_SOLUTION_HPX
    COMPONENT_DEPENDENCIES iostreams
  )

  target_include_directories(IHDG_MANUFACTURED_SOLUTION_HPX PRIVATE ${YAML_CPP_INCLUDE_DIR} ${HPX_INCLUDE_DIRS})  
  target_link_libraries(IHDG_MANUFACTURED_SOLUTION_HPX ${YAML_CPP_LIBRARIES} ${HPX_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

  install(TARGETS IHDG_MANUFACTURED_SOLUTION_HPX DESTINATION bin)
endif()

install(TARGETS IHDG_MANUFACTURED_SOLUTION_SERIAL DESTINATION bin)]]
