add_executable(
  IHDG_SWE_SERIAL
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_IHDG/main_files/ihdg_main_swe.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

set_target_properties(IHDG_SWE_SERIAL PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
install(TARGETS IHDG_SWE_SERIAL DESTINATION bin)

#[[if(USE_OMPI)
  add_executable(
    IHDG_SWE_OMPI
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_IHDG/main_files/ihdg_ompi_main_swe.cpp
    ${IHDG_SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  set_target_properties(IHDG_SWE_OMPI PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")

  target_include_directories(IHDG_SWE_OMPI PRIVATE ${YAML_CPP_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH})  
  target_link_libraries(IHDG_SWE_OMPI ${YAML_CPP_LIBRARIES} ${MPI_CXX_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${OpenMP_CXX_FLAGS})

  install(TARGETS IHDG_SWE_OMPI DESTINATION bin)
endif()
]]