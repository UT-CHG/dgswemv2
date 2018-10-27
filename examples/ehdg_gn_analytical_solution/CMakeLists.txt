add_executable(
  EHDG_GN_ANALYTICAL_SOLUTION_SERIAL
  ${PROJECT_SOURCE_DIR}/examples/ehdg_gn_analytical_solution/analytical_main_gn.cpp
  ${SOURCES}
  ${GN_SOURCES}
)

set_target_properties(EHDG_GN_ANALYTICAL_SOLUTION_SERIAL PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
install(TARGETS EHDG_GN_ANALYTICAL_SOLUTION_SERIAL DESTINATION bin)

if(USE_OMPI)
  add_executable(
    EHDG_GN_ANALYTICAL_SOLUTION_OMPI
    ${PROJECT_SOURCE_DIR}/examples/ehdg_gn_analytical_solution/analytical_ompi_main_gn.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  set_target_properties(EHDG_GN_ANALYTICAL_SOLUTION_OMPI PROPERTIES COMPILE_FLAGS "${COMPILE_FLGS} -DHAS_PETSC")
  install(TARGETS EHDG_GN_ANALYTICAL_SOLUTION_OMPI DESTINATION bin)
endif()
