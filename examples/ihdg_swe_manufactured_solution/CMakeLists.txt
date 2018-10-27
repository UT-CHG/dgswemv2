add_executable(
  IHDG_MANUFACTURED_SOLUTION_SERIAL
  ${PROJECT_SOURCE_DIR}/examples/ihdg_swe_manufactured_solution/manufactured_main_swe.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

set_target_properties(IHDG_MANUFACTURED_SOLUTION_SERIAL PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
install(TARGETS IHDG_MANUFACTURED_SOLUTION_SERIAL DESTINATION bin)

if(USE_OMPI)
  add_executable(
    IHDG_MANUFACTURED_SOLUTION_OMPI
    ${PROJECT_SOURCE_DIR}/examples/ihdg_swe_manufactured_solution/manufactured_ompi_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  set_target_properties(IHDG_MANUFACTURED_SOLUTION_OMPI PROPERTIES COMPILE_FLAGS "${COMPILE_FLGS} -DHAS_PETSC")
  install(TARGETS IHDG_MANUFACTURED_SOLUTION_OMPI DESTINATION bin)
endif()
