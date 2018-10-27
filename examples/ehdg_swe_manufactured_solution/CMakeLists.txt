add_executable(
  EHDG_MANUFACTURED_SOLUTION_SERIAL
  ${PROJECT_SOURCE_DIR}/examples/ehdg_swe_manufactured_solution/manufactured_main_swe.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

set_target_properties(EHDG_MANUFACTURED_SOLUTION_SERIAL PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
install(TARGETS EHDG_MANUFACTURED_SOLUTION_SERIAL DESTINATION bin)

if(USE_OMPI)
  add_executable(
    EHDG_MANUFACTURED_SOLUTION_OMPI
    ${PROJECT_SOURCE_DIR}/examples/ehdg_swe_manufactured_solution/manufactured_ompi_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  set_target_properties(EHDG_MANUFACTURED_SOLUTION_OMPI PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
  install(TARGETS EHDG_MANUFACTURED_SOLUTION_OMPI DESTINATION bin)
endif()

if(USE_HPX)
  add_executable(
    EHDG_MANUFACTURED_SOLUTION_HPX
    ${PROJECT_SOURCE_DIR}/examples/ehdg_swe_manufactured_solution/manufactured_hpx_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/hpx_communicator.cpp
  )

  set_target_properties(EHDG_MANUFACTURED_SOLUTION_HPX PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
  hpx_setup_target(
    EHDG_MANUFACTURED_SOLUTION_HPX
    COMPONENT_DEPENDENCIES iostreams
  )
  install(TARGETS EHDG_MANUFACTURED_SOLUTION_HPX DESTINATION bin)
endif()