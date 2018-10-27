add_executable(
  EHDG_SWE_SERIAL
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_EHDG/main_files/ehdg_main_swe.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

set_target_properties(EHDG_SWE_SERIAL PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
install(TARGETS EHDG_SWE_SERIAL DESTINATION bin)

if(USE_OMPI)
  add_executable(
    EHDG_SWE_OMPI
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_EHDG/main_files/ehdg_ompi_main_swe.cpp
    ${SOURCES}
    ${SWE_SOURCES}
    ${PROJECT_SOURCE_DIR}/source/communication/ompi_communicator.cpp
  )

  set_target_properties(EHDG_SWE_OMPI PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
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

  set_target_properties(EHDG_SWE_HPX PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
  hpx_setup_target(
    EHDG_SWE_HPX
    COMPONENT_DEPENDENCIES iostreams
  )
  install(TARGETS EHDG_SWE_HPX DESTINATION bin)
endif()