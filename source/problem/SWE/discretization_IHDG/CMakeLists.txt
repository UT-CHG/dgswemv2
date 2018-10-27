add_executable(
  IHDG_SWE_SERIAL
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/discretization_IHDG/main_files/ihdg_main_swe.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

set_target_properties(IHDG_SWE_SERIAL PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_include_directories(IHDG_SWE_SERIAL PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(IHDG_SWE_SERIAL ${YAML_CPP_LIBRARIES} )
install(TARGETS IHDG_SWE_SERIAL DESTINATION bin)