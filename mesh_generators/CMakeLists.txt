add_executable(
  rectangular_mesh_generator
  ${PROJECT_SOURCE_DIR}/mesh_generators/rectangular_mesh_generator.cpp
)

set_target_properties(rectangular_mesh_generator PROPERTIES COMPILE_FLAGS ${COMPILE_FLGS})
install(TARGETS rectangular_mesh_generator DESTINATION bin/utilities/mesh_generators)
