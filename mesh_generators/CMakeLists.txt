add_executable(
  rectangular_mesh_generator
  ${PROJECT_SOURCE_DIR}/mesh_generators/rectangular_mesh_generator.cpp
)

target_compile_definitions(rectangular_mesh_generator PRIVATE ${LINALG_DEFINITION})
target_include_directories(rectangular_mesh_generator PRIVATE ${YAML_CPP_INCLUDE_DIR})
target_link_libraries(rectangular_mesh_generator ${YAML_CPP_LIBRARIES})
install(TARGETS rectangular_mesh_generator DESTINATION bin/utilities/mesh_generators)
