add_executable(
  quad_mesh_generator
  ${PROJECT_SOURCE_DIR}/mesh_generators/quad_mesh_generator.cpp
)

target_compile_definitions(quad_mesh_generator PRIVATE ${LINALG_DEFINITION})
target_include_directories(quad_mesh_generator PRIVATE ${YAML_CPP_INCLUDE_DIR})
target_link_libraries(quad_mesh_generator ${YAML_CPP_LIBRARIES})
install(TARGETS quad_mesh_generator DESTINATION bin/utilities/mesh_generators)
