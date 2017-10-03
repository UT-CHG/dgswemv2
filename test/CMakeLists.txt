set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall")

add_executable(
  test_tuple_helpers_exe
  test_tuple_helpers.cpp
)

add_test(
  Unit_tuple_helpers
  test_tuple_helpers_exe
)

add_executable(
  test_heterogeneous_containers_exe
  test_heterogeneous_containers.cpp
)

add_test(
  Unit_heterogeneous_containers
  test_heterogeneous_containers_exe
)

add_executable(
  test_adcirc_format_exe
  test_adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
)

add_test(
  Unit_adcirc_format
  test_adcirc_format_exe
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/sample_fort.14
)

add_executable(
  test_basis_dubiner_2d_exe
  test_basis_dubiner_2d.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp
)

add_test(
  Unit_basis_dubiner_2d
  test_basis_dubiner_2d_exe
)

add_executable(
  test_integ_gausslegendre_1d_exe
  test_integ_gausslegendre_1d.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
)

add_test(
  Unit_integ_gausslegendre_1d
  test_integ_gausslegendre_1d_exe
)

add_executable(
  test_integ_dunavant_2d_exe
  test_integ_dunavant_2d.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
)

add_test(
  Unit_integ_dunavant_2d
  test_integ_dunavant_2d_exe
)

add_executable(
  test_shape_straight_triangle_exe
  test_shape_straight_triangle.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
)

add_test(
  Unit_shape_straight_triangle
  test_shape_straight_triangle_exe
)

add_executable(
  test_element_triangle_exe
  test_element_triangle.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
)

add_test(
  Unit_element_triangle
  test_element_triangle_exe
)

add_executable(
  test_boundary_interface_exe
  test_boundary_interface.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
)

add_test(
  Unit_boundary_interface
  test_boundary_interface_exe
)

add_executable(
  test_rk_stepper_exe
  test_rk_stepper.cpp
  ${PROJECT_SOURCE_DIR}/source/simulation/stepper.cpp
)

add_test(
  Unit_RK_stepper
  test_rk_stepper_exe
)

#add_executable(
#  test_triangulation_exe
#  test_triangulation.cpp
#  ${PROJECT_SOURCE_DIR}/src/writer/triangulation.cpp
#)

#add_test(
#  Unit_triangulation
#  test_triangulation_exe
#)

add_executable(
  test_llf_flux_exe
  test_llf_flux.cpp
)

add_test(
  Unit_llf_flux
  test_llf_flux_exe
)

#add_executable(
#  test_swe_boundaries_exe
#  test_swe_boundaries.cpp
#)

#add_test(
#  Uint_test_land_boundary.cpp
#  test_swe_boundaries_exe
#)

add_executable(
  test_mesh_metadata_exe
  test_mesh_metadata.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
)

add_test(
  Unit_mesh_metadata
  test_mesh_metadata_exe
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/sample_fort.14
)

include_directories("${YAML_CPP_INCLUDE_DIR}")

add_executable(
  test_input_parameters_exe
  test_input_parameters.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/input_parameters.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/swe_inputs.cpp
)

target_link_libraries(test_input_parameters_exe ${YAML_CPP_LIBRARIES})

add_test(
  Unit_input_parameters
  test_input_parameters_exe
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/correct_input.15
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/input_that_doesnt_exist.15
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/missing_field_input.15
)

#Everytime Cmake is run, we generate an input file with correct build parameters
file(WRITE files_for_testing/correct_input.15
"mesh:\n\
  format: Adcirc\n\
  file_name: ${PROJECT_SOURCE_DIR}/test/files_for_testing/sample_fort.14\n\n\
timestepping:\n\
  dt: 1                           #in seconds\n\
  end_time: 86400                 #in seconds\n\
  order: 2\n\
  nstages: 2\n\n\
polynomial_order: 2"
)

add_executable(
  test_partition_exe
  test_partition.cpp
  ${PROJECT_SOURCE_DIR}/partitioner/partition.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
)

target_include_directories(test_partition_exe PRIVATE ${METIS_INCLUDE_DIRS})
target_link_libraries(test_partition_exe ${METIS_LIBRARIES})

add_test(
  Unit_partition
  test_partition_exe
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/sample_fort.14
)
