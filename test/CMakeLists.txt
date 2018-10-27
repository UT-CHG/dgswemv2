add_executable(
  test_linear_algebra_exe
  test_linear_algebra.cpp
)

set_target_properties(test_linear_algebra_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

#[[target_include_directories(test_linear_algebra_exe PRIVATE ${PETSC_INCLUDES})  
target_link_libraries(test_linear_algebra_exe ${PETSC_LIBRARIES})

add_executable(
  test_jacobians_ehdg_exe
  test_jacobians_ehdg.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

target_include_directories(test_jacobians_ehdg_exe PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(test_jacobians_ehdg_exe ${YAML_CPP_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_test(
  Unit_test_jacobians_ehdg
  test_jacobians_ehdg_exe
)

add_executable(
  test_jacobians_ihdg_exe
  test_jacobians_ihdg.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

target_include_directories(test_jacobians_ihdg_exe PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(test_jacobians_ihdg_exe ${YAML_CPP_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_test(
  Unit_test_jacobians_ihdg
  test_jacobians_ihdg_exe
)

add_executable(
  test_global_jacobian_ihdg_exe
  test_global_jacobian_ihdg.cpp
  ${SOURCES}
  ${SWE_SOURCES}
)

target_include_directories(test_global_jacobian_ihdg_exe PRIVATE ${YAML_CPP_INCLUDE_DIR})  
target_link_libraries(test_global_jacobian_ihdg_exe ${YAML_CPP_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

add_test(
  Unit_test_global_jacobian_ihdg
  test_global_jacobian_ihdg_exe
)]]

add_executable(
  test_tuple_helpers_exe
  test_tuple_helpers.cpp
)

set_target_properties(test_tuple_helpers_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

add_test(
  Unit_tuple_helpers
  test_tuple_helpers_exe
)

add_executable(
  test_heterogeneous_containers_exe
  test_heterogeneous_containers.cpp
)

set_target_properties(test_heterogeneous_containers_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

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

set_target_properties(test_adcirc_format_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

add_test(
  Unit_adcirc_format
  test_adcirc_format_exe
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/sample_fort.14
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/weir/weir.14
)

add_executable(
  test_basis_legendre_1d_exe
  test_basis_legendre_1d.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_1D/basis_legendre_1D.cpp
)

set_target_properties(test_basis_legendre_1d_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

add_test(
  Unit_basis_legendre_1d
  test_basis_legendre_1d_exe
)

add_executable(
  test_basis_dubiner_2d_exe
  test_basis_dubiner_2d.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp
)

set_target_properties(test_basis_dubiner_2d_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

add_test(
  Unit_basis_dubiner_2d
  test_basis_dubiner_2d_exe
)

add_executable(
  test_integ_gausslegendre_1d_exe
  test_integ_gausslegendre_1d.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
)

set_target_properties(test_integ_gausslegendre_1d_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

add_test(
  Unit_integ_gausslegendre_1d
  test_integ_gausslegendre_1d_exe
)

add_executable(
  test_integ_dunavant_2d_exe
  test_integ_dunavant_2d.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
)

set_target_properties(test_integ_dunavant_2d_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

add_test(
  Unit_integ_dunavant_2d
  test_integ_dunavant_2d_exe
)

add_executable(
  test_shape_straight_triangle_exe
  test_shape_straight_triangle.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
)

set_target_properties(test_shape_straight_triangle_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

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

set_target_properties(test_element_triangle_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_link_libraries(test_element_triangle_exe ${YAML_CPP_LIBRARIES})

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

set_target_properties(test_boundary_interface_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_link_libraries(test_boundary_interface_exe ${YAML_CPP_LIBRARIES})

add_test(
  Unit_boundary_interface
  test_boundary_interface_exe
)

add_executable(
  test_edges_exe
  test_edges.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_1D/basis_legendre_1D.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
)

set_target_properties(test_edges_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_link_libraries(test_edges_exe ${YAML_CPP_LIBRARIES})

add_test(
  Unit_edges
  test_edges_exe
)

add_executable(
  test_rk_stepper_exe
  test_rk_stepper.cpp
  ${PROJECT_SOURCE_DIR}/source/simulation/stepper/rk_stepper.cpp
)

set_target_properties(test_rk_stepper_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_link_libraries(test_rk_stepper_exe ${YAML_CPP_LIBRARIES})

add_test(
  Unit_RK_stepper
  test_rk_stepper_exe
)

add_executable(
  test_llf_flux_exe
  test_llf_flux.cpp
)

set_target_properties(test_llf_flux_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

add_test(
  Unit_llf_flux
  test_llf_flux_exe
)

add_executable(
  test_wetting_drying_exe
  test_wetting_drying.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
  ${PROJECT_SOURCE_DIR}/source/simulation/stepper/rk_stepper.cpp
)

set_target_properties(test_wetting_drying_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_link_libraries(test_wetting_drying_exe ${YAML_CPP_LIBRARIES})

add_test(
  Unit_wetting_drying
  test_wetting_drying_exe
)

add_executable(
  test_mesh_metadata_exe
  test_mesh_metadata.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
)

set_target_properties(test_mesh_metadata_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

add_test(
  Unit_mesh_metadata
  test_mesh_metadata_exe
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/sample_fort.14
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/weir/weir.14
)

add_executable(
  test_swe_inputs_exe
  test_swe_inputs.cpp
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_input/swe_inputs.cpp
)

set_target_properties(test_swe_inputs_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_link_libraries(test_swe_inputs_exe ${YAML_CPP_LIBRARIES})

add_test(
  Unit_swe_inputs
  test_swe_inputs_exe
)

add_executable(
  test_input_parameters_exe
  test_input_parameters.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_input/swe_inputs.cpp
)

set_target_properties(test_input_parameters_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_link_libraries(test_input_parameters_exe ${YAML_CPP_LIBRARIES})

add_test(
  Unit_input_parameters
  test_input_parameters_exe
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/correct_input.15
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/input_that_doesnt_exist.15
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/missing_field_input.15
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/correct_input_no_output.15

)

#Everytime Cmake is run, we generate an input file with correct build parameters
file(WRITE files_for_testing/correct_input.15
"mesh:\n\
  format: Adcirc\n\
  file_name: ${PROJECT_SOURCE_DIR}/test/files_for_testing/sample_fort.14\n\
  coordinate_system: cartesian\n\n\
timestepping:\n\
  start_time: 11-05-2015 12:00              #dd-mm-yyyy hh:mm\n\
  end_time: 12-05-2015 00:00                #dd-mm-yyyy hh:mm\n\
  dt: 1                           #in seconds\n\
  order: 2\n\
  nstages: 2\n\
  ramp_duration: 1.5\n\n\
polynomial_order: 2\n\n\
problem:\n\
  name: swe\n\
  gravity: 9.81\n\
  density_air: 1.225\n\
  density_water: 1000\n\
  h_o: 0.001\n\
  initial_conditions:\n\
    type: Constant\n\
    initial_surface_height: 10.0\n\
    initial_momentum_x: 0.0\n\
    initial_momentum_y: 0.0\n\
  function_source: Enable\n\
  bottom_friction:\n\
    type: Manning\n\
    coefficient: 0.001\n\
    input_file: manning.test\n\
  meteo_forcing:\n\
    type: HWIND\n\
    raw_input_file: weather.test\n\
    input_file: weather.test\n\
    frequency: 900\n\
  tide_potential:\n\
    type: Test\n\
  coriolis: Enable\n\
  wetting_drying:\n\
    h_o: 0.1\n\
  slope_limiting:\n\
    type: Cockburn-Shu\n\
    M: 50\n\
    nu: 1\n\n\
output:\n\
  logfile:\n\
    name: log\n\
    verbose: true\n\
  path: output/\n\
  vtk:\n\
    frequency: 100\n\n\
load_balancer:\n\
  name: random\n\
  rebalance_frequency: 10"
)

file(WRITE files_for_testing/correct_input_no_output.15
"mesh:\n\
  format: Adcirc\n\
  file_name: ${PROJECT_SOURCE_DIR}/test/files_for_testing/sample_fort.14\n\
  coordinate_system: cartesian\n\n\
timestepping:\n\
  start_time: 11-05-2015 12:00              #dd-mm-yyyy hh:mm\n\
  end_time: 12-05-2015 00:00                #dd-mm-yyyy hh:mm\n\
  dt: 1                           #in seconds\n\
  order: 2\n\
  nstages: 2\n\n\
polynomial_order: 2\n\n\
problem:\n\
  name: swe\n\
  gravity: 9.81\n\
  density_air: 1.225\n\
  density_water: 1000\n\
  h_o: 0.001\n\
  initial_conditions:\n\
    type: Constant\n\
    initial_surface_height: 10.0\n\
    initial_momentum_x: 0.0\n\
    initial_momentum_y: 0.0\n\
  function_source: Enable\n\
  bottom_friction:\n\
    type: Manning\n\
    coefficient: 0.001\n\
    input_file: manning.test\n\
  meteo_forcing:\n\
    type: HWIND\n\
    raw_input_file: weather.test\n\
    input_file: weather.test\n\
    frequency: 900\n\
  tide_potential:\n\
    type: Test\n\
  coriolis: Enable\n\
  wetting_drying:\n\
    h_o: 0.1\n\
  slope_limiting:\n\
    type: Cockburn-Shu\n\
    M: 50\n\
    nu: 1\n"
)

add_executable(
  test_partition_exe
  test_partition.cpp
  ${PROJECT_SOURCE_DIR}/partitioner/partition.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
)

set_target_properties(test_partition_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})
target_include_directories(test_partition_exe PRIVATE ${METIS_INCLUDE_DIRS})
target_link_libraries(test_partition_exe ${METIS_LIBRARIES})

add_test(
  Unit_partition
  test_partition_exe
  ${PROJECT_SOURCE_DIR}/test/files_for_testing/sample_fort.14
)

add_executable(
  test_heartbeat_exe
  test_heartbeat.cpp
)

set_target_properties(test_heartbeat_exe PROPERTIES COMPILE_FLAGS ${LINALG_FLAG})

add_test(
  Unit_heartbeat
  test_heartbeat_exe
)

if(USE_HPX)
  #[[add_executable(
    test_rkdg_swe_serialization_exe
    test_rkdg_swe_serialization.cpp
  )

  set_target_properties(test_rkdg_swe_serialization_exe PROPERTIES COMPILE_FLAGS "-DHAS_HPX")

  target_include_directories(test_rkdg_swe_serialization_exe PRIVATE ${HPX_INCLUDE_DIRS})
  target_link_libraries(test_rkdg_swe_serialization_exe ${HPX_LIBRARIES})

  add_test(
    Unit_rkdg_swe_serialization
    test_rkdg_swe_serialization_exe
  )]]

  add_executable(
    test_shape_serialization_exe
    test_shape_serialization.cpp
    ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
  )

  set_target_properties(test_shape_serialization_exe PROPERTIES COMPILE_FLAGS "-DHAS_HPX ${LINALG_FLAG}")
  target_include_directories(test_shape_serialization_exe PRIVATE ${HPX_INCLUDE_DIRS})
  target_link_libraries(test_shape_serialization_exe ${HPX_LIBRARIES} ${YAML_CPP_LIBRARIES})

  add_test(
    Unit_shape_serialization
    test_shape_serialization_exe
  )

  add_executable(
    test_element_serialization_exe
    test_element_serialization.cpp
    ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
    ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp
    ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
    ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
    ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
  )

  set_target_properties(test_element_serialization_exe PROPERTIES COMPILE_FLAGS "-DHAS_HPX ${LINALG_FLAG}")
  target_include_directories(test_element_serialization_exe PRIVATE ${HPX_INCLUDE_DIRS})
  target_link_libraries(test_element_serialization_exe ${HPX_LIBRARIES} ${YAML_CPP_LIBRARIES})

  add_test(
    Unit_element_serialization
    test_element_serialization_exe
  )

  #[[add_executable(
    test_mesh_serialization_exe
    test_mesh_serialization.cpp
    ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
    ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp
    ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
    ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp
    ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
    ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
    ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_input/swe_inputs.cpp
  )

  set_target_properties(test_mesh_serialization_exe PROPERTIES COMPILE_FLAGS "-DHAS_HPX")

  target_include_directories(test_mesh_serialization_exe PRIVATE ${HPX_INCLUDE_DIRS})
  target_link_libraries(test_mesh_serialization_exe ${HPX_LIBRARIES} ${YAML_CPP_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

  add_test(
    Unit_mesh_serialization
    test_mesh_serialization_exe
    ${PROJECT_SOURCE_DIR}/test/files_for_testing/correct_input.15
  )]]

  add_executable(
    test_rk_stepper_serialization_exe
    test_rk_stepper_serialization.cpp
    ${PROJECT_SOURCE_DIR}/source/simulation/stepper/rk_stepper.cpp
  )

  set_target_properties(test_rk_stepper_serialization_exe PROPERTIES COMPILE_FLAGS "-DHAS_HPX ${LINALG_FLAG}")
  target_include_directories(test_rk_stepper_serialization_exe PRIVATE ${HPX_INCLUDE_DIRS})
  target_link_libraries(test_rk_stepper_serialization_exe ${HPX_LIBRARIES} ${YAML_CPP_LIBRARIES})

  
  add_test(
    Unit_rk_stepper_serialization
    test_rk_stepper_serialization_exe
  )

  add_executable(
    test_serialize_chrono_exe
    test_serialize_chrono.cpp
  )

  set_target_properties(test_serialize_chrono_exe PROPERTIES COMPILE_FLAGS "-DHAS_HPX ${LINALG_FLAG}")
  target_include_directories(test_serialize_chrono_exe PRIVATE ${HPX_INCLUDE_DIRS})
  target_link_libraries(test_serialize_chrono_exe ${HPX_LIBRARIES} ${YAML_CPP_LIBRARIES})

  add_test(
    Unit_serialize_chrono
    test_serialize_chrono_exe
  )

  if(USE_BLAZE)
    add_executable(
      test_blaze_serialization_exe
      test_blaze_serialization.cpp
    )

    set_target_properties(test_blaze_serialization_exe PROPERTIES COMPILE_FLAGS "-DHAS_HPX -DUSE_BLAZE")
    target_include_directories(test_blaze_serialization_exe PRIVATE ${HPX_INCLUDE_DIRS})
    target_link_libraries(test_blaze_serialization_exe ${HPX_LIBRARIES} ${YAML_CPP_LIBRARIES})
  
    add_test(
      Unit_blaze_serialization
      test_blaze_serialization_exe
    )
  endif()

  if(USE_EIGEN)
    add_executable(
      test_eigen_serialization_exe
      test_eigen_serialization.cpp
    )

    set_target_properties(test_eigen_serialization_exe PROPERTIES COMPILE_FLAGS "-DHAS_HPX -DUSE_EIGEN")
    target_include_directories(test_eigen_serialization_exe PRIVATE ${HPX_INCLUDE_DIRS})
    target_link_libraries(test_eigen_serialization_exe ${HPX_LIBRARIES} ${YAML_CPP_LIBRARIES})
  
    add_test(
      Unit_eigen_serialization
      test_eigen_serialization_exe
    )
  endif()
endif()
