cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project (DG_HYPER CXX C)

option(USE_OMPI "Use MPI+OpenMP" OFF)
option(USE_HPX "Use HPX" OFF)
option(COMPILER_WARNINGS "Enable Compiler Warnings" OFF)
option(SET_VERBOSE "Verbose MakeFile" OFF)
option(BUILD_EXAMPLES "Build Examples" OFF)

option(RKDG_SWE "Build RKDG SWE simulations" ON)
option(EHDG_SWE "Build EHDG SWE simulations" OFF)
option(IHDG_SWE "Build EHDG SWE simulations" OFF)

option(EHDG_GN "Build EHDG GN simulations" OFF)
option(IHDG_GN "Build IHDG GN simulations" OFF)

enable_testing()

find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # require at least intel 17.0
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17.0)
    message(FATAL_ERROR "Intel version must be at least 17.0")
  endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # require at least 5.2
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.2)
    message(FATAL_ERROR "GCC version must be at least 5.2")
  endif()
else()
  message(WARNING "Using untested compiler")
endif()

if(COMPILER_WARNINGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

set(CMAKE_VERBOSE_MAKEFILE OFF)
if(SET_VERBOSE)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

get_filename_component (default_prefix "../install" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in install."
FORCE)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

find_package(METIS REQUIRED)
find_package(yaml-cpp REQUIRED)

find_package(Eigen3 REQUIRED)
find_package(blaze REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(USE_OMPI)
  find_package(OpenMP REQUIRED)
  find_package(MPI REQUIRED)
endif()

if(USE_HPX)
  find_package(HPX REQUIRED)
endif()

list (APPEND SOURCES
  ${PROJECT_SOURCE_DIR}/source/preprocessor/ADCIRC_reader/adcirc_format.cpp
  ${PROJECT_SOURCE_DIR}/source/preprocessor/mesh_metadata.cpp

  ${PROJECT_SOURCE_DIR}/source/simulation/stepper/rk_stepper.cpp

  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_1D/basis_legendre_1D.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp

  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp

  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
)

list(APPEND SWE_SOURCES
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_input/swe_inputs.cpp
  ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_parser/swe_parser.cpp
)

list(APPEND GN_SOURCES
${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_input/swe_inputs.cpp
${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_parser/swe_parser.cpp
)

include_directories(source)

add_subdirectory(test)
add_subdirectory(mesh_generators)
add_subdirectory(partitioner)

if(RKDG_SWE)
  add_subdirectory(source/problem/SWE/discretization_RKDG)
  if(BUILD_EXAMPLES)
    add_subdirectory(examples/rkdg_swe_manufactured_solution)
  endif()
endif()

if(EHDG_SWE)
  add_subdirectory(source/problem/SWE/discretization_EHDG)
  if(BUILD_EXAMPLES)
    add_subdirectory(examples/ehdg_swe_manufactured_solution)
  endif()
endif()

if(IHDG_SWE)
  find_package(PETSc REQUIRED)

  add_subdirectory(source/problem/SWE/discretization_IHDG)
  if(BUILD_EXAMPLES)
    add_subdirectory(examples/ihdg_swe_manufactured_solution)
  endif()
endif()

if(EHDG_GN)
  find_package(PETSc REQUIRED)

  add_subdirectory(source/problem/Green-Naghdi/discretization_EHDG)
  if(BUILD_EXAMPLES)
    add_subdirectory(examples/ehdg_gn_analytical_solution)
  endif()
endif()

if(IHDG_GN)
  #[[add_subdirectory(source/problem/Green-Naghdi/discretization_IHDG)
  if(BUILD_EXAMPLES)
    add_subdirectory(examples/ihdg_gn_analytical_solution)
  endif()]]
endif()