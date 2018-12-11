cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project (DGSWEMv2 CXX C)

option(USE_EIGEN "Use Eigen" ON)
option(USE_BLAZE "Use Blaze" OFF)
option(USE_OMPI "Use MPI+OpenMP" OFF)
option(USE_HPX "Use HPX" OFF)
option(COMPILER_WARNINGS "Enable Compiler Warnings" OFF)
option(SET_VERBOSE "Verbose MakeFile" OFF)

option(RKDG "Build with RKDG discretization" ON)
option(EHDG "Build with explicit HDG discretization" OFF)
option(IHDG "Build with implicit HDG discretization" OFF)

option(SWE "Build with shallow water equations support" ON)
option(GN "Build with Green-Nahgdi equations support" OFF)

enable_testing()

find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # require at least intel 17.0
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17.0)
    message(FATAL_ERROR "Intel version must be at least 17.0")
  endif()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # require at least 6.0
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
    message(FATAL_ERROR "GCC version must be at least 6.0")
  endif()
else()
  message(WARNING "Using untested compiler")
endif()

add_definitions(-std=c++14 -g)

if(COMPILER_WARNINGS)
  add_definitions(-Wall)
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

if(NOT USE_EIGEN AND NOT USE_BLAZE) 
  message(WARNING "Linear algebra package not specified! Using default: Eigen!")
  set(USE_EIGEN ON)
endif()

if(USE_EIGEN)
  find_package(Eigen3 REQUIRED)
  list(APPEND LINALG_DEFINITION USE_EIGEN)
  link_libraries(Eigen3::Eigen)
endif()

if(USE_BLAZE)
  find_package(blaze REQUIRED)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  list(APPEND LINALG_DEFINITION USE_BLAZE)
  link_libraries(blaze::blaze ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

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

  ${PROJECT_SOURCE_DIR}/source/simulation/stepper/explicit_ssp_rk_stepper.cpp

  ${PROJECT_SOURCE_DIR}/source/basis/polynomials/basis_polynomials.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_1D/basis_legendre_1D.cpp
  ${PROJECT_SOURCE_DIR}/source/basis/bases_2D/basis_dubiner_2D.cpp

  ${PROJECT_SOURCE_DIR}/source/shape/shapes_2D/shape_straighttriangle.cpp

  ${PROJECT_SOURCE_DIR}/source/integration/integrations_2D/integration_dunavant_2D.cpp
  ${PROJECT_SOURCE_DIR}/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
)

#If solving the Green-Naghdi equations, SWE_SUPPORT is implicitly implied
if (GN)
  set(SWE ON)
  set(EHDG ON)
  set(IHDG ON) #implicitly requires PETSc for OMPI targets
endif()

if (SWE)
  list(APPEND SOURCES
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_input/swe_inputs.cpp
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_function_files/swe_initial_condition_functions.cpp
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_function_files/swe_source_functions.cpp
    ${PROJECT_SOURCE_DIR}/source/problem/SWE/problem_function_files/swe_true_solution_functions.cpp
  )
  list( APPEND PROBLEM_DEFINITIONS SWE_SUPPORT )
endif()

if (GN)
  list(APPEND PROBLEM_DEFINITIONS GN_SUPPORT)
endif()

if (RKDG)
  list(APPEND PROBLEM_DEFINITIONS RKDG_SUPPORT)
endif()
if (EHDG)
  list(APPEND PROBLEM_DEFINITIONS EHDG_SUPPORT)
endif()
if (IHDG)
  if(USE_OMPI)
    find_package(PETSc REQUIRED)
  endif()
  list(APPEND PROBLEM_DEFINITIONS IHDG_SUPPORT)
endif()

include_directories(source)

add_subdirectory(test)
add_subdirectory(mesh_generators)
add_subdirectory(partitioner)

add_subdirectory(source)