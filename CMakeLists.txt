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

get_filename_component (default_prefix "../install" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in install."
FORCE)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

find_package(METIS REQUIRED)
find_package(yaml-cpp REQUIRED)

if(USE_OMPI)
  find_package(OpenMP REQUIRED)
  find_package(MPI REQUIRED)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(USE_HPX)
  find_package(HPX REQUIRED)
endif()

if(COMPILER_WARNINGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

set(CMAKE_VERBOSE_MAKEFILE OFF)
if(SET_VERBOSE)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

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
  find_package(Eigen3 REQUIRED)

  find_package( blaze )
  if( blaze_FOUND )
     add_library( blaze_target INTERFACE )
     target_link_libraries( blaze_target INTERFACE blaze::blaze )
  endif()

  add_subdirectory(source/problem/SWE/discretization_EHDG)
  if(BUILD_EXAMPLES)
    add_subdirectory(examples/ehdg_swe_manufactured_solution)
  endif()
endif()

if(IHDG_SWE)
  find_package(Eigen3 REQUIRED)

  add_subdirectory(source/problem/SWE/discretization_IHDG)
  #if(BUILD_EXAMPLES)
  #  add_subdirectory(examples/ehdg_swe_manufactured_solution)
  #endif()
endif()
