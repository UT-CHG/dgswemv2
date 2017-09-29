cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project (DG_HYPER CXX C)

option(USE_OMPI "Use MPI+OpenMP" OFF)
option(USE_HPX "Use HPX" OFF)
option(USE_CUDA "Use CUDA" ON)
option(SET_VERBOSE "Verbose MakeFile" OFF)

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

if(USE_CUDA)
  find_package(CUDA REQUIRED)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_60,code=sm_60 -std=c++11)
endif()

if(USE_OMPI)
  find_package(OpenMP REQUIRED)
  find_package(MPI REQUIRED)

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(USE_HPX)
  find_package(HPX REQUIRED)
endif()

set(CMAKE_VERBOSE_MAKEFILE OFF)
if(SET_VERBOSE)
  set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

include_directories(source)

add_subdirectory(source)
add_subdirectory(test)
add_subdirectory(mesh_generators)
add_subdirectory(partitioner)
