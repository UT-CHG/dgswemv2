cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project (DG_HYPER CXX C)

option(CMAKE_VERBOSE_MAKEFILE "Verbose CMake" ON)
option(USE_OMPI "Use MPI+OpenMP" OFF)
option(USE_HPX "Use HPX" OFF)
enable_testing()

find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)

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

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

find_package(METIS REQUIRED)
find_package(yaml-cpp REQUIRED)

if(USE_OMPI)
  find_package(MPI REQUIRED)
  find_package(OpenMP REQUIRED)

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(source)
add_subdirectory(test)
add_subdirectory(mesh_generators)
add_subdirectory(partitioner)