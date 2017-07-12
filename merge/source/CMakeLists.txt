cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(DG_HYPER_SWE CXX)

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
  messsage(WARNING "Using untested compiler")
endif()

set(CMAKE_VERBOSE_MAKEFILE off)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++14 -Wall")

list (APPEND DG_HYPER_SWE_SOURCES
  /h2/kazbek/UT-CHG/dgswemv2/merge/source/mesh_metadata.cpp
  /h2/kazbek/UT-CHG/dgswemv2/merge/source/ADCIRC_reader/adcirc_format.cpp
  
  /h2/kazbek/UT-CHG/dgswemv2/merge/source/stepper.cpp

  /h2/kazbek/UT-CHG/dgswemv2/merge/source/basis/polynomials/basis_polynomials.cpp
  /h2/kazbek/UT-CHG/dgswemv2/merge/source/basis/bases_2D/basis_dubiner_2D.cpp
   
  /h2/kazbek/UT-CHG/dgswemv2/merge/source/shape/shapes_2D/shape_straighttriangle.cpp

  /h2/kazbek/UT-CHG/dgswemv2/merge/source/integration/integrations_2D/integration_dunavant_2D.cpp
  /h2/kazbek/UT-CHG/dgswemv2/merge/source/integration/integrations_1D/integration_gausslegendre_1D.cpp
)
      
add_executable(
  DG_HYPER_SWE
  /h2/kazbek/UT-CHG/dgswemv2/merge/source/problem/SWE/main_SWE.cpp
  ${DG_HYPER_SWE_SOURCES}
)