cmake_minimum_required (VERSION 2.8 FATAL_ERROR)
project (DG_HYPER CXX)
 
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

set(CMAKE_VERBOSE_MAKEFILE on)

# default installation
get_filename_component (default_prefix "../install" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in install."
FORCE)

include_directories(source)

enable_testing()

add_subdirectory(source)
add_subdirectory(test)
