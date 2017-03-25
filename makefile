SHELL = /bin/bash

# if a better c++ compiler can't be found
cxx_default = g++

# look for CC, else use cxx_default
cxx = $(shell (which CC > /dev/null && echo CC) || echo $(cxx_default))

#dirs = \
  source/ \
  source/basis_functions \
  source/integration_rules

hdrs = \
  source/general_definitions.h \
  source/class_problem.h \
  source/problem/problem_SWE_2D.h \
  source/class_mesh_v2.h \
  source/class_integration.h \
  source/class_basis.h \
  source/integration_rules/integration_rules_area.h \
  source/integration_rules/integration_rules_line.h \
  source/basis_functions/basis_functions.h \
  source/class_basis_geometry.h \
  source/class_element.h \
  source/class_interface.h \
  source/elements/class_element_2D.h \
  source/interfaces/class_interface_2D.h

srcs = \
  source/problem/problem_SWE_2D.cpp \
  source/class_mesh_v2.cpp \
  source/class_integration.cpp \
  source/integration_rules/integration_area_dunavant.cpp \
  source/integration_rules/integration_line_gausslegendre.cpp \
  source/class_basis.cpp \
  source/basis_functions/basis_polynomials.cpp \
  source/basis_functions/basis_triangle_dubiner.cpp \
  source/class_basis_geometry.cpp \
  source/elements/class_element_2D.cpp \
  source/elements/elements_2D/element_tri.cpp \
  source/interfaces/class_interface_2D.cpp 
#  source/mesh_processor/read_adcirc_mesh.cpp

exe/%.o2:  flags += -O2 -DDEBUG=1
exe/%.o3:  flags += -O3 -DDEBUG=0
exe/%.gdb: flags += -O0 -g -DDEBUG=1 #-D_GLIBCXX_DEBUG
exe/%.vg:  flags += -O1 -g -DDEBUG=1
exe/%.gp:  flags += -O3 -DDEBUG=0 -pg

exe/%.o2 \
exe/%.o3 \
exe/%.gdb \
exe/%.vg \
exe/%.gp: $(srcs) $(hdrs) source/%.cpp
	mkdir -p exe
	$(cxx) -std=c++1y  -Wall $(flags) -o $@ $(srcs) source/$*.cpp
