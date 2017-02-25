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
  source/class_basis.h \
  source/class_basis_geometry.h \
  source/class_integration.h \
  source/class_element.h \
  source/basis_functions/basis_functions.h \
  source/integration_rules/integration_rules_area.h \
  source/integration_rules/integration_rules_line.h

srcs = \
  source/class_basis.cpp \
  source/class_basis_geometry.cpp \
  source/class_element.cpp \
  source/basis_functions/basis_triangle_dubiner.cpp \
  source/basis_functions/basis_polynomials.cpp \
  source/class_integration.cpp \
  source/integration_rules/integration_area_dunavant.cpp \
  source/integration_rules/integration_line_gausslegendre.cpp \
  source/mesh_processor/read_adcirc_mesh.cpp

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
	$(cxx) -std=c++11  -Wall $(flags) -o $@ $(srcs) source/$*.cpp
