#ifndef GENERAL_DEFINITIONS
#define GENERAL_DEFINITIONS

#include <iostream>
#include <fstream>
#include <string>

#include <cmath>

#include <array>
#include <vector>
#include <map>

#define PROBLEM_SWE_2D_LINEAR

#ifdef PROBLEM_SWE_2D
#include "problem/SWE_2D/definitions_SWE_2D.h"
#endif
#ifdef PROBLEM_SWE_2D_LINEAR
#include "problem/SWE_2D_LINEAR/definitions_SWE_2D_LINEAR.h"
#endif

template<int DIM>
using Point = std::array<double, DIM>;

template<class type>
using Array2D = std::vector<std::vector<type>>;

template<class type>
using Array3D = std::vector<std::vector<std::vector<type>>>;

template<class type>
using Array4D = std::vector<std::vector<std::vector<std::vector<type>>>>;

#define PI 3.14159265359	

#define DEFAULT_ID 4294967295
#define DEFAULT_BOUNDARY 255

#define X 0
#define Y 1
#define Z 2

#define Z1 0
#define Z2 1
#define Z3 2

#define N1 0
#define N2 1
#define N3 2

// timestepper types
#define EEULER 0
#define RK2 1
#define RK3 2
#define RK4 3

// element types (as VTK cell types)
#define TRIANGLE 5

// postprocessor element divisions
#define N_DIV 1

// 2D bases
#define DUBINER_2D 16
// 2D geometric Bases

// 2D integration rules
#define DUNAVANT_2D 16
// 1D integration rules
#define GAUSS_LEGENDRE_1D 0

#endif