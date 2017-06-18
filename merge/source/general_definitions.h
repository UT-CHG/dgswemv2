#ifndef GENERAL_DEFINITIONS
#define GENERAL_DEFINITIONS

#include <iostream>
#include <fstream>
#include <string>

#include <cmath>
#include <functional>

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

namespace Basis {
	template<int dim>
	class Basis {
	public:
		virtual Array2D<double> get_phi(int, const std::vector<Point<dim>>&) = 0;
		virtual Array3D<double> get_dphi(int, const std::vector<Point<dim>>&) = 0;
		virtual std::pair<bool, Array2D<double>> get_m_inv(int) = 0;
	};
}

namespace Integration {
	template<int dim>
	class Integration {
	public:
		virtual std::pair<std::vector<double>, std::vector<Point<dim>>> get_rule(int) = 0;
	};
}

namespace Shape {
	template<int dim>
	class Shape {
	public:
		virtual std::vector<double> get_J_det(const std::vector<Point<dim>>&) = 0;
		virtual Array3D<double> get_J_inv(const std::vector<Point<dim>>&) = 0;
		virtual Array2D<double> get_surface_J(const std::vector<Point<dim>>&) = 0;
		virtual Array3D<double> get_surface_normal(const std::vector<Point<dim>>&) = 0;

		virtual std::vector<double> get_surface_J_(int, const std::vector<Point<dim>>&) = 0;
		virtual Array2D<double> get_surface_normal_(int, const std::vector<Point<dim>>&) = 0;

		virtual void get_VTK(std::vector<Point<3>>&, Array2D<unsigned int>&, const std::vector<Point<dim>>&) = 0;
	};
}

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

//TIMESTEPPING FLAGS FOR MULTISTEP METHODS
#define U0 0U
#define K1 1U
#define K2 2U
#define K3 3U
#define K4 4U

#endif