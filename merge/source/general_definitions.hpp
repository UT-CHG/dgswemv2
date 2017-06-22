#ifndef GENERAL_DEFINITIONS
#define GENERAL_DEFINITIONS

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <functional>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <assert.h>
#include <array>
#include <vector>
#include <map>

typedef unsigned int uint;

template<uint DIM>
using Point = std::array<double, DIM>;

template<class type>
using Array2D = std::vector<std::vector<type>>;

template<class type>
using Array3D = std::vector<std::vector<std::vector<type>>>;

template<class type>
using Array4D = std::vector<std::vector<std::vector<std::vector<type>>>>;

namespace Basis {
	template<uint dim>
	class Basis {
	public:
		virtual Array2D<double> get_phi(uint, const std::vector<Point<dim>>&) = 0;
		virtual Array3D<double> get_dphi(uint, const std::vector<Point<dim>>&) = 0;
		virtual std::pair<bool, Array2D<double>> get_m_inv(uint) = 0;
	};
}

namespace Integration {
	template<uint dim>
	class Integration {
	public:
		virtual std::pair<std::vector<double>, std::vector<Point<dim>>> get_rule(uint) = 0;
	};
}

namespace Master {
	template<uint dim>
	class Master {
	public:
		virtual std::vector<Point<dim>> boundary_to_master(uint, const std::vector<Point<dim - 1>>&) = 0;
	};
}

namespace Shape {
	template<uint dim>
	class Shape {
	public:
		virtual std::vector<double> get_J_det(const std::vector<Point<dim>>&) = 0;
		virtual Array3D<double> get_J_inv(const std::vector<Point<dim>>&) = 0;
		virtual std::vector<double> get_surface_J(uint, const std::vector<Point<dim>>&) = 0;
		virtual Array2D<double> get_surface_normal(uint, const std::vector<Point<dim>>&) = 0;
		virtual void get_VTK(std::vector<Point<3>>&, Array2D<uint>&, const std::vector<Point<dim>>&) = 0;
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

#define OCEAN 0
#define LAND 1
#define INTERNAL 2

// element types (as VTK cell types)
#define TRIANGLE 5

// postprocessor element divisions
#define N_DIV 1

#endif