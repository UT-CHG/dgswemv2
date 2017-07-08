#ifndef GENERAL_DEFINITIONS
#define GENERAL_DEFINITIONS

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
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
#include <unordered_map>
#include <tuple>

typedef unsigned int uint;

template<uint dim>
using Point = std::array<double, dim>;

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
		virtual Array2D<double> GetPhi(uint, const std::vector<Point<dim>>&) = 0;
		virtual Array3D<double> GetDPhi(uint, const std::vector<Point<dim>>&) = 0;
		virtual std::pair<bool, Array2D<double>> GetMinv(uint) = 0;
	};
}

namespace Integration {
	template<uint dim>
	class Integration {
	public:
		virtual std::pair<std::vector<double>, std::vector<Point<dim>>> GetRule(uint) = 0;
	};
}

namespace Master {
	template<uint dim>
	class Master {
	public:
		uint p;

		Array2D<double> phi_gp;
		Array3D<double> dphi_gp;

		Array2D<double> int_fact_phi;
		Array3D<double> int_fact_dphi;

		std::pair<bool, Array2D<double>> m_inv;

		Array2D<double> phi_postprocessor_cell;
		Array2D<double> phi_postprocessor_point;

	public:
		Master(uint p) : p(p) {}

		virtual std::vector<Point<dim>> BoundaryToMasterCoordinates(uint, const std::vector<Point<dim - 1>>&) = 0;
	};
}

namespace Shape {
	template<uint dim>
	class Shape {
	public:
		std::vector<Point<dim>> nodal_coordinates;

	public:
		Shape(std::vector<Point<dim>>& nodal_coordinates) : nodal_coordinates(nodal_coordinates) {}

		virtual std::vector<double> GetJdet(const std::vector<Point<dim>>&) = 0;
		virtual Array3D<double> GetJinv(const std::vector<Point<dim>>&) = 0;
		virtual std::vector<double> GetSurfaceJ(uint, const std::vector<Point<dim>>&) = 0;
		virtual Array2D<double> GetSurfaceNormal(uint, const std::vector<Point<dim>>&) = 0;
		virtual void GetVTK(std::vector<Point<3>>&, Array2D<uint>&) = 0;
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
#define INTERNAL 255

// element types (as VTK cell types)
#define TRIANGLE 5

// postprocessor element divisions
#define N_DIV 2

#endif