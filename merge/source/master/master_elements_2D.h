#ifndef CLASS_MASTER_ELEMENT_H
#define CLASS_MASTER_ELEMENT_H

#include "../general_definitions.h"

#include "../basis/bases_2D.h"

#include "../integration/integrations_1D.h"
#include "../integration/integrations_2D.h"

namespace Master {
	template<class basis_type, class integration_type>
	class Triangle : Master<2> {
	public:
		int p;

		basis_type basis;

		Array2D<double> phi_gp;
		Array3D<double> dphi_gp;

		Array2D<double> int_fact_phi;
		Array3D<double> int_fact_dphi;

		std::pair<bool, Array2D<double>> m_inv;

		Array2D<double> phi_postprocessor_cell;
		Array2D<double> phi_postprocessor_point;

	public:
		Triangle(int);
		std::vector<Point<2>> BoundaryToMasterCoordinates(int, const std::vector<Point<1>>&);

	private:
		std::vector<Point<2>> VTKPostCell();
		std::vector<Point<2>> VTKPostPoint();
	};
}

#include "elements_2D/master_triangle.tpp"

#endif