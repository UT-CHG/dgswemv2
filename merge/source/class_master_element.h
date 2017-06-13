#ifndef CLASS_MASTER_ELEMENT_H
#define CLASS_MASTER_ELEMENT_H

#include "general_definitions.h"

#include "basis_functions/bases_2D.h"

#include "integration_rules/integration_rules_1D.h"
#include "integration_rules/integration_rules_2D.h"

//template<int element_type = TRIANGLE, class Basis = Dubiner_2D, 
//	class Integration_int = Dunavant_2D, class Integration_bound = GaussLegendre_1D>
class MasterElement {
public:
	Dubiner_2D basis;
	
	GaussLegendre_1D integration_boundary;
	Dunavant_2D integration_internal;

	std::pair<bool, Array2D<double>> m_inv;

	Array2D<double> phi_internal;
	Array3D<double> dphi_internal;
	Array3D<double> phi_boundary;

	Array2D<double> internal_int_fac_phi;
	Array3D<double> internal_int_fac_dphi;
	Array3D<double> boundary_int_fac_phi;

	Array2D<double> phi_postprocessor_cell;
	Array2D<double> phi_postprocessor_point;
public:
	MasterElement(int p) {
		int element_type = TRIANGLE; //for testing

		switch (element_type) {
		case TRIANGLE: this->MasterTriangle(p); break;
		default:
			printf("\n");
			printf("MASTER ELEMENT CONSTRUCTOR - Fatal error!\n");
			printf("Undefined element type = %d\n", element_type);
			exit(1);
		}
	}
private:
	void MasterTriangle(int);
};

#endif