#ifndef CLASS_MASTER_ELEMENT_H
#define CLASS_MASTER_ELEMENT_H

#include "general_definitions.h"

#include "basis/bases_2D.h"

#include "integration/integrations_1D.h"
#include "integration/integrations_2D.h"

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type>
class MasterElement {
public:
	unsigned char number_boundaries;	

	basis_type basis;
	int number_bf;

	integration_bound_type integration_boundary;
	int number_gp_boundary;
	
	integration_int_type integration_internal;
	int number_gp_internal;

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

#include "class_master_element.tpp"

#endif