#ifndef CLASS_MASTER_ELEMENT_H
#define CLASS_MASTER_ELEMENT_H

#include "../general_definitions.h"

#include "../basis/bases_2D.h"

#include "../integration/integrations_1D.h"
#include "../integration/integrations_2D.h"

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type>
class MasterElement {
public:
	int p;
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
	MasterElement(int);
	std::vector<Point<dimension>> TriangleBoundaryToMasterCoordinates(int, const std::vector<Point<dimension - 1>>&);

private:
	std::vector<Point<dimension>> TriangleVTKPostCell();
	std::vector<Point<dimension>> TriangleVTKPostPoint();
};

template<int dimension, int element_type, class basis_type,
	class integration_int_type, class integration_bound_type>
MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>
::MasterElement(int p) : p(p) {
	std::pair<std::vector<double>, std::vector<Point<dimension - 1>>> integration_rule_boundary = this->integration_boundary.get_rule(2 * p);
	std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule_internal = this->integration_internal.get_rule(2 * p);

	this->number_gp_boundary = integration_rule_boundary.first.size();
	this->number_gp_internal = integration_rule_internal.first.size();

	this->m_inv = this->basis.get_m_inv(p);

	this->number_bf = this->m_inv.second[0].size();

	this->phi_internal = this->basis.get_phi(p, integration_rule_internal.second);
	this->basis.basis_test(p, this->phi_internal, integration_rule_internal);

	this->dphi_internal = this->basis.get_dphi(p, integration_rule_internal.second);

	//this block requires template specialization or switch statement
	{
		this->number_boundaries = 3;
		this->phi_boundary.reserve(this->number_boundaries);
		std::vector<Point<dimension>> z_boundary;

		for (int i = 0; i < this->number_boundaries; i++) {
			z_boundary = this->TriangleBoundaryToMasterCoordinates(i, integration_rule_boundary.second);
			this->phi_boundary.push_back(basis.get_phi(p, z_boundary));
		}

		std::vector<Point<dimension>> z_postprocessor_cell = this->TriangleVTKPostCell();
		this->phi_postprocessor_cell = this->basis.get_phi(p, z_postprocessor_cell);

		std::vector<Point<dimension>> z_postprocessor_point = this->TriangleVTKPostPoint();
		this->phi_postprocessor_point = this->basis.get_phi(p, z_postprocessor_point);
	}
	//end block

	this->internal_int_fac_phi = this->phi_internal;
	for (int i = 0; i < this->internal_int_fac_phi.size(); i++) { //iterate through basis functions
		for (int j = 0; j < internal_int_fac_phi[i].size(); j++) { //iterate through internal GPs
			this->internal_int_fac_phi[i][j] *= integration_rule_internal.first[j]; //apply weight
		}
	}

	this->internal_int_fac_dphi = this->dphi_internal;
	for (int i = 0; i < internal_int_fac_dphi.size(); i++) { //iterate through basis functions
		for (int j = 0; j < internal_int_fac_dphi[i].size(); j++) { //iterate through differentiation directions
			for (int k = 0; k < internal_int_fac_dphi[i][j].size(); k++) { //iterate through internal GPs
				this->internal_int_fac_dphi[i][j][k] *= integration_rule_internal.first[k]; //apply weight
			}
		}
	}

	this->boundary_int_fac_phi = this->phi_boundary;
	for (int i = 0; i < boundary_int_fac_phi.size(); i++) { //iterate thorough boundaries
		for (int j = 0; j < boundary_int_fac_phi[i].size(); j++) { //iterate through basis funtions
			for (int k = 0; k < boundary_int_fac_phi[i][j].size(); k++) { //iterate through boundary GPs
				this->boundary_int_fac_phi[i][j][k] *= integration_rule_boundary.first[k]; //apply weight
			}
		}
	}
}


#include "elements_2D/master_triangle.tpp"

#endif