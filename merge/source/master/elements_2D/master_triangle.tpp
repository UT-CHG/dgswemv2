#include "../master_elements_2D.h"

namespace Master {
	template<class basis_type, class integration_type>
	Triangle<basis_type, integration_type>::Triangle(int p) : p(p) {
		integration_type integration;
		std::pair<std::vector<double>, std::vector<Point<2>>> integration_rule = integration.get_rule(2 * p);

		this->phi_gp = this->basis.get_phi(p, integration_rule.second);
		this->dphi_gp = this->basis.get_dphi(p, integration_rule.second);

		std::vector<Point<2>> z_postprocessor_cell = this->VTKPostCell();
		this->phi_postprocessor_cell = this->basis.get_phi(p, z_postprocessor_cell);

		std::vector<Point<2>> z_postprocessor_point = this->VTKPostPoint();
		this->phi_postprocessor_point = this->basis.get_phi(p, z_postprocessor_point);

		this->int_fact_phi = this->phi_gp;
		for (int i = 0; i < this->int_fact_phi.size(); i++) { //iterate through basis functions
			for (int j = 0; j < int_fact_phi[i].size(); j++) { //iterate through internal GPs
				this->int_fact_phi[i][j] *= integration_rule.first[j]; //apply weight
			}
		}

		this->int_fact_dphi = this->dphi_gp;
		for (int i = 0; i < int_fact_dphi.size(); i++) { //iterate through basis functions
			for (int j = 0; j < int_fact_dphi[i].size(); j++) { //iterate through differentiation directions
				for (int k = 0; k < int_fact_dphi[i][j].size(); k++) { //iterate through internal GPs
					this->int_fact_dphi[i][j][k] *= integration_rule.first[k]; //apply weight
				}
			}
		}

		this->m_inv = this->basis.get_m_inv(p);
	}

	template<class basis_type, class integration_type>
	std::vector<Point<2>> Triangle<basis_type, integration_type>::boundary_to_master(int boundary, const std::vector<Point<1>>& z_boundary) {
		std::vector<Point<2>> z_master(z_boundary.size());

		if (boundary == 0) {
			for (int i = 0; i < z_master.size(); i++) {
				z_master[i][Z1] = -z_boundary[i][Z1];
				z_master[i][Z2] = z_boundary[i][Z1];
			}
		}
		else if (boundary == 1) {
			for (int i = 0; i < z_master.size(); i++) {
				z_master[i][Z1] = -1;
				z_master[i][Z2] = -z_boundary[i][Z1];
			}
		}
		else if (boundary == 2) {
			for (int i = 0; i < z_master.size(); i++) {
				z_master[i][Z1] = z_boundary[i][Z1];
				z_master[i][Z2] = -1;
			}
		}

		return z_master;
	}

	template<class basis_type, class integration_type>
	std::vector<Point<2>> Triangle<basis_type, integration_type>::VTKPostCell() {
		std::vector<Point<2>> z_postprocessor_cell(N_DIV*N_DIV);

		double dz = 2.0 / N_DIV;

		int n_pt = 0;
		for (int i = 0; i < N_DIV; i++) {
			for (int j = 0; j < N_DIV - i; j++) {
				z_postprocessor_cell[n_pt][Z1] = -1.0 + dz*j + dz / 3.0; //CENTROID
				z_postprocessor_cell[n_pt][Z2] = -1.0 + dz*i + dz / 3.0;
				n_pt++;
			}
		}
		for (int i = 1; i < N_DIV; i++) {
			for (int j = 0; j < N_DIV - i; j++) {
				z_postprocessor_cell[n_pt][Z1] = -1.0 + dz*j + 2 * dz / 3.0; //CENTROID
				z_postprocessor_cell[n_pt][Z2] = -1.0 + dz*i - dz / 3.0;
				n_pt++;
			}
		}

		return z_postprocessor_cell;
	}

	template<class basis_type, class integration_type>
	std::vector<Point<2>> Triangle<basis_type, integration_type>::VTKPostPoint() {
		std::vector<Point<2>> z_postprocessor_point((N_DIV + 1)*(N_DIV + 2) / 2);

		double dz = 2.0 / N_DIV;

		int n_pt = 0;
		for (int i = 0; i < N_DIV; i++) {
			for (int j = 0; j <= N_DIV - i; j++) {
				z_postprocessor_point[n_pt][Z1] = -1.0 + dz*j;
				z_postprocessor_point[n_pt][Z2] = -1.0 + dz*i;
				n_pt++;
			}
		}

		return z_postprocessor_point;
	}
}