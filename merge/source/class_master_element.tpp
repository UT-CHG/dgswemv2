template<int dimension, int element_type, class basis_type, class integration_int_type, class integration_bound_type>
void MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>::MasterTriangle(int p) {
	this->number_boundaries = 3;

	std::pair<std::vector<double>, std::vector<Point<1>>> integration_rule_boundary = this->integration_boundary.get_rule(2 * p);
	std::pair<std::vector<double>, std::vector<Point<2>>> integration_rule_internal = this->integration_internal.get_rule(2 * p);

	this->number_gp_boundary = integration_rule_boundary.first.size();
	this->number_gp_internal = integration_rule_internal.first.size();

	this->m_inv = this->basis.get_m_inv(p);

	this->number_bf = this->m_inv.second[0].size();

	this->phi_internal = this->basis.get_phi(p, integration_rule_internal.second);
	this->basis.basis_test(p, this->phi_internal, integration_rule_internal); //TEST BASIS

	this->dphi_internal = this->basis.get_dphi(p, integration_rule_internal.second);

	this->phi_boundary.reserve(3);
	std::vector<Point<2>> z_boundary(integration_rule_boundary.first.size());
	for (int i = 0; i < 3; i++) {
		if (i == 0) {
			for (size_t j = 0; j < z_boundary.size(); j++) {
				z_boundary[j][Z1] = -integration_rule_boundary.second[j][Z1];
				z_boundary[j][Z2] = integration_rule_boundary.second[j][Z1];
			}
		}
		else if (i == 1) {
			for (size_t j = 0; j < z_boundary.size(); j++) {
				z_boundary[j][Z1] = -1;
				z_boundary[j][Z2] = -integration_rule_boundary.second[j][Z1];
			}
		}
		else if (i == 2) {
			for (size_t j = 0; j < z_boundary.size(); j++) {
				z_boundary[j][Z1] = integration_rule_boundary.second[j][Z1];
				z_boundary[j][Z2] = -1;
			}
		}
		this->phi_boundary.push_back(basis.get_phi(p, z_boundary));
	}

	this->internal_int_fac_phi = this->phi_internal;
	for (size_t i = 0; i < this->internal_int_fac_phi.size(); i++) { //iterate through basis functions
		for (size_t j = 0; j < internal_int_fac_phi[i].size(); j++) { //iterate through internal GPs
			this->internal_int_fac_phi[i][j] *= integration_rule_internal.first[j]; //apply weight
		}
	}

	this->internal_int_fac_dphi = this->dphi_internal;
	for (size_t i = 0; i < internal_int_fac_dphi.size(); i++) { //iterate through basis functions
		for (size_t j = 0; j < internal_int_fac_dphi[i].size(); j++) { //iterate through differentiation directions
			for (size_t k = 0; k < internal_int_fac_dphi[i][j].size(); k++) { //iterate through internal GPs
				this->internal_int_fac_dphi[i][j][k] *= integration_rule_internal.first[k]; //apply weight
			}
		}
	}

	this->boundary_int_fac_phi = this->phi_boundary;
	for (size_t i = 0; i < boundary_int_fac_phi.size(); i++) { //iterate thorough boundaries
		for (size_t j = 0; j < boundary_int_fac_phi[i].size(); j++) { //iterate through basis funtions
			for (size_t k = 0; k < boundary_int_fac_phi[i][j].size(); k++) { //iterate through boundary GPs
				this->boundary_int_fac_phi[i][j][k] *= integration_rule_boundary.first[k]; //apply weight
			}
		}
	}

	double dz = 2.0 / N_DIV;

	std::vector<Point<2>> z_postprocessor_cell(N_DIV*N_DIV);
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

	this->phi_postprocessor_cell = this->basis.get_phi(p, z_postprocessor_cell);

	std::vector<Point<2>> z_postprocessor_point((N_DIV + 1)*(N_DIV + 2) / 2);
	n_pt = 0;
	for (int i = 0; i < N_DIV; i++) {
		for (int j = 0; j <= N_DIV - i; j++) {
			z_postprocessor_point[n_pt][Z1] = -1.0 + dz*j;
			z_postprocessor_point[n_pt][Z2] = -1.0 + dz*i;
			n_pt++;
		}
	}

	this->phi_postprocessor_point = this->basis.get_phi(p, z_postprocessor_point);
}