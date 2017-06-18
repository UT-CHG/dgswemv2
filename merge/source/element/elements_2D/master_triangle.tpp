template<int dimension, int element_type, class basis_type, class integration_int_type, class integration_bound_type>
std::vector<Point<dimension>> MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>
::TriangleBoundaryToMasterCoordinates(int boundary, const std::vector<Point<dimension - 1>>& z_boundary) {
	std::vector<Point<dimension>> z_master(z_boundary.size());

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

template<int dimension, int element_type, class basis_type, class integration_int_type, class integration_bound_type>
std::vector<Point<dimension>> MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>::TriangleVTKPostCell() {
	std::vector<Point<dimension>> z_postprocessor_cell(N_DIV*N_DIV);

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

template<int dimension, int element_type, class basis_type, class integration_int_type, class integration_bound_type>
std::vector<Point<dimension>> MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>::TriangleVTKPostPoint() {
	std::vector<Point<dimension>> z_postprocessor_point((N_DIV + 1)*(N_DIV + 2) / 2);
		
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
