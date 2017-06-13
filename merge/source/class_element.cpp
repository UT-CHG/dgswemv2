#include "class_element.h"

ELEMENT::ELEMENT(MasterElement& master,
	unsigned int ID, unsigned int* neighbor_ID, unsigned char* boundary_type, 
	double** nodal_coordinates, BASIS_GEOM* basis_geom) :

	master(master), m_inv(master.m_inv),
	phi_internal(master.phi_internal), phi_boundary(master.phi_boundary),
	phi_postprocessor_cell(master.phi_postprocessor_cell),
	phi_postprocessor_point(master.phi_postprocessor_point),

	ID(ID), nodal_coordinates(nodal_coordinates),
	neighbor_ID(neighbor_ID), boundary_type(boundary_type)
{
    if (basis_geom != nullptr) {
        this->basis_geom = basis_geom;
    }

	this->element_type = TRIANGLE; //TEST
	switch (element_type){
	case TRIANGLE: this->Triangle(); break;
	default:
		printf("\n");
		printf("ELEMENT CONSTRUCTOR - Fatal error!\n");
		printf("Undefined element type = %d\n", element_type);
		exit(1);
	}

    //COMPUTE DIFFERENTIATION AND NUMERICAL INTEGRATION FACTORS
	this->ComputeDifferentiationFactors();
    this->ComputeIntegrationFactors();

	//SET INITIAL CONDITIONS FOR TESTING
	this->u[SP][0] = 1; //NO SPHERICAL CORRECTION
	this->u[ZB][0] = 3; //FLAT BED
	this->u[2][0] = 0; // ZERO VELOCITY FIELD
	this->u[3][0] = 0;
	this->u[4][0] = 0; //FLAT SURFACE

	for (int i = 1; i < this->number_bf; i++) {
		this->u[SP][i] = 0;
		this->u[ZB][i] = 0; 
		this->u[2][i] = 0;
		this->u[3][i] = 0;
		this->u[4][i] = 0;
	}

	double L = 90000;
	double w = 2 * PI / 43200;
	double beta = w * sqrt(1 / (this->u[ZB][0] * GRAVITY));

	double h_true[3] = {
		0.3*cos(beta * this->nodal_coordinates[X][0]) / cos(beta * L),
		0.3*cos(beta * this->nodal_coordinates[X][1]) / cos(beta * L),
		0.3*cos(beta * this->nodal_coordinates[X][2]) / cos(beta * L),
	};

	this->u[4][0] = h_true[0] / 3.0 + h_true[1] / 3.0 + h_true[2] / 3.0;
	this->u[4][1] = -h_true[0] / 6.0 - h_true[1] / 6.0 + h_true[2] / 3.0;
	this->u[4][2] = -h_true[0] / 2.0 + h_true[1] / 2.0;

	//INITIALIZE SP AND ZB AT GPs
	this->ComputeBoundaryU(SP);
	this->ComputeInternalU(SP);
	this->ComputeBoundaryU(ZB);
	this->ComputeInternalU(ZB);
}

void ELEMENT::allocate_memory() {
	//INITIALIZE ARRAYS TO STORE INTERFACE DATA
	this->interfaces = new INTERFACE*[this->number_boundaries];
	this->interface_owner = new bool[this->number_boundaries];

	//INITIALIZE ARRAYS TO STORE Us and RHS
	this->RHS = new double[this->number_bf];

	this->u = new double*[SIZE_U];
	for (int i = 0; i < SIZE_U; i++) {
		this->u[i] = new double[this->number_bf];
	}

	this->u_internal = new double*[SIZE_U_INTERNAL];
	for (int i = 0; i < SIZE_U_INTERNAL; i++) {
		this->u_internal[i] = new double[this->number_gp_internal];
	}

	this->u_boundary = new double**[this->number_boundaries];
	for (int i = 0; i < this->number_boundaries; i++) {
		this->u_boundary[i] = new double*[SIZE_U_BOUNDARY];
		for (int j = 0; j < SIZE_U_BOUNDARY; j++) {
			this->u_boundary[i][j] = new double[this->number_gp_boundary];
		}
	}

	//INITIALIZE ARRAYS TO STORE INTEGRATION AND DIFFERENTIATION FACORS
	this->internal_dphi_fac.resize(this->number_bf);
	for (int i = 0; i < this->number_bf; i++) {
		this->internal_dphi_fac[i].resize(this->dimension);
		for (int j = 0; j < this->dimension; j++) {
			this->internal_dphi_fac[i][j].reserve(this->number_gp_internal);
		}
	}

	this->internal_int_fac_phi.resize(this->number_bf);
	for (int i = 0; i < this->number_bf; i++) {
		this->internal_int_fac_phi[i].reserve(this->number_gp_internal);
	}

	this->internal_int_fac_dphi.resize(this->number_bf);
	for (int i = 0; i < this->number_bf; i++) {
		this->internal_int_fac_dphi[i].resize(this->dimension);
		for (int j = 0; j < this->dimension; j++) {
			this->internal_int_fac_dphi[i][j].reserve(this->number_gp_internal);
		}
	}

	this->boundary_int_fac_phi.resize(this->number_boundaries);
	for (int i = 0; i < this->number_boundaries; i++) {
		this->boundary_int_fac_phi[i].resize(this->number_bf);
		for (int j = 0; j < this->number_bf; j++) {
			this->boundary_int_fac_phi[i][j].reserve(this->number_gp_boundary);
		}
	}

	//INITIALIZE ARRAYS TO STORE GEOMETRIC DATA
	if (this->basis_geom == nullptr) {
		this->det_J_internal.reserve(1);

		this->J_inv_internal.resize(this->dimension);
		for (int i = 0; i < this->dimension; i++) {
			this->J_inv_internal[i].resize(this->dimension);
			for (int j = 0; j < this->dimension; j++) {
				this->J_inv_internal[i][j].reserve(1);
			}
		}

		this->surface_J_boundary.resize(this->number_boundaries);
		for (int i = 0; i < this->number_boundaries; i++) {
			this->surface_J_boundary[i].reserve(1);
		}

		this->normal_boundary.resize(this->number_boundaries);
		for (int i = 0; i < this->number_boundaries; i++) {
			this->normal_boundary[i].resize(this->dimension);
			for (int j = 0; j < this->dimension; j++) {
				this->normal_boundary[i][j].reserve(1);
			}
		}
	}
	else {
		//Placeholder for cases p_geom > 1
	}
}

void ELEMENT::ComputeDifferentiationFactors() {
	//This loop can be rearranged for better efficiency,
	//however that will require two nested loops
	double d_phi;
	if (this->basis_geom == nullptr) {
		for (int i = 0; i < this->number_bf; i++) {
			for (int j = 0; j < this->dimension; j++) {
				for (int k = 0; k < this->number_gp_internal; k++) {
					d_phi = 0;
					for (int l = 0; l < this->dimension; l++) {
						d_phi += this->master.dphi_internal[i][l][k] * this->J_inv_internal[l][j][0];//FIX
					}
					this->internal_dphi_fac[i][j][k] = d_phi;
				}
			}
		}
	}
	else {
		//Placeholder for cases p_geom > 1
	}
}

void ELEMENT::ComputeIntegrationFactors() 
{
	if (this->basis_geom == nullptr) {
		this->internal_int_fac_phi = this->master.internal_int_fac_phi;

		double d_phi;
		for (int i = 0; i < this->number_bf; i++) {
			for (int j = 0; j < this->dimension; j++) {
				for (int k = 0; k < this->number_gp_internal; k++) {
					d_phi = 0;
					for (int l = 0; l < this->dimension; l++) {
						d_phi += this->master.internal_int_fac_dphi[i][l][k] * this->J_inv_internal[l][j][0];//FIX
					}
					this->internal_int_fac_dphi[i][j][k] = d_phi;
				}
			}
		}

		this->boundary_int_fac_phi = this->master.boundary_int_fac_phi;
		for (int i = 0; i < this->number_boundaries; i++) {
			for (int j = 0; j < this->number_bf; j++) {
				for (int k = 0; k < this->number_gp_boundary; k++) {
					this->boundary_int_fac_phi[i][j][k] *= this->surface_J_boundary[i][0] / abs(this->det_J_internal[0]);
				}
			}
		}
	}
	else {
		//Placeholder for cases p_geom > 1
	}
}

std::map<unsigned int, INTERFACE*> ELEMENT::CreateInterfaces() {
	std::map<unsigned int, INTERFACE*> internal_interfaces;

	bool boundary;
	for (int i = 0; i < this->number_boundaries; i++) {
		if (this->interfaces[i] == nullptr) {
			if (this->neighbor_ID[i] == DEFAULT_ID) boundary = true;
			else if (this->neighbor_ID[i] != DEFAULT_ID) boundary = false;
			
			double** normal = new double*[this->dimension];
			for (int j = 0; j < this->dimension; j++) {
				normal[j] = new double[this->number_gp_boundary];
			}

			if (this->basis_geom == nullptr) {
				for (int j = 0; j < this->dimension; j++) {
					for (int k = 0; k < this->number_gp_boundary; k++) {
						normal[j][k] = this->normal_boundary[i][j][0];
					}
				}
			}
			else {
				for (int j = 0; j < this->dimension; j++) {
					for (int k = 0; k < this->number_gp_boundary; k++) {
						normal[j][k] = this->normal_boundary[i][j][k];
					}
				}
			}

			this->interfaces[i] = new INTERFACE(this->dimension, this->number_gp_boundary,
				this->u_boundary[i], normal, boundary);

			this->interface_owner[i] = true;

			if (this->neighbor_ID[i] != DEFAULT_ID) {
				internal_interfaces[this->neighbor_ID[i]] = this->interfaces[i];
			}
		}
	}

	return internal_interfaces;
}

void ELEMENT::AppendInterface(unsigned int neighbor_ID, INTERFACE* interface_ptr) {
    for (int i = 0; i < this->number_boundaries; i++) {
        if (this->neighbor_ID[i] == neighbor_ID) {
            this->interfaces[i] = interface_ptr;

            this->interfaces[i]->SetPointerEX(this->u_boundary[i]);
        }
    }
}

std::vector<std::pair<unsigned char, INTERFACE*>> ELEMENT::GetOwnInterfaces() {
    std::vector<std::pair<unsigned char, INTERFACE*>> own_interfaces;

    for (int i = 0; i < this->number_boundaries; i++) {
        if (this->interface_owner[i]) {
            std::pair<unsigned char, INTERFACE*> own_interface;

            own_interface.first = this->boundary_type[i];
            own_interface.second = this->interfaces[i];
            own_interfaces.push_back(own_interface);
        }
    }

    return own_interfaces;
}


void ELEMENT::ComputeInternalU(int u_flag) {
	for (int i = 0; i < this->number_gp_internal; i++) {
		this->u_internal[u_flag][i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < this->number_gp_internal; j++) {
			this->u_internal[u_flag][j] += this->u[u_flag][i] * this->phi_internal[i][j];
		}
	}
}

void ELEMENT::ComputeInternalDU(int dim, int u_flag, int u_flag_store) {
	for (int i = 0; i < this->number_gp_boundary; i++) {
		this->u_internal[u_flag_store][i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < this->number_gp_internal; j++) {
			this->u_internal[u_flag_store][j] += this->u[u_flag][i] * this->internal_dphi_fac[i][dim][j];
		}
	}
}

void ELEMENT::ComputeBoundaryU(int u_flag) {
	for (int i = 0; i < this->number_boundaries; i++) {
		for (int j = 0; j < this->number_gp_boundary; j++) {
			this->u_boundary[i][u_flag][j] = 0.0;
		}

		for (int j = 0; j < this->number_bf; j++) {
			for (int k = 0; k < this->number_gp_boundary; k++) {
				this->u_boundary[i][u_flag][k] += this->u[u_flag][j] * this->phi_boundary[i][j][k];
			}
		}
	}
}

double ELEMENT::IntegrationInternalPhi(int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->number_gp_internal; i++) {
		integral += this->u_internal[u_flag][i] * this->internal_int_fac_phi[phi_n][i];
	}

	return integral;
}

double ELEMENT::IntegrationInternalDPhi(int dim, int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->number_gp_internal; i++) {
		integral += this->u_internal[u_flag][i] * this->internal_int_fac_dphi[phi_n][dim][i];
	}

	return integral;
}

double ELEMENT::IntegrationBoundaryPhi(int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->number_boundaries; i++) {
		for (int j = 0; j < this->number_gp_boundary; j++) {
			integral += this->u_boundary[i][u_flag][j] * this->boundary_int_fac_phi[i][phi_n][j];
		}
	}

	return integral;
}

void ELEMENT::SolveLSE(int u_flag) {
	if (this->m_inv.first) { //diagonal
		for (int i = 0; i < this->number_bf; i++) {
			this->u[u_flag][i] = this->m_inv.second[0][i] * this->RHS[i];
		}
	}
	else if (!(this->m_inv.first)) { //not diagonal
		for (int i = 0; i < this->number_bf; i++) {
			this->u[u_flag][i] = 0;
			for (int j = 0; j < this->number_bf; j++) {
				this->u[u_flag][i] += this->m_inv.second[i][j] * this->RHS[j];
			}
		}
	}
}

void ELEMENT::InitializeVTK(std::vector<double*>& points, std::vector<unsigned int*>& cells) {
	switch (this->element_type) {
	case TRIANGLE: this->InitializeVTKTriangle(points, cells); break;
	default:
		printf("\n");
		printf("ELEMENT_2D InitializeVTK - Fatal error!\n");
		printf("Undefined element type = %d\n", this->element_type);
		exit(1);
	}
}

void ELEMENT::WriteCellDataVTK(std::vector<double>& cell_data, int u_flag) {
	switch (element_type) {
	case TRIANGLE: this->WriteCellDataVTKTriangle(cell_data, u_flag); break;
	default:
		printf("\n");
		printf("ELEMENT_2D WriteCellDataVTK - Fatal error!\n");
		printf("Undefined element type = %d\n", this->element_type);
		exit(1);
	}
}

void ELEMENT::WritePointDataVTK(std::vector<double>& point_data, int u_flag) {
	switch (element_type) {
	case TRIANGLE: this->WritePointDataVTKTriangle(point_data, u_flag); break;
	default:
		printf("\n");
		printf("ELEMENT_2D WritePointDataVTK - Fatal error!\n");
		printf("Undefined element type = %d\n", this->element_type);
		exit(1);
	}
}