#include "class_element_2D.h"

ELEMENT_2D::ELEMENT_2D(int element_type, unsigned int ID, unsigned int* neighbor_ID, unsigned char* boundary_type,
    double* nodal_coordinates_x, double* nodal_coordinates_y,
    BASIS* basis, BASIS_GEOM* basis_geom) : ELEMENT(ID) {
    this->basis = basis;

    if (basis_geom != nullptr) {
        this->basis_geom = basis_geom;
    }

	switch (element_type){
	case TRIANGLE: this->Triangle(neighbor_ID, boundary_type, nodal_coordinates_x, nodal_coordinates_y); break;
	default:
		printf("\n");
		printf("ELEMENT_2D CONSTRUCTOR - Fatal error!\n");
		printf("Undefined element type = %d\n", element_type);
		exit(1);
	}

    //COMPUTE DIFFERENTIATION AND NUMERICAL INTEGRATION FACTORS
	this->ComputeDPhi();
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

ELEMENT_2D::~ELEMENT_2D() {
	for (int i = 0; i < this->dimension; i++) {
		delete[] this->nodal_coordinates[i];
	}
	delete[] this->nodal_coordinates;

	delete[] this->boundary_type;
	delete[] this->neighbor_ID;
    
    /*for (int i = 0; i < this->number_interfaces; i++) {
    if (this->interface_owner[i]) {
    delete this->interfaces[i];
    }
    }*///at present it deletes all interfaces, however it should keep some for its neighbors
    //delete through mesh class

    delete[] this->interfaces;
    delete[] this->interface_owner;

	delete[] this->RHS;

    for (int i = 0; i < SIZE_U; i++) {
        delete[] this->u[i];
    }
	delete[] this->u;

    for (int i = 0; i < SIZE_U_INTERNAL; i++) {
        delete[] this->u_internal[i];
    }
	delete[] this->u_internal;

    for (int i = 0; i < this->number_boundaries; i++) {
        for (int j = 0; j < SIZE_U_BOUNDARY; j++) {
            delete[] this->u_boundary[i][j];
        }
        delete[] this->u_boundary[i];
    }
    delete[] this->u_boundary;
    
	for (int i = 0; i < this->number_bf; i++) {
		delete[] this->internal_int_fac_phi[i];
	}
	delete[] this->internal_int_fac_phi;

	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->number_bf; j++) {
			delete[] this->internal_int_fac_dphi[i][j];
			delete[] this->dphi_internal[i][j];
		}
		delete[] this->internal_int_fac_dphi[i];
		delete[] this->dphi_internal[i];
	}
	delete[] this->internal_int_fac_dphi;
	delete[] this->dphi_internal;

	for (int i = 0; i < this->number_boundaries; i++) {
		for (int j = 0; j < this->number_bf; j++) {
			delete[] this->boundary_int_fac_phi[i][j];
		}
		delete[] this->boundary_int_fac_phi[i];
    }
	delete[] this->boundary_int_fac_phi;

	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->number_boundaries; j++) {
			for (int k = 0; k < this->number_bf; k++) {
				delete[] this->boundary_int_fac_n[i][j][k];
			}
			delete[] this->boundary_int_fac_n[i][j];
		}
		delete[] this->boundary_int_fac_n[i];
	}
    delete[] this->boundary_int_fac_n;
}

void ELEMENT_2D::allocate_memory() {
	//INITIALIZE ARRAYS TO STORE ESSENTIAL ELEMENT DATA
	this->nodal_coordinates = new double*[this->dimension];
	for (int i = 0; i < this->dimension; i++) {
		this->nodal_coordinates[i] = new double[this->number_bf_geom];
	}

	this->boundary_type = new unsigned char[this->number_boundaries];
	this->neighbor_ID = new unsigned int[this->number_boundaries];
	
	this->interfaces = new INTERFACE_2D*[this->number_boundaries];
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

	//INITIALIZE ARRAYS TO STORE INTEGRATION AND DIFFERENTIATION FACTORS
	this->internal_int_fac_phi = new double*[this->number_bf];
	for (int i = 0; i < this->number_bf; i++) {
		this->internal_int_fac_phi[i] = new double[this->number_gp_internal];
	}

	this->internal_int_fac_dphi = new double**[this->dimension];
	this->dphi_internal = new double**[this->dimension];
	for (int i = 0; i < this->dimension; i++) {
		this->internal_int_fac_dphi[i] = new double*[this->number_bf];
		this->dphi_internal[i] = new double*[this->number_bf];
		for (int j = 0; j < this->number_bf; j++) {
			this->internal_int_fac_dphi[i][j] = new double[this->number_gp_internal];
			this->dphi_internal[i][j] = new double[this->number_gp_internal];
		}
	}

	this->boundary_int_fac_phi = new double**[this->number_boundaries];
	for (int i = 0; i < this->number_boundaries; i++) {
		this->boundary_int_fac_phi[i] = new double*[this->number_bf];
		for (int j = 0; j < this->number_bf; j++) {
			this->boundary_int_fac_phi[i][j] = new double[this->number_gp_boundary];
		}
	}

	this->boundary_int_fac_n = new double***[this->dimension];
	for (int i = 0; i < this->dimension; i++) {
		this->boundary_int_fac_n[i] = new double**[this->number_boundaries];
		for (int j = 0; j < this->number_boundaries; j++) {
			this->boundary_int_fac_n[i][j] = new double*[this->number_bf];
			for (int k = 0; k < this->number_bf; k++) {
				this->boundary_int_fac_n[i][j][k] = new double[this->number_gp_boundary];
			}
		}
	}
	//INITIALIZE ARRAYS TO STORE GEOMETRIC DATA
	if (this->basis_geom == nullptr) {
		this->det_J_internal = new double[1];

		this->J_inv_t_internal = new double**[this->dimension];
		for (int i = 0; i < this->dimension; i++) {
			this->J_inv_t_internal[i] = new double*[this->dimension];
			for (int j = 0; j < this->dimension; j++) {
				this->J_inv_t_internal[i][j] = new double[1];
			}
		}

		this->surface_J_boundary = new double*[this->number_boundaries];
		for (int i = 0; i < this->number_boundaries; i++) {
			this->surface_J_boundary[i] = new double[1];
		}

		this->normal_boundary = new double**[this->dimension];
		for (int i = 0; i < this->dimension; i++) {
			this->normal_boundary[i] = new double*[this->number_boundaries];
			for (int j = 0; j < this->number_boundaries; j++) {
				this->normal_boundary[i][j] = new double[1];
			}
		}
	}
	else {
		//Placeholder for cases p_geom > 1
	}
}

void ELEMENT_2D::ComputeDPhi() {
	double*** dphi_dz_internal = this->basis->GetDPhiDZInternal();

	if (this->basis_geom == nullptr) {
		for (int i = 0; i < this->dimension; i++) {
			for (int j = 0; j < this->number_bf; j++) {
				for (int k = 0; k < this->number_gp_internal; k++) {
					this->dphi_internal[i][j][k] = 0;
					for (int l = 0; l < this->dimension; l++) {
						this->dphi_internal[i][j][k] += this->J_inv_t_internal[i][l][0] * dphi_dz_internal[l][j][k];
					}
				}
			}
		}
	}
	else {
		//Placeholder for cases p_geom > 1
	}
}

void ELEMENT_2D::ComputeIntegrationFactors() {
	this->phi_internal = this->basis->GetPhiInternal();
	this->phi_boundary = this->basis->GetPhiBoundary();

    double* w_internal = this->basis->GetIntegrationRuleInternal()->GetWeight();
    double* w_boundary = this->basis->GetIntegrationRuleBoundary()->GetWeight();

    if (this->basis_geom == nullptr) {
        for (int i = 0; i < this->number_bf; i++) {
            for (int j = 0; j < this->number_gp_internal; j++) {
                this->internal_int_fac_phi[i][j] = this->phi_internal[i][j] * w_internal[j];
            }
        }

		for (int i = 0; i < this->dimension; i++) {
			for (int j = 0; j < this->number_bf; j++) {
				for (int k = 0; k < this->number_gp_internal; k++) {
					this->internal_int_fac_dphi[i][j][k] = this->dphi_internal[i][j][k] * w_internal[j];
				}
			}
		}

        for (int i = 0; i < this->number_boundaries; i++) {
            for (int j = 0; j < this->number_bf; j++) {
                for (int k = 0; k < this->number_gp_boundary; k++) {
                    this->boundary_int_fac_phi[i][j][k] = this->phi_boundary[i][j][k] *
                        w_boundary[k] * this->surface_J_boundary[i][0] / abs(this->det_J_internal[0]);
                }
            }
        }

		for (int i = 0; i < this->dimension; i++) {
			for (int j = 0; j < this->number_boundaries; j++) {
				for (int k = 0; k < this->number_bf; k++) {
					for (int l = 0; l < this->number_gp_boundary; l++) {
						this->boundary_int_fac_n[i][j][k][l] = this->phi_boundary[j][k][l] * this->normal_boundary[i][j][0] *
							w_boundary[l] * this->surface_J_boundary[j][0] / abs(this->det_J_internal[0]);
					}
				}
			}
		}
	}
    else {
        //Placeholder for cases p_geom > 1
    }
}

std::map<unsigned int, INTERFACE*> ELEMENT_2D::CreateInterfaces() {
	std::map<unsigned int, INTERFACE*> internal_interfaces;

	for (int i = 0; i < this->number_boundaries; i++) {
		if (this->interfaces[i] == nullptr) {
			bool straight;
			if (this->basis_geom == nullptr) straight = true;
			else if (this->basis_geom != nullptr) straight = false;

			bool boundary;
			if (this->neighbor_ID[i] == DEFAULT_ID) boundary = true;
			else if (this->neighbor_ID[i] != DEFAULT_ID) boundary = false;

			this->interfaces[i] = new INTERFACE_2D(this->number_gp_boundary, this->u_boundary[i],
				this->normal_boundary[X][i], this->normal_boundary[Y][i], straight, boundary);

			this->interface_owner[i] = true;

			if (this->neighbor_ID[i] != DEFAULT_ID) {
				internal_interfaces[this->neighbor_ID[i]] = this->interfaces[i];
			}
		}
	}

	//DELETE NO LONGER NEEDED GEOMETRIC DATA
	delete[] this->det_J_internal;

	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->dimension; j++) {
			delete[] this->J_inv_t_internal[i][j];
		}
		delete[] this->J_inv_t_internal[i];
	}
	delete[] this->J_inv_t_internal;

	for (int i = 0; i < this->number_boundaries; i++) {
		delete[] this->surface_J_boundary[i];
	}
	delete[] this->surface_J_boundary;

	for (int i = 0; i < this->dimension; i++) {
		for (int j = 0; j < this->number_boundaries; j++) {
			delete[] this->normal_boundary[i][j];
		}
		delete[] this->normal_boundary[i];
	}
	delete[] this->normal_boundary;

	return internal_interfaces;
}

void ELEMENT_2D::AppendInterface(unsigned int neighbor_ID, INTERFACE* interface_ptr) {
    for (int i = 0; i < this->number_boundaries; i++) {
        if (this->neighbor_ID[i] == neighbor_ID) {
            this->interfaces[i] = (INTERFACE_2D*)interface_ptr;

            this->interfaces[i]->SetPointerEX(this->u_boundary[i]);
        }
    }
}

std::vector<std::pair<unsigned char, INTERFACE*>> ELEMENT_2D::GetOwnInterfaces() {
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

void ELEMENT_2D::ComputeInternalU(int u_flag) {
    for (int i = 0; i < this->number_gp_internal; i++) {
        this->u_internal[u_flag][i] = 0.0;
    }

    for (int i = 0; i < this->number_bf; i++) {
        for (int j = 0; j < this->number_gp_internal; j++) {
            this->u_internal[u_flag][j] += this->u[u_flag][i] * this->phi_internal[i][j];
        }
    }
}

void ELEMENT_2D::ComputeBoundaryU(int u_flag) {
    for (int k = 0; k < this->number_boundaries; k++) {
        for (int i = 0; i < this->number_gp_boundary; i++) {
            this->u_boundary[k][u_flag][i] = 0.0;
        }

        for (int i = 0; i < this->number_bf; i++) {
            for (int j = 0; j < this->number_gp_boundary; j++) {
                this->u_boundary[k][u_flag][j] += this->u[u_flag][i] * this->phi_boundary[k][i][j];
            }
        }
    }
}

void ELEMENT_2D::ComputeInternalDUDX(int u_flag, int u_flag_store) {
	for (int i = 0; i < this->number_gp_internal; i++) {
		this->u_internal[u_flag_store][i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < this->number_gp_internal; j++) {
			this->u_internal[u_flag_store][j] += this->u[u_flag][i] * this->dphi_internal[X][i][j];
		}
	}
}

void ELEMENT_2D::ComputeInternalDUDY(int u_flag, int u_flag_store) {
	for (int i = 0; i < this->number_gp_internal; i++) {
		this->u_internal[u_flag_store][i] = 0.0;
	}

	for (int i = 0; i < this->number_bf; i++) {
		for (int j = 0; j < this->number_gp_internal; j++) {
			this->u_internal[u_flag_store][j] += this->u[u_flag][i] * this->dphi_internal[Y][i][j];
		}
	}
}

double ELEMENT_2D::IntegrationInternalPhi(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_internal; i++) {
        integral += this->u_internal[u_flag][i] * this->internal_int_fac_phi[phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationInternalDPhiDX(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_internal; i++) {
        integral += this->u_internal[u_flag][i] * this->internal_int_fac_dphi[X][phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationInternalDPhiDY(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_gp_internal; i++) {
        integral += this->u_internal[u_flag][i] * this->internal_int_fac_dphi[Y][phi_n][i];
    }

    return integral;
}

double ELEMENT_2D::IntegrationBoundaryPhi(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_boundaries; i++) {
        for (int j = 0; j < this->number_gp_boundary; j++) {
            integral += this->u_boundary[i][u_flag][j] * this->boundary_int_fac_phi[i][phi_n][j];
        }
    }

    return integral;
}

double ELEMENT_2D::IntegrationBoundaryNX(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_boundaries; i++) {
        for (int j = 0; j < this->number_gp_boundary; j++) {
            integral += this->u_boundary[i][u_flag][j] * this->boundary_int_fac_n[X][i][phi_n][j];
        }
    }

    return integral;
}

double ELEMENT_2D::IntegrationBoundaryNY(int u_flag, int phi_n) {
    double integral = 0;

    for (int i = 0; i < this->number_boundaries; i++) {
        for (int j = 0; j < this->number_gp_boundary; j++) {
            integral += this->u_boundary[i][u_flag][j] * this->boundary_int_fac_n[Y][i][phi_n][j];
        }
    }

    return integral;
}

void ELEMENT_2D::SolveLSE(int u_flag) {
    if (this->orthogonal) {
        for (int i = 0; i < this->number_bf; i++) {
            this->u[u_flag][i] = this->m_inv[0][i] * this->RHS[i];
        }
    }
    else if (!(this->orthogonal)) {
        for (int i = 0; i < this->number_bf; i++) {
            this->u[u_flag][i] = 0;
            for (int j = 0; j < this->number_bf; j++) {
                this->u[u_flag][i] += this->m_inv[i][j] * this->RHS[j];
            }
        }
    }
}

void ELEMENT_2D::InitializeVTK(std::vector<double*>& points, std::vector<unsigned int*>& cells) {
	switch (this->element_type) {
	case TRIANGLE: this->InitializeVTKTriangle(points, cells); break;
	default:
		printf("\n");
		printf("ELEMENT_2D InitializeVTK - Fatal error!\n");
		printf("Undefined element type = %d\n", this->element_type);
		exit(1);
	}
}

void ELEMENT_2D::WriteCellDataVTK(std::vector<double>& cell_data, int u_flag) {
	switch (this->element_type) {
	case TRIANGLE: this->WriteCellDataVTKTriangle(cell_data, u_flag); break;
	default:
		printf("\n");
		printf("ELEMENT_2D WriteCellDataVTK - Fatal error!\n");
		printf("Undefined element type = %d\n", this->element_type);
		exit(1);
	}
}

void ELEMENT_2D::WritePointDataVTK(std::vector<double>& point_data, int u_flag) {
	switch (this->element_type) {
	case TRIANGLE: this->WritePointDataVTKTriangle(point_data, u_flag); break;
	default:
		printf("\n");
		printf("ELEMENT_2D WritePointDataVTK - Fatal error!\n");
		printf("Undefined element type = %d\n", this->element_type);
		exit(1);
	}
}