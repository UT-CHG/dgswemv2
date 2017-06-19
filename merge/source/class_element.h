#ifndef CLASS_ELEMENT_H
#define CLASS_ELEMENT_H

#include "general_definitions.h"

#include "class_interface.h"
#include "element/class_master_element.h"
#include "shape/shapes_2D.h"
#include "boundary/class_boundary.h"

template<int dimension = 2, int element_type = TRIANGLE, 
	class basis_type = Basis::Dubiner_2D, 
	class integration_int_type = Integration::Dunavant_2D, 
	class integration_bound_type = Integration::GaussLegendre_1D,
	class shape_type = Shape::StraightTriangle>
class Element{
    friend class PROBLEM;

private:
	MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>& master;
	shape_type& shape;

    unsigned int ID;
	unsigned char number_boundaries;

   	std::vector<unsigned char> boundary_type;
	std::vector<unsigned int> neighbor_ID;
	std::pair<std::vector<bool>, std::vector<INTERFACE*>> interfaces;

	std::vector<Point<dimension>> nodal_coordinates;
    
	int number_bf;

    int number_gp_internal;
    int number_gp_boundary;

	Array3D<double> internal_dphi_fac;	
    Array2D<double> internal_int_fac_phi;
	Array3D<double> internal_int_fac_dphi;
	Array3D<double> boundary_int_fac_phi;
	
	SWE::Data data;

	Array2D<double> u;
    Array3D<double> u_substep;

    Array2D<double> u_gp;
	Array3D<double> u_boundary_;
    double*** u_boundary;

    std::vector<double> RHS;
	std::pair<bool, Array2D<double>> m_inv;

public:
	Element(MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>&,
		shape_type&, unsigned int, std::vector<unsigned int>&, std::vector<unsigned char>&, std::vector<Point<dimension>>&);

    std::map<unsigned int, INTERFACE*> CreateInterfaces();
    void AppendInterface(unsigned int, INTERFACE*);
    std::vector<std::pair<unsigned char, INTERFACE*>> GetOwnInterfaces();

	std::vector<RawBoundary<dimension - 1>*> CreateBoundaries();

    void ComputeInternalU(int);
	void ComputeInternalDU(int, int, int);
    void ComputeBoundaryU(int);
   
    double IntegrationInternalPhi(int, int);
    double IntegrationInternalDPhi(int, int, int);

    double IntegrationBoundaryPhi(int, int);

    void SolveLSE(int);

    void InitializeVTK(std::vector<Point<3>>&, Array2D<unsigned int>&);
    void WriteCellDataVTK(std::vector<double>&, int);
	void WritePointDataVTK(std::vector<double>&, int);

private:
	void allocate_memory();
};

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::Element(MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>& master,
	shape_type& shape, unsigned int ID, std::vector<unsigned int>& neighbor_ID, std::vector<unsigned char>& boundary_type,
	std::vector<Point<dimension>>& nodal_coordinates) :

	master(master), shape(shape),

	//these values are stored at element to avoid repeat access to the values from master element
	number_boundaries(master.number_boundaries), number_bf(master.number_bf),
	number_gp_boundary(master.number_gp_boundary), number_gp_internal(master.number_gp_internal),

	ID(ID), nodal_coordinates(nodal_coordinates),
	neighbor_ID(neighbor_ID), boundary_type(boundary_type)
{
	this->allocate_memory();

	for (int i = 0; i < this->number_boundaries; i++) {
		this->interfaces.first[i] = false;
		this->interfaces.second[i] = nullptr;
	}

	//DEFORMATION
	std::vector<double> det_J = shape.get_J_det(this->nodal_coordinates);
	Array3D<double> J_inv = shape.get_J_inv(this->nodal_coordinates);
	Array2D<double> surface_J = shape.get_surface_J(this->nodal_coordinates);;

	if (det_J.size() == 1) { //constant Jacobian
							 //DIFFERENTIATION FACTORS
		this->internal_dphi_fac.resize(this->master.dphi_internal.size());
		for (int i = 0; i < this->master.dphi_internal.size(); i++) {
			this->internal_dphi_fac[i].resize(dimension);
			for (int j = 0; j < dimension; j++) {
				this->internal_dphi_fac[i][j].reserve(this->master.dphi_internal[i][j].size());
				for (int k = 0; k < this->master.dphi_internal[i][j].size(); k++) {
					double dphi = 0;
					for (int l = 0; l < dimension; l++) {
						dphi += this->master.dphi_internal[i][l][k] * J_inv[l][j][0];
					}
					this->internal_dphi_fac[i][j].push_back(dphi);
				}
			}
		}

		//INTEGRATION FACTORS
		this->internal_int_fac_phi = this->master.internal_int_fac_phi;

		this->internal_int_fac_dphi.resize(this->master.internal_int_fac_dphi.size());
		for (int i = 0; i < this->master.internal_int_fac_dphi.size(); i++) {
			this->internal_int_fac_dphi[i].resize(dimension);
			for (int j = 0; j < dimension; j++) {
				this->internal_int_fac_dphi[i][j].reserve(this->master.internal_int_fac_dphi[i][j].size());
				for (int k = 0; k < this->master.internal_int_fac_dphi[i][j].size(); k++) {
					double int_dphi = 0;
					for (int l = 0; l < dimension; l++) {
						int_dphi += this->master.internal_int_fac_dphi[i][l][k] * J_inv[l][j][0];
					}
					this->internal_int_fac_dphi[i][j].push_back(int_dphi);
				}
			}
		}

		this->boundary_int_fac_phi = this->master.boundary_int_fac_phi;
		for (int i = 0; i < this->master.boundary_int_fac_phi.size(); i++) {
			for (int j = 0; j < this->master.boundary_int_fac_phi[i].size(); j++) {
				for (int k = 0; k < this->master.boundary_int_fac_phi[i][j].size(); k++) {
					this->boundary_int_fac_phi[i][j][k] *= surface_J[i][0] / abs(det_J[0]);
				}
			}
		}

		//MASS MATRIX
		this->m_inv = this->master.m_inv;
	}
	else {
		//Placeholder for cases p_geom > 1
	}

	this->data.state_data.push_back(SWE::State(this->internal_int_fac_phi.size()));
	this->data.internal_data = SWE::Internal((*this->internal_int_fac_phi.begin()).size());

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
		0.3*cos(beta * this->nodal_coordinates[0][X]) / cos(beta * L),
		0.3*cos(beta * this->nodal_coordinates[1][X]) / cos(beta * L),
		0.3*cos(beta * this->nodal_coordinates[2][X]) / cos(beta * L),
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

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
void Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::allocate_memory() {
	//INITIALIZE ARRAYS TO STORE INTERFACE DATA
	this->interfaces.first.resize(this->number_boundaries);
	this->interfaces.second.resize(this->number_boundaries);

	//INITIALIZE ARRAYS TO STORE Us and RHS
	this->RHS.resize(this->number_bf);

	this->u.resize(SIZE_U);
	for (int i = 0; i < SIZE_U; i++) {
		this->u[i].resize(this->number_bf);
	}

	this->u_gp.resize(SIZE_U_INTERNAL);
	for (int i = 0; i < SIZE_U_INTERNAL; i++) {
		this->u_gp[i].resize(this->number_gp_internal);
	}

	this->u_boundary_.resize(this->number_boundaries);
	for (int i = 0; i < this->number_boundaries; i++) {
		this->u_boundary_[i].resize(SIZE_U_BOUNDARY);
		for (int j = 0; j < SIZE_U_BOUNDARY; j++) {
			this->u_boundary_[i][j].resize(this->number_gp_boundary);
		}
	}

	this->u_boundary = new double**[this->number_boundaries];
	for (int i = 0; i < this->number_boundaries; i++) {
		this->u_boundary[i] = new double*[SIZE_U_BOUNDARY];
		for (int j = 0; j < SIZE_U_BOUNDARY; j++) {
			this->u_boundary[i][j] = new double[this->number_gp_boundary];
		}
	}
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
std::vector<RawBoundary<dimension - 1>*> 
Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::CreateBoundaries() {
	std::vector<RawBoundary<dimension - 1>*> my_raw_boundaries;

	Basis::Basis<dimension>* basis;

	std::function<std::vector<Point<dimension>>(const std::vector<Point<dimension - 1>>&)> boundary_to_master;
	std::function<Array2D<double>()> get_surface_normal;
	std::function<std::vector<double>()> get_surface_J;

	for (int i = 0; i < this->boundary_type.size(); i++) {
		basis = (Basis::Basis<dimension>*)(&master.basis);

		boundary_to_master = std::bind(&MasterElement<dimension, element_type, basis_type,
			integration_int_type, integration_bound_type>::TriangleBoundaryToMasterCoordinates,
			master, i, std::placeholders::_1); //need specialization in master element

		get_surface_normal = std::bind(&Shape::StraightTriangle::get_surface_normal_,
			shape, i, nodal_coordinates);

		get_surface_J = std::bind(&Shape::StraightTriangle::get_surface_J_,
			shape, i, nodal_coordinates);

		my_raw_boundaries.push_back(new RawBoundary<dimension-1>(this->boundary_type[i], this->neighbor_ID[i], master.p, 
			this->data, basis, boundary_to_master, get_surface_normal, get_surface_J));
	}

	return my_raw_boundaries;
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
std::map<unsigned int, INTERFACE*> 
Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::CreateInterfaces() {
	std::map<unsigned int, INTERFACE*> internal_interfaces;

	std::vector<double> det_J = shape.get_J_det(this->nodal_coordinates);
	Array3D<double> surface_normal = shape.get_surface_normal(this->nodal_coordinates);

	bool boundary;
	for (int i = 0; i < this->interfaces.second.size(); i++) {
		if (this->interfaces.second[i] == nullptr) {
			if (this->neighbor_ID[i] == DEFAULT_ID) boundary = true;
			else if (this->neighbor_ID[i] != DEFAULT_ID) boundary = false;

			double** normal = new double*[dimension];
			for (int j = 0; j < dimension; j++) {
				normal[j] = new double[this->number_gp_boundary];
			}

			if (det_J.size() == 1) {
				for (int j = 0; j < dimension; j++) {
					for (int k = 0; k < this->number_gp_boundary; k++) {
						normal[j][k] = surface_normal[i][j][0];
					}
				}
			}
			else {
				for (int j = 0; j < dimension; j++) {
					for (int k = 0; k < this->number_gp_boundary; k++) {
						normal[j][k] = surface_normal[i][j][k];
					}
				}
			}

			this->interfaces.first[i] = true;

			this->interfaces.second[i] = new INTERFACE(dimension, this->number_gp_boundary,
				this->u_boundary[i], normal, boundary);

			if (this->neighbor_ID[i] != DEFAULT_ID) {
				internal_interfaces[this->neighbor_ID[i]] = this->interfaces.second[i];
			}
		}
	}

	return internal_interfaces;
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
void Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::AppendInterface(unsigned int neighbor_ID, INTERFACE* interface_ptr) {
	for (int i = 0; i < this->number_boundaries; i++) {
		if (this->neighbor_ID[i] == neighbor_ID) {
			this->interfaces.second[i] = interface_ptr;

			this->interfaces.second[i]->SetPointerEX(this->u_boundary[i]);
		}
	}
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
std::vector<std::pair<unsigned char, INTERFACE*>> 
Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::GetOwnInterfaces() {
	std::vector<std::pair<unsigned char, INTERFACE*>> own_interfaces;

	for (int i = 0; i < this->number_boundaries; i++) {
		if (this->interfaces.first[i]) {
			std::pair<unsigned char, INTERFACE*> own_interface;

			own_interface.first = this->boundary_type[i];
			own_interface.second = this->interfaces.second[i];
			own_interfaces.push_back(own_interface);
		}
	}

	return own_interfaces;
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
void Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::ComputeInternalU(int u_flag) {
	for (int i = 0; i < this->u_gp[u_flag].size(); i++) {
		this->u_gp[u_flag][i] = 0.0;
	}

	for (int i = 0; i < this->u[u_flag].size(); i++) {
		for (int j = 0; j < this->u_gp[u_flag].size(); j++) {
			this->u_gp[u_flag][j] += this->u[u_flag][i] * this->master.phi_internal[i][j];
		}
	}
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
void Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::ComputeInternalDU(int dir, int u_flag, int u_flag_store) {
	for (int i = 0; i < this->u_gp[u_flag_store].size(); i++) {
		this->u_gp[u_flag_store][i] = 0.0;
	}

	for (int i = 0; i < this->u[u_flag].size(); i++) {
		for (int j = 0; j < this->u_gp[u_flag_store].size(); j++) {
			this->u_gp[u_flag_store][j] += this->u[u_flag][i] * this->internal_dphi_fac[i][dir][j];
		}
	}
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
void Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::ComputeBoundaryU(int u_flag) {
	for (int i = 0; i < this->u_boundary_.size(); i++) {
		for (int j = 0; j < this->u_boundary_[i][u_flag].size(); j++) {
			this->u_boundary[i][u_flag][j] = 0.0;
		}

		for (int j = 0; j < this->u[u_flag].size(); j++) {
			for (int k = 0; k < this->u_boundary_[i][u_flag].size(); k++) {
				this->u_boundary[i][u_flag][k] += this->u[u_flag][j] * this->master.phi_boundary[i][j][k];
			}
		}
	}
}

template<int dimension, int element_type, class basis_type, class integration_int_type,
	class integration_bound_type, class shape_type>
double Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::IntegrationInternalPhi(int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->internal_int_fac_phi[phi_n].size(); i++) {
		integral += this->u_gp[u_flag][i] * this->internal_int_fac_phi[phi_n][i];
	}

	return integral;
}

template<int dimension, int element_type, class basis_type, class integration_int_type, 
	class integration_bound_type, class shape_type>
double Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::IntegrationInternalDPhi(int dir, int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->internal_int_fac_dphi[phi_n][dir].size(); i++) {
		integral += this->u_gp[u_flag][i] * this->internal_int_fac_dphi[phi_n][dir][i];
	}

	return integral;
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
double Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::IntegrationBoundaryPhi(int u_flag, int phi_n) {
	double integral = 0;

	for (int i = 0; i < this->boundary_int_fac_phi.size(); i++) {
		for (int j = 0; j < this->boundary_int_fac_phi[i][phi_n].size(); j++) {
			integral += this->u_boundary[i][u_flag][j] * this->boundary_int_fac_phi[i][phi_n][j];
		}
	}

	return integral;
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
void Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::SolveLSE(int u_flag) {
	if (this->m_inv.first) { //diagonal
		for (int i = 0; i < this->u[u_flag].size(); i++) {
			this->u[u_flag][i] = this->m_inv.second[0][i] * this->RHS[i];
		}
	}
	else if (!(this->m_inv.first)) { //not diagonal
		for (int i = 0; i < this->u[u_flag].size(); i++) {
			this->u[u_flag][i] = 0;
			for (int j = 0; j < this->u[u_flag].size(); j++) {
				this->u[u_flag][i] += this->m_inv.second[i][j] * this->RHS[j];
			}
		}
	}
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
void Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::InitializeVTK(std::vector<Point<3>>& points, Array2D<unsigned int>& cells) {
	this->shape.get_VTK(points, cells, nodal_coordinates);
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
void Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::WriteCellDataVTK(std::vector<double>& cell_data, int u_flag) {
	Array2D<double> temp = this->master.phi_postprocessor_cell;

	for (int i = 0; i < temp.size(); i++) {
		for (int j = 0; j < temp[i].size(); j++) {
			temp[i][j] *= this->u[u_flag][i];
		}
	}

	for (int i = 1; i < temp.size(); i++) {
		for (int j = 0; j < temp[i].size(); j++) {
			temp[0][j] += temp[i][j];
		}
	}

	cell_data.insert(cell_data.end(), temp[0].begin(), temp[0].end());
}

template<int dimension, int element_type, class basis_type, 
	class integration_int_type, class integration_bound_type, class shape_type>
void Element<dimension, element_type, basis_type, integration_int_type, integration_bound_type, shape_type>
::WritePointDataVTK(std::vector<double>& point_data, int u_flag) {
	Array2D<double> temp = this->master.phi_postprocessor_point;

	for (int i = 0; i < temp.size(); i++) {
		for (int j = 0; j < temp[i].size(); j++) {
			temp[i][j] *= this->u[u_flag][i];
		}
	}

	for (int i = 1; i < temp.size(); i++) {
		for (int j = 0; j < temp[i].size(); j++) {
			temp[0][j] += temp[i][j];
		}
	}

	point_data.insert(point_data.end(), temp[0].begin(), temp[0].end());
}

#endif