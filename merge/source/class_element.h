#ifndef CLASS_ELEMENT_H
#define CLASS_ELEMENT_H

#include "general_definitions.h"

#include "class_boundary.h"

#include "master/master_elements_2D.h"
#include "shape/shapes_2D.h"

template<int dimension = 2, 
	class master_type = Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>,
	class shape_type = Shape::StraightTriangle>
class Element {
public:
	SWE::Data data;

private:
	unsigned int ID;

	master_type& master;
	shape_type& shape;

	std::vector<unsigned char> boundary_type;
	std::vector<unsigned int> neighbor_ID;

	std::vector<Point<dimension>> nodal_coordinates;

	Array3D<double> dphi_fact;
	Array2D<double> int_fact_phi;
	Array3D<double> int_fact_dphi;

	std::pair<bool, Array2D<double>> m_inv;

public:
	Element(master_type&, shape_type&, unsigned int, std::vector<unsigned int>&, 
		std::vector<unsigned char>&, std::vector<Point<dimension>>&);

	std::vector<RawBoundary<dimension - 1>*> CreateBoundaries();

	void ComputeUgp(const std::vector<double>&, std::vector<double>&);
	void ComputeDUgp(int, const std::vector<double>&, std::vector<double>&);
	double IntegrationPhi(int, const std::vector<double>&);
	double IntegrationDPhi(int, int, const std::vector<double>&);

	std::vector<double> SolveLSE(const std::vector<double>&);

	void InitializeVTK(std::vector<Point<3>>&, Array2D<unsigned int>&);
	void WriteCellDataVTK(std::vector<double>&, int);
	void WritePointDataVTK(std::vector<double>&, int);
};

template<int dimension, class master_type, class shape_type>
Element<dimension, master_type, shape_type>::Element(master_type& master, shape_type& shape, unsigned int ID, 
	std::vector<unsigned int>& neighbor_ID, std::vector<unsigned char>& boundary_type, std::vector<Point<dimension>>& nodal_coordinates) :

	ID(ID), master(master), shape(shape), nodal_coordinates(nodal_coordinates),
	neighbor_ID(neighbor_ID), boundary_type(boundary_type)
{
	//DEFORMATION
	std::vector<double> det_J = shape.get_J_det(this->nodal_coordinates);
	Array3D<double> J_inv = shape.get_J_inv(this->nodal_coordinates);

	if (det_J.size() == 1) { //constant Jacobian
		//DIFFERENTIATION FACTORS
		this->dphi_fact.resize(this->master.dphi_gp.size());
		for (int i = 0; i < this->master.dphi_gp.size(); i++) {
			this->dphi_fact[i].resize(dimension);
			for (int j = 0; j < dimension; j++) {
				this->dphi_fact[i][j].reserve(this->master.dphi_gp[i][j].size());
				for (int k = 0; k < this->master.dphi_gp[i][j].size(); k++) {
					double dphi = 0;
					for (int l = 0; l < dimension; l++) {
						dphi += this->master.dphi_gp[i][l][k] * J_inv[l][j][0];
					}
					this->dphi_fact[i][j].push_back(dphi);
				}
			}
		}

		//INTEGRATION FACTORS
		this->int_fact_phi = this->master.int_fact_phi;
		for (int i = 0; i < this->int_fact_phi.size(); i++) {
			for (int j = 0; j < this->int_fact_phi[i].size(); j++) {
				this->int_fact_phi[i][j] *= abs(det_J[0]);
			}
		}

		this->int_fact_dphi.resize(this->master.int_fact_dphi.size());
		for (int i = 0; i < this->master.int_fact_dphi.size(); i++) {
			this->int_fact_dphi[i].resize(dimension);
			for (int j = 0; j < dimension; j++) {
				this->int_fact_dphi[i][j].reserve(this->master.int_fact_dphi[i][j].size());
				for (int k = 0; k < this->master.int_fact_dphi[i][j].size(); k++) {
					double int_dphi = 0;
					for (int l = 0; l < dimension; l++) {
						int_dphi += this->master.int_fact_dphi[i][l][k] * J_inv[l][j][0];
					}
					int_dphi *= abs(det_J[0]);
					this->int_fact_dphi[i][j].push_back(int_dphi);
				}
			}
		}

		//MASS MATRIX
		this->m_inv = this->master.m_inv;
		for (int i = 0; i < this->m_inv.second.size(); i++) {
			for (int j = 0; j < this->m_inv.second[i].size(); j++) {
				this->m_inv.second[i][j] /= abs(det_J[0]);
			}
		}
	}
	else {
		//Placeholder for nonconstant Jacobian
	}

	this->data.state.push_back(SWE::State(this->int_fact_phi.size()));
	this->data.internal = SWE::Internal((*this->int_fact_phi.begin()).size());
}

template<int dimension, class master_type, class shape_type>
std::vector<RawBoundary<dimension - 1>*> Element<dimension, master_type, shape_type>::CreateBoundaries() {
	std::vector<RawBoundary<dimension - 1>*> my_raw_boundaries;

	Basis::Basis<dimension>* basis;

	std::function<std::vector<Point<dimension>>(const std::vector<Point<dimension - 1>>&)> boundary_to_master;
	std::function<Array2D<double>()> get_surface_normal;
	std::function<std::vector<double>()> get_surface_J;

	for (int i = 0; i < this->boundary_type.size(); i++) {
		basis = (Basis::Basis<dimension>*)(&master.basis);

		boundary_to_master = std::bind(&master_type::BoundaryToMasterCoordinates,
			master, i, std::placeholders::_1);

		get_surface_normal = std::bind(&shape_type::get_surface_normal,
			shape, i, nodal_coordinates);

		get_surface_J = std::bind(&shape_type::get_surface_J,
			shape, i, nodal_coordinates);

		my_raw_boundaries.push_back(new RawBoundary<dimension - 1>(this->boundary_type[i], this->neighbor_ID[i], master.p,
			this->data, basis, boundary_to_master, get_surface_normal, get_surface_J));
	}

	return my_raw_boundaries;
}

template<int dimension, class master_type, class shape_type>
void Element<dimension, master_type, shape_type>::ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp) {
	for (int i = 0; i < u_gp.size(); i++) {
		u_gp[i] = 0.0;
	}

	for (int i = 0; i < u.size(); i++) {
		for (int j = 0; j < u_gp.size(); j++) {
			u_gp[j] += u[i] * this->master.phi_gp[i][j];
		}
	}
}

template<int dimension, class master_type, class shape_type>
void Element<dimension, master_type, shape_type>::ComputeDUgp(int dir, const std::vector<double>& u, std::vector<double>& du_gp) {
	for (int i = 0; i < this->u_gp.size(); i++) {
		this->du_gp[i] = 0.0;
	}

	for (int i = 0; i < this->u.size(); i++) {
		for (int j = 0; j < this->du_gp.size(); j++) {
			this->du_gp[j] += this->u[i] * this->internal_dphi_fact[i][dir][j];
		}
	}
}


template<int dimension, class master_type, class shape_type>
double Element<dimension, master_type, shape_type>::IntegrationPhi(int phi_n, const std::vector<double>& u_gp) {
	double integral = 0;

	for (int i = 0; i < this->int_fact_phi[phi_n].size(); i++) {
		integral += u_gp[i] * this->int_fact_phi[phi_n][i];
	}
	
	return integral;
}

template<int dimension, class master_type, class shape_type>
double Element<dimension, master_type, shape_type>::IntegrationDPhi(int dir, int phi_n, const std::vector<double>& u_gp) {
	double integral = 0;

	for (int i = 0; i < this->int_fact_dphi[phi_n][dir].size(); i++) {
		integral += u_gp[i] * this->int_fact_dphi[phi_n][dir][i];
	}

	return integral;
}

template<int dimension, class master_type, class shape_type>
std::vector<double> Element<dimension, master_type, shape_type>::SolveLSE(const std::vector<double>& rhs) {
	std::vector<double> solution;

	if (this->m_inv.first) { //diagonal
		for (int i = 0; i < rhs.size(); i++) {
			solution.push_back(this->m_inv.second[0][i] * rhs[i]);
		}
	}
	else if (!(this->m_inv.first)) { //not diagonal
		for (int i = 0; i < this->m_inv.second.size(); i++) {
			solution.push_back(0);
			for (int j = 0; j < rhs.size(); j++) {
				solution[i] += this->m_inv.second[i][j] * rhs[j];
			}
		}
	}

	return solution;
}

template<int dimension, class master_type, class shape_type>
void Element<dimension, master_type, shape_type>::InitializeVTK(std::vector<Point<3>>& points, Array2D<unsigned int>& cells) {
	this->shape.get_VTK(points, cells, nodal_coordinates);
}

template<int dimension, class master_type, class shape_type>
void Element<dimension, master_type, shape_type>::WriteCellDataVTK(std::vector<double>& cell_data, int u_flag) {
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

template<int dimension, class master_type, class shape_type>
void Element<dimension, master_type, shape_type>::WritePointDataVTK(std::vector<double>& point_data, int u_flag) {
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