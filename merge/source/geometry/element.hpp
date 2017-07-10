#ifndef CLASS_ELEMENT_HPP
#define CLASS_ELEMENT_HPP

#include "boundary.hpp"

namespace Geometry {
	template<uint dimension, class master_type, class shape_type, class data_type>
	class Element {
	public:
		data_type data;

	private:
		uint ID;

		master_type& master;
		shape_type shape;

		std::vector<uint> neighbor_ID;
		std::vector<unsigned char> boundary_type;

		Array3D<double> dphi_fact;
		Array2D<double> int_fact_phi;
		Array3D<double> int_fact_dphi;

		std::pair<bool, Array2D<double>> m_inv;

	public:
		Element(uint, master_type&, std::vector<Point<dimension>>&,
			std::vector<uint>&, std::vector<unsigned char>&);

		void CreateRawBoundaries(std::map<uint, std::map<uint, RawBoundary<dimension - 1, data_type>>>&,
			std::map<unsigned char, std::vector<RawBoundary<dimension - 1, data_type>>>&);

		void ComputeUgp(const std::vector<double>&, std::vector<double>&);
		void ComputeDUgp(uint, const std::vector<double>&, std::vector<double>&);
		double IntegrationPhi(uint, const std::vector<double>&);
		double IntegrationDPhi(uint, uint, const std::vector<double>&);
	
		std::vector<double> SolveLSE(const std::vector<double>&);

		void InitializeVTK(std::vector<Point<3>>&, Array2D<uint>&);
		void WriteCellDataVTK(const std::vector<double>&, std::vector<double>&);
		void WritePointDataVTK(const std::vector<double>&, std::vector<double>&);

		typedef master_type master_element_type;
	};

	template<uint dimension, class master_type, class shape_type, class data_type>
	Element<dimension, master_type, shape_type, data_type>::Element(uint ID, master_type& master, std::vector<Point<dimension>>& nodal_coordinates,
		std::vector<uint>& neighbor_ID, std::vector<unsigned char>& boundary_type) :
		ID(ID), master(master), shape(shape_type(nodal_coordinates)),
		neighbor_ID(neighbor_ID), boundary_type(boundary_type)
	{
		//DEFORMATION
		std::vector<Point<dimension>> gp = this->master.integration.GetRule(2 * this->master.p).second;

		std::vector<double> det_J = this->shape.GetJdet(gp);
		Array3D<double> J_inv = this->shape.GetJinv(gp);

		if (det_J.size() == 1) { //constant Jacobian
			//DIFFERENTIATION FACTORS
			this->dphi_fact.resize(this->master.dphi_gp.size());
			for (uint dof = 0; dof < this->master.dphi_gp.size(); dof++) {
				this->dphi_fact[dof].resize(dimension);
				for (uint dir = 0; dir < dimension; dir++) {
					this->dphi_fact[dof][dir].reserve(this->master.dphi_gp[dof][dir].size());
					for (uint gp = 0; gp < this->master.dphi_gp[dof][dir].size(); gp++) {
						double dphi = 0;
						for (uint z = 0; z < dimension; z++) {
							dphi += this->master.dphi_gp[dof][z][gp] * J_inv[z][dir][0];
						}
						this->dphi_fact[dof][dir].push_back(dphi);
					}
				}
			}

			//INTEGRATION FACTORS
			this->int_fact_phi = this->master.int_fact_phi;
			for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {
				for (uint gp = 0; gp < this->int_fact_phi[dof].size(); gp++) {
					this->int_fact_phi[dof][gp] *= std::abs(det_J[0]);
				}
			}

			this->int_fact_dphi.resize(this->master.int_fact_dphi.size());
			for (uint dof = 0; dof < this->master.int_fact_dphi.size(); dof++) {
				this->int_fact_dphi[dof].resize(dimension);
				for (uint dir = 0; dir < dimension; dir++) {
					this->int_fact_dphi[dof][dir].reserve(this->master.int_fact_dphi[dof][dir].size());
					for (uint gp = 0; gp < this->master.int_fact_dphi[dof][dir].size(); gp++) {
						double int_dphi = 0;
						for (uint z = 0; z < dimension; z++) {
							int_dphi += this->master.int_fact_dphi[dof][z][gp] * J_inv[z][dir][0];
						}
						int_dphi *= std::abs(det_J[0]);
						this->int_fact_dphi[dof][dir].push_back(int_dphi);
					}
				}
			}

			//MASS MATRIX
			this->m_inv = this->master.m_inv;
			for (uint i = 0; i < this->m_inv.second.size(); i++) {
				for (uint j = 0; j < this->m_inv.second[i].size(); j++) {
					this->m_inv.second[i][j] /= std::abs(det_J[0]);
				}
			}
		}
		else {
			//Placeholder for nonconstant Jacobian
		}

		this->data.set_ndof(this->master.phi_gp.size());
		this->data.set_ngp_internal((*this->master.phi_gp.begin()).size());
	}

	template<uint dimension, class master_type, class shape_type, class data_type>
	void Element<dimension, master_type, shape_type, data_type>::CreateRawBoundaries(
		std::map<uint, std::map<uint, RawBoundary<dimension - 1, data_type>>>& pre_interfaces,
		std::map<unsigned char, std::vector<RawBoundary<dimension - 1, data_type>>>& pre_boundaries)
	{
		Basis::Basis<dimension>* my_basis = (Basis::Basis<dimension>*)(&this->master.basis);
		Master::Master<dimension>* my_master = (Master::Master<dimension>*)(&this->master);
		Shape::Shape<dimension>* my_shape = (Shape::Shape<dimension>*)(&this->shape);

		for (uint i = 0; i < this->boundary_type.size(); i++) {
			if (this->boundary_type[i] == INTERNAL) {
				pre_interfaces[this->ID].emplace(
					std::make_pair(this->neighbor_ID[i],
						RawBoundary<dimension - 1, data_type>(this->master.p, i, this->data, *my_basis, *my_master, *my_shape)));
			}
			else {
				pre_boundaries[this->boundary_type[i]].emplace_back(
					RawBoundary<dimension - 1, data_type>(this->master.p, i, this->data, *my_basis, *my_master, *my_shape));
			}
		}
	}

	template<uint dimension, class master_type, class shape_type, class data_type>
	void Element<dimension, master_type, shape_type, data_type>::ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp) {
		std::fill(u_gp.begin(), u_gp.end(), 0.0);

		for (uint dof = 0; dof < u.size(); dof++) {
			for (uint gp = 0; gp < u_gp.size(); gp++) {
				u_gp[gp] += u[dof] * this->master.phi_gp[dof][gp];
			}
		}
	}

	template<uint dimension, class master_type, class shape_type, class data_type>
	void Element<dimension, master_type, shape_type, data_type>::ComputeDUgp(uint dir, const std::vector<double>& u, std::vector<double>& du_gp) {
		std::fill(du_gp.begin(), du_gp.end(), 0.0);

		for (uint dof = 0; dof < u.size(); dof++) {
			for (uint gp = 0; gp < du_gp.size(); gp++) {
				du_gp[gp] += u[dof] * this->dphi_fact[dof][dir][gp];
			}
		}
	}


	template<uint dimension, class master_type, class shape_type, class data_type>
	double Element<dimension, master_type, shape_type, data_type>::IntegrationPhi(uint phi_n, const std::vector<double>& u_gp) {
		double integral = 0;

		for (uint gp = 0; gp < u_gp.size(); gp++) {
			integral += u_gp[gp] * this->int_fact_phi[phi_n][gp];
		}

		return integral;
	}

	template<uint dimension, class master_type, class shape_type, class data_type>
	double Element<dimension, master_type, shape_type, data_type>::IntegrationDPhi(uint dir, uint phi_n, const std::vector<double>& u_gp) {
		double integral = 0;

		for (uint gp = 0; gp < u_gp.size(); gp++) {
			integral += u_gp[gp] * this->int_fact_dphi[phi_n][dir][gp];
		}

		return integral;
	}

	template<uint dimension, class master_type, class shape_type, class data_type>
	std::vector<double> Element<dimension, master_type, shape_type, data_type>::SolveLSE(const std::vector<double>& rhs) {
		std::vector<double> solution;

		if (this->m_inv.first) { //diagonal
			for (uint i = 0; i < rhs.size(); i++) {
				solution.push_back(this->m_inv.second[0][i] * rhs[i]);
			}
		}
		else if (!(this->m_inv.first)) { //not diagonal
			for (uint i = 0; i < this->m_inv.second.size(); i++) {
				solution.push_back(0);
				for (uint j = 0; j < rhs.size(); j++) {
					solution[i] += this->m_inv.second[i][j] * rhs[j];
				}
			}
		}

		return solution;
	}

	template<uint dimension, class master_type, class shape_type, class data_type>
	void Element<dimension, master_type, shape_type, data_type>::InitializeVTK(std::vector<Point<3>>& points, Array2D<uint>& cells) {
		this->shape.GetVTK(points, cells);
	}

	template<uint dimension, class master_type, class shape_type, class data_type>
	void Element<dimension, master_type, shape_type, data_type>::WriteCellDataVTK(const std::vector<double>& u, std::vector<double>& cell_data) {
		Array2D<double> temp = this->master.phi_postprocessor_cell;

		for (uint dof = 0; dof < temp.size(); dof++) {
			for (uint cell = 0; cell < temp[dof].size(); cell++) {
				temp[dof][cell] *= u[dof];
			}
		}

		for (uint dof = 1; dof < temp.size(); dof++) {
			for (uint cell = 0; cell < temp[dof].size(); cell++) {
				temp[0][cell] += temp[dof][cell];
			}
		}

		cell_data.insert(cell_data.end(), temp[0].begin(), temp[0].end());
	}

	template<uint dimension, class master_type, class shape_type, class data_type>
	void Element<dimension, master_type, shape_type, data_type>::WritePointDataVTK(const std::vector<double>& u, std::vector<double>& point_data) {
		Array2D<double> temp = this->master.phi_postprocessor_point;

		for (uint dof = 0; dof < temp.size(); dof++) {
			for (uint pt = 0; pt < temp[dof].size(); pt++) {
				temp[dof][pt] *= u[dof];
			}
		}

		for (uint dof = 1; dof < temp.size(); dof++) {
			for (uint pt = 0; pt < temp[dof].size(); pt++) {
				temp[0][pt] += temp[dof][pt];
			}
		}

		point_data.insert(point_data.end(), temp[0].begin(), temp[0].end());
	}
}

#endif