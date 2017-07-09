#ifndef CLASS_BOUNDARY_HPP
#define CLASS_BOUNDARY_HPP

#include "../general_definitions.hpp"

namespace Geometry {
	template<uint dimension, class data_type>
	class RawBoundary {
	public:
		uint p;
		uint n_bound;

		data_type& data;

		Basis::Basis<dimension + 1>& basis;
		Master::Master<dimension + 1>& master;
		Shape::Shape<dimension + 1>& shape;

		RawBoundary(uint p, uint n_bound, data_type& data,
			Basis::Basis<dimension + 1>& basis, Master::Master<dimension + 1>& master, Shape::Shape<dimension + 1>& shape) :
			p(p), n_bound(n_bound), data(data),
			basis(basis), master(master), shape(shape) {}
	};

	template<uint dimension, class integration_type, class data_type, class boundary_type>
	class Boundary {
	public:
		data_type& data;
		boundary_type boundary_condition;

		Array2D<double> surface_normal;

	private:
		Array2D<double> phi_gp;
		Array2D<double> int_fact_phi;

	public:
		Boundary(const RawBoundary<dimension, data_type>&);

		void ComputeUgp(const std::vector<double>&, std::vector<double>&);
		double IntegrationPhi(uint, const std::vector<double>&);
	};

	template<uint dimension, class integration_type, class data_type, class boundary_type>
	Boundary<dimension, integration_type, data_type, boundary_type>::Boundary(const RawBoundary<dimension, data_type>& raw_boundary) : data(raw_boundary.data) {
		integration_type integration;
		std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * raw_boundary.p);

		std::vector<Point<dimension + 1>> z_master =
			raw_boundary.master.BoundaryToMasterCoordinates(raw_boundary.n_bound, integration_rule.second);

		this->phi_gp = raw_boundary.basis.GetPhi(raw_boundary.p, z_master);

		std::vector<double> surface_J =
			raw_boundary.shape.GetSurfaceJ(raw_boundary.n_bound, z_master);

		if (surface_J.size() == 1) { //constant Jacobian
			this->int_fact_phi = this->phi_gp;
			for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {
				for (uint gp = 0; gp < this->int_fact_phi[dof].size(); gp++) {
					this->int_fact_phi[dof][gp] *= integration_rule.first[gp] * surface_J[0];
				}
			}

			this->surface_normal = Array2D<double>(integration_rule.first.size(),
				*raw_boundary.shape.GetSurfaceNormal(raw_boundary.n_bound, z_master).begin());
		}

		this->data.set_ngp_boundary(integration_rule.first.size());
	}

	template<uint dimension, class integration_type, class data_type, class boundary_type>
	void Boundary<dimension, integration_type, data_type, boundary_type>::ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp) {
		std::fill(u_gp.begin(), u_gp.end(), 0.0);

		for (uint dof = 0; dof < u.size(); dof++) {
			for (uint gp = 0; gp < u_gp.size(); gp++) {
				u_gp[gp] += u[dof] * this->phi_gp[dof][gp];
			}
		}
	}

	template<uint dimension, class integration_type, class data_type, class boundary_type>
	double Boundary<dimension, integration_type, data_type, boundary_type>::IntegrationPhi(uint phi_n, const std::vector<double>& u_gp) {
		double integral = 0;

		for (uint gp = 0; gp < u_gp.size(); gp++) {
			integral += u_gp[gp] * this->int_fact_phi[phi_n][gp];
		}

		return integral;
	}

	template<uint dimension, class integration_type, class data_type>
	class Interface {
	public:
		data_type& data_in;
		data_type& data_ex;
		Array2D<double> surface_normal;

	private:
		Array2D<double> phi_gp_in;
		Array2D<double> phi_gp_ex;
		Array2D<double> int_fact_phi_in;
		Array2D<double> int_fact_phi_ex;

	public:
		Interface(const RawBoundary<dimension, data_type>&, const RawBoundary<dimension, data_type>&);

		void ComputeUgpIN(const std::vector<double>&, std::vector<double>&);
		void ComputeUgpEX(const std::vector<double>&, std::vector<double>&);
		double IntegrationPhiIN(uint, const std::vector<double>&);
		double IntegrationPhiEX(uint, const std::vector<double>&);
	};

	template<uint dimension, class integration_type, class data_type>
	Interface<dimension, integration_type, data_type>::Interface
		(const RawBoundary<dimension, data_type>& raw_boundary_in, const RawBoundary<dimension, data_type>& raw_boundary_ex) :
		data_in(raw_boundary_in.data), data_ex(raw_boundary_ex.data)
	{
		uint p = std::max(raw_boundary_in.p, raw_boundary_ex.p);

		integration_type integration;
		std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * p);

		std::vector<Point<dimension + 1>> z_master =
			raw_boundary_ex.master.BoundaryToMasterCoordinates(raw_boundary_ex.n_bound, integration_rule.second);
		this->phi_gp_ex = raw_boundary_ex.basis.GetPhi(raw_boundary_ex.p, z_master);

		z_master =
			raw_boundary_in.master.BoundaryToMasterCoordinates(raw_boundary_in.n_bound, integration_rule.second);
		this->phi_gp_in = raw_boundary_in.basis.GetPhi(raw_boundary_in.p, z_master);

		std::vector<double> surface_J =
			raw_boundary_in.shape.GetSurfaceJ(raw_boundary_in.n_bound, z_master);

		if (surface_J.size() == 1) { //constant Jacobian
			this->int_fact_phi_in = this->phi_gp_in;
			for (uint dof = 0; dof < this->int_fact_phi_in.size(); dof++) {
				for (uint gp = 0; gp < this->int_fact_phi_in[dof].size(); gp++) {
					this->int_fact_phi_in[dof][gp] *= integration_rule.first[gp] * surface_J[0];
				}
			}

			this->int_fact_phi_ex = this->phi_gp_ex;
			for (uint dof = 0; dof < this->int_fact_phi_ex.size(); dof++) {
				for (uint gp = 0; gp < this->int_fact_phi_ex[dof].size(); gp++) {
					this->int_fact_phi_ex[dof][gp] *= integration_rule.first[gp] * surface_J[0];
				}
			}

			this->surface_normal = Array2D<double>(integration_rule.first.size(),
				*raw_boundary_in.shape.GetSurfaceNormal(raw_boundary_in.n_bound, z_master).begin());
		}

		this->data_in.set_ngp_boundary(integration_rule.first.size());
		this->data_ex.set_ngp_boundary(integration_rule.first.size());
	}

	template<uint dimension, class integration_type, class data_type>
	void Interface<dimension, integration_type, data_type>::ComputeUgpIN(const std::vector<double>& u, std::vector<double>& u_gp) {
		std::fill(u_gp.begin(), u_gp.end(), 0.0);

		for (uint dof = 0; dof < u.size(); dof++) {
			for (uint gp = 0; gp < u_gp.size(); gp++) {
				u_gp[gp] += u[dof] * this->phi_gp_in[dof][gp];
			}
		}
	}

	template<uint dimension, class integration_type, class data_type>
	double Interface<dimension, integration_type, data_type>::IntegrationPhiIN(uint phi_n, const std::vector<double>& u_gp) {
		double integral = 0;

		for (uint gp = 0; gp < u_gp.size(); gp++) {
			integral += u_gp[gp] * this->int_fact_phi_in[phi_n][gp];
		}

		return integral;
	}

	template<uint dimension, class integration_type, class data_type>
	void Interface<dimension, integration_type, data_type>::ComputeUgpEX(const std::vector<double>& u, std::vector<double>& u_gp) {
		std::fill(u_gp.begin(), u_gp.end(), 0.0);

		for (uint dof = 0; dof < u.size(); dof++) {
			for (uint gp = 0; gp < u_gp.size(); gp++) {
				u_gp[gp] += u[dof] * this->phi_gp_ex[dof][gp];
			}
		}
	}

	template<uint dimension, class integration_type, class data_type>
	double Interface<dimension, integration_type, data_type>::IntegrationPhiEX(uint phi_n, const std::vector<double>& u_gp) {
		double integral = 0;

		for (uint gp = 0; gp < u_gp.size(); gp++) {
			integral += u_gp[gp] * this->int_fact_phi_ex[phi_n][gp];
		}

		return integral;
	}
}

#endif