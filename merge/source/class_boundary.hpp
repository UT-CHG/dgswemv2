#ifndef CLASS_BOUNDARY_HPP
#define CLASS_BOUNDARY_HPP

#include "general_definitions.hpp"

#include "integration/integrations_1D.hpp"

template<uint dimension = 1, class data_type = SWE::Data>
class RawBoundary {
public:
	unsigned char type;
	uint neighbor_ID;

	uint p;

	Basis::Basis<dimension + 1>* basis;
	std::function<std::vector<Point<dimension + 1>>(const std::vector<Point<dimension>>&)> boundary_to_master;
	std::function<Array2D<double>()> get_surface_normal;
	std::function<std::vector<double>()> get_surface_J;

	data_type& data;

	RawBoundary(unsigned char type, uint neighbor_ID, uint p, data_type& data, Basis::Basis<dimension + 1>* basis,
		std::function<std::vector<Point<dimension + 1>>(const std::vector<Point<dimension>>&)>& boundary_to_master,
		std::function<Array2D<double>()>& get_surface_normal, std::function<std::vector<double>()>& get_surface_J) :

		type(type), neighbor_ID(neighbor_ID), p(p), data(data), basis(basis), boundary_to_master(boundary_to_master),
		get_surface_normal(get_surface_normal), get_surface_J(get_surface_J) {};
};

template<uint dimension = 1, class integration_type = Integration::GaussLegendre_1D, class data_type = SWE::Data>
class Boundary {
public:
	data_type& data;
	Array2D<double> surface_normal;

private:
	Array2D<double> phi_gp;
	Array2D<double> int_fact_phi;

public:
	Boundary(RawBoundary<dimension>&);

	void ComputeUgp(const std::vector<double>&, std::vector<double>&);
	double IntegrationPhi(uint, const std::vector<double>&);
};

template<uint dimension, class integration_type, class data_type>
Boundary<dimension, integration_type, data_type>::Boundary(RawBoundary<dimension>& raw_boundary) : data(raw_boundary.data) {
	integration_type integration;
	std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.get_rule(2 * raw_boundary.p);

	std::vector<Point<dimension + 1>> z_master = raw_boundary.boundary_to_master(integration_rule.second);

	this->phi_gp = raw_boundary.basis->get_phi(raw_boundary.p, z_master);

	std::vector<double> surface_J = raw_boundary.get_surface_J();

	if (surface_J.size() == 1) { //constant Jacobian
		this->int_fact_phi = this->phi_gp;
		for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {
			for (uint gp = 0; gp < this->int_fact_phi[dof].size(); gp++) {
				this->int_fact_phi[dof][gp] *= integration_rule.first[gp] * surface_J[0];
			}
		}

		this->surface_normal = Array2D<double>(integration_rule.first.size(), *raw_boundary.get_surface_normal().begin());
	}

	this->data.set_ngp_boundary(integration_rule.first.size());
}

template<uint dimension, class integration_type, class data_type>
void Boundary<dimension, integration_type, data_type>::ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp) {
	std::fill(u_gp.begin(), u_gp.end(), 0.0);

	for (uint dof = 0; dof < u.size(); dof++) {
		for (uint gp = 0; gp < u_gp.size(); gp++) {
			u_gp[gp] += u[dof] * this->phi_gp[dof][gp];
		}
	}
}

template<uint dimension, class integration_type, class data_type>
double Boundary<dimension, integration_type, data_type>::IntegrationPhi(uint phi_n, const std::vector<double>& u_gp) {
	double integral = 0;

	for (uint gp = 0; gp < this->int_fact_phi[phi_n].size(); gp++) {
		integral += u_gp[gp] * this->int_fact_phi[phi_n][gp];
	}

	return integral;
}

template<uint dimension = 1, class integration_type = Integration::GaussLegendre_1D, class data_type = SWE::Data>
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
	Interface(RawBoundary<dimension>&, RawBoundary<dimension>&);

	void ComputeUgpIN(const std::vector<double>&, std::vector<double>&);
	void ComputeUgpEX(const std::vector<double>&, std::vector<double>&);
	double IntegrationPhiIN(uint, const std::vector<double>&);
	double IntegrationPhiEX(uint, const std::vector<double>&);
};

template<uint dimension, class integration_type, class data_type>
Interface<dimension, integration_type, data_type>::Interface(RawBoundary<dimension>& raw_boundary_in, RawBoundary<dimension>& raw_boundary_ex) : 
	data_in(raw_boundary_in.data), data_ex(raw_boundary_ex.data) 
{
	uint p = std::max(raw_boundary_in.p, raw_boundary_ex.p);

	integration_type integration;
	std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.get_rule(2 * p);

	std::vector<Point<dimension + 1>> z_master = raw_boundary_in.boundary_to_master(integration_rule.second);
	this->phi_gp_in = raw_boundary_in.basis->get_phi(raw_boundary_in.p, z_master);

	z_master = raw_boundary_ex.boundary_to_master(integration_rule.second);
	this->phi_gp_ex = raw_boundary_ex.basis->get_phi(raw_boundary_ex.p, z_master);

	std::vector<double> surface_J = raw_boundary_in.get_surface_J();
	
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
	
		this->surface_normal = Array2D<double>(integration_rule.first.size(), *raw_boundary_in.get_surface_normal().begin());
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

	for (uint gp = 0; gp < this->int_fact_phi_in[phi_n].size(); gp++) {
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

	for (uint gp = 0; gp < this->int_fact_phi_ex[phi_n].size(); gp++) {
		integral += u_gp[gp] * this->int_fact_phi_ex[phi_n][gp];
	}

	return integral;
}

#endif