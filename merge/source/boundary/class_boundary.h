#ifndef CLASS_BOUNDARY_H
#define CLASS_BOUNDARY_H

#include "../general_definitions.h"

#include "../integration/integrations_1D.h"

template<int dimension = 1>
class RawBoundary {
public:
	unsigned char type;
	unsigned int neighbor_ID;

	int p;

	Basis::Basis<dimension + 1>* basis;
	std::function<std::vector<Point<dimension + 1>>(const std::vector<Point<dimension>>&)> boundary_to_master;
	std::function<Array2D<double>()> get_surface_normal;
	std::function<std::vector<double>()> get_surface_J;

	SWE::Data& data;

	RawBoundary(unsigned char type, unsigned int neighbor_ID, int p, SWE::Data& data, Basis::Basis<dimension + 1>* basis,
		std::function<std::vector<Point<dimension + 1>>(const std::vector<Point<dimension>>&)>& boundary_to_master,
		std::function<Array2D<double>()>& get_surface_normal, std::function<std::vector<double>()>& get_surface_J) :
		
		type(type), neighbor_ID(neighbor_ID), p(p), data(data), basis(basis), boundary_to_master(boundary_to_master),
		get_surface_normal(get_surface_normal), get_surface_J(get_surface_J) {};
};

template<int dimension = 1, class integration_type = Integration::GaussLegendre_1D>
class Boundary {
private:
	Array2D<double> phi_gp;
	Array2D<double> int_fact_phi;
	
	Array2D<double> surface_normal;

	SWE::Data& data;
public:
	Boundary(RawBoundary<dimension>&);

	void ComputeUgp(std::vector<double>&, const std::vector<double>&);
	double IntegrationPhi(int, const std::vector<double>&);
};

template<int dimension, class integration_type>
Boundary<dimension, integration_type>::Boundary(RawBoundary<dimension>& raw_boundary) : data(raw_boundary.data) {
	integration_type integration;
	std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.get_rule(2 * raw_boundary.p);

	std::vector<Point<dimension + 1>> z_master = raw_boundary.boundary_to_master(integration_rule.second);
	this->phi_gp = raw_boundary.basis->get_phi(raw_boundary.p, z_master);

	std::vector<double> surface_J = raw_boundary.get_surface_J();

	if (surface_J.size() == 1) { //constant Jacobian
		this->int_fact_phi = this->phi_gp;
		for (int i = 0; i < this->int_fact_phi.size(); i++) {
			for (int j = 0; j < this->int_fact_phi[i].size(); j++) {
				this->int_fact_phi[i][j] *= integration_rule.first[j] * surface_J[0];
			}
		}
	}

	this->surface_normal = raw_boundary.get_surface_normal();

	this->data.boundary_data = SWE::Boundary((*this->int_fact_phi.begin()).size());
}

template<int dimension, class integration_type>
void Boundary<dimension, integration_type>
::ComputeUgp(std::vector<double>& u_gp, const std::vector<double>& u) {
	for (int i = 0; i < this->u_gp.size(); i++) {
		this->u_gp[i] = 0.0;
	}

	for (int i = 0; i < this->u.size(); i++) {
		for (int j = 0; j < this->u_gp.size(); j++) {
			this->u_gp[j] += this->u[i] * this->phi_gp[i][j];
		}
	}
}

template<int dimension, class integration_type>
double Boundary<dimension, integration_type>
::IntegrationPhi(int phi_n, const std::vector<double>& u_gp) {
	double integral = 0;

	for (int i = 0; i < this->int_fac_phi[phi_n].size(); i++) {
		integral += u_gp[i] * this->int_fac_phi[phi_n][i];
	}
	
	return integral;
}

template<int dimension = 1, class integration_type = Integration::GaussLegendre_1D>
class Interface {
private:
	Array2D<double> phi_gp_in;
	Array2D<double> phi_gp_ex;
	Array2D<double> int_fact_phi_in;
	Array2D<double> int_fact_phi_ex;
	
	Array2D<double> surface_normal;

	SWE::Data& data_in;
	SWE::Data& data_ex;
public:
	Interface(RawBoundary<dimension>&, RawBoundary<dimension>&);

	void ComputeUgpIN(std::vector<double>&, const std::vector<double>&);
	void ComputeUgpEX(std::vector<double>&, const std::vector<double>&);
	double IntegrationPhiIN(int, const std::vector<double>&);
	double IntegrationPhiEX(int, const std::vector<double>&);
};

template<int dimension, class integration_type>
Interface<dimension, integration_type>
::Interface(RawBoundary<dimension>& raw_boundary_in, RawBoundary<dimension>& raw_boundary_ex) : 
data_in(raw_boundary_in.data), data_ex(raw_boundary_ex.data) {
	int p = std::max(raw_boundary_in.p,raw_boundary_ex.p);

	integration_type integration;
	std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.get_rule(2 * p);

	std::vector<Point<dimension + 1>> z_master = raw_boundary_in.boundary_to_master(integration_rule.second);
	this->phi_gp_in = raw_boundary_in.basis->get_phi(p, z_master);

	z_master = raw_boundary_ex.boundary_to_master(integration_rule.second);
	this->phi_gp_ex = raw_boundary_ex.basis->get_phi(p, z_master);
	
	std::vector<double> surface_J = raw_boundary_in.get_surface_J();

	if (surface_J.size() == 1) { //constant Jacobian
		this->int_fact_phi_in = this->phi_gp_in;
		this->int_fact_phi_ex = this->phi_gp_ex;
		for (int i = 0; i < this->int_fact_phi_in.size(); i++) {
			for (int j = 0; j < this->int_fact_phi_in[i].size(); j++) {
				this->int_fact_phi_in[i][j] *= integration_rule.first[j] * surface_J[0];
				this->int_fact_phi_ex[i][j] *= integration_rule.first[j] * surface_J[0];
			}
		}
	}

	this->surface_normal = raw_boundary_in.get_surface_normal();

	this->data_in.boundary_data = SWE::Boundary((*this->int_fact_phi_in.begin()).size());
	this->data_ex.boundary_data = SWE::Boundary((*this->int_fact_phi_ex.begin()).size());
}

template<int dimension, class integration_type>
void Interface<dimension, integration_type>
::ComputeUgpIN(std::vector<double>& u_gp, const std::vector<double>& u) {
	for (int i = 0; i < this->u_gp.size(); i++) {
		this->u_gp[i] = 0.0;
	}

	for (int i = 0; i < this->u.size(); i++) {
		for (int j = 0; j < this->u_gp.size(); j++) {
			this->u_gp[j] += this->u[i] * this->phi_gp_in[i][j];
		}
	}
}

template<int dimension, class integration_type>
double Interface<dimension, integration_type>
::IntegrationPhiIN(int phi_n, const std::vector<double>& u_gp) {
	double integral = 0;

	for (int i = 0; i < this->int_fac_phi_in[phi_n].size(); i++) {
		integral += u_gp[i] * this->int_fac_phi_in[phi_n][i];
	}
	
	return integral;
}

template<int dimension, class integration_type>
void Interface<dimension, integration_type>
::ComputeUgpEX(std::vector<double>& u_gp, const std::vector<double>& u) {
	for (int i = 0; i < this->u_gp.size(); i++) {
		this->u_gp[i] = 0.0;
	}

	for (int i = 0; i < this->u.size(); i++) {
		for (int j = 0; j < this->u_gp.size(); j++) {
			this->u_gp[j] += this->u[i] * this->phi_gp_ex[i][j];
		}
	}
}

template<int dimension, class integration_type>
double Interface<dimension, integration_type>
::IntegrationPhiEX(int phi_n, const std::vector<double>& u_gp) {
	double integral = 0;

	for (int i = 0; i < this->int_fac_phi_ex[phi_n].size(); i++) {
		integral += u_gp[i] * this->int_fac_phi_ex[phi_n][i];
	}
	
	return integral;
}

#endif