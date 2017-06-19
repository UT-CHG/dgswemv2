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

	Array2D<double>& u;

	RawBoundary(unsigned char type, unsigned int neighbor_ID, int p, Array2D<double>& u, Basis::Basis<dimension + 1>* basis,
		std::function<std::vector<Point<dimension + 1>>(const std::vector<Point<dimension>>&)>& boundary_to_master,
		std::function<Array2D<double>()>& get_surface_normal, std::function<std::vector<double>()>& get_surface_J) :
		type(type), neighbor_ID(neighbor_ID), p(p), u(u), basis(basis), boundary_to_master(boundary_to_master),
		get_surface_normal(get_surface_normal), get_surface_J(get_surface_J) {};
};

template<int dimension = 1, class integration_type = Integration::GaussLegendre_1D>
class Boundary {
private:
	Array2D<double> phi_gp;
	Array2D<double> int_fact_phi;
	
	Array2D<double> surface_normal;

	Array2D<double>& u;
	Array2D<double> u_gp;
public:
	Boundary(RawBoundary<dimension>& raw_boundary);

	void ComputeUgp(int);
	double IntegrationPhi(int, int);
};

template<int dimension, class integration_type>
Boundary<dimension, integration_type>::Boundary(RawBoundary<dimension>& raw_boundary) : u(raw_boundary.u) {
	integration_type integration;
	std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.get_rule(2 * raw_boundary.p);

	std::vector<Point<dimension + 1>> z_master = raw_boundary.boundary_to_master(integration_rule.second);
	this->phi_gp = raw_boundary.basis->get_phi(raw_boundary.p, z_master);

	this->surface_normal = raw_boundary.get_surface_normal();
	printf("%f %f\n", surface_normal[0][X], surface_normal[0][Y]);

	std::vector<double> surface_J = raw_boundary.get_surface_J();
	printf("%f\n", surface_J[0]);

	if (surface_J.size() == 1) { //constant Jacobian
		this->int_fact_phi = this->phi_gp;
		for (int i = 0; i < this->int_fact_phi.size(); i++) {
			for (int j = 0; j < this->int_fact_phi[i].size(); j++) {
				this->int_fact_phi[i][j] *= integration_rule.first[j] * surface_J[0];
			}
		}
	}

	this->u_gp.resize(SIZE_U_BOUNDARY);
	for (int i = 0; i < SIZE_U_BOUNDARY; i++) {
		this->u_gp[i].resize(integration_rule.first.size());
	}
}

template<int dimension, class integration_type>
void Boundary<dimension, integration_type>::ComputeUgp(int u_flag) {
	for (int i = 0; i < this->u_gp[u_flag].size(); i++) {
		this->u_gp[u_flag][i] = 0.0;
	}

	for (int i = 0; i < this->u[u_flag].size(); i++) {
		for (int j = 0; j < this->u_gp[u_flag].size(); j++) {
			this->u_gp[u_flag][j] += this->u[u_flag][i] * this->phi_gp[i][j];
		}
	}
}

template<int dimension, class integration_type>
double Boundary<dimension, integration_type>::IntegrationPhi(int u_flag, int phi_n) {
	double integral = 0;

	for (int j = 0; j < this->int_fac_phi[phi_n].size(); j++) {
		integral += this->u_gp[u_flag][j] * this->int_fac_ph[phi_n][j];
	}
	
	return integral;
}

#endif