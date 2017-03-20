#include "class_mesh_v2.h"

MESH::MESH(int p, int p_geom) {
	this->p = p;
	this->p_geom = p_geom;

	this->InitializeElements();
	this->InitializeInterfaces();
}

MESH::~MESH() {
	for (auto it = this->interfaces.begin(); it != this->interfaces.end(); it++) {
		for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
			delete *itt;
		}
		it->second.clear();
	}
	this->interfaces.clear();

	for (auto it = this->elements.begin(); it != this->elements.end(); it++) {
		delete it->second;
	}
	this->elements.clear();

	for (auto it = this->bases_2D.begin(); it != this->bases_2D.end(); it++) {
		delete it->second;
	}
	this->bases_2D.clear();

	for (auto it = this->geometric_bases_2D.begin(); it != this->geometric_bases_2D.end(); it++) {
		delete it->second;
	}
	this->geometric_bases_2D.clear();

	for (auto it = this->line_rules.begin(); it != this->line_rules.end(); it++) {
		delete it->second;
	}
	this->line_rules.clear();

	for (auto it = this->area_rules.begin(); it != this->area_rules.end(); it++) {
		delete it->second;
	}
	this->area_rules.clear();
}

void MESH::InitializeElements() {
	this->line_rules[GAUSS_LEGENDRE] = new INTEGRATION_1D(GAUSS_LEGENDRE, 2 * this->p);
	this->area_rules[DUNAVANT] = new INTEGRATION_2D(DUNAVANT, 2 * this->p);

	INTEGRATION_1D* line_rule;

	if (this->line_rules.find(GAUSS_LEGENDRE)->first == GAUSS_LEGENDRE) {
		line_rule = this->line_rules.find(GAUSS_LEGENDRE)->second;
	}
	else {
		this->line_rules[GAUSS_LEGENDRE] = new INTEGRATION_1D(GAUSS_LEGENDRE, 2 * p);
		line_rule = this->line_rules.find(GAUSS_LEGENDRE)->second;
	}

	INTEGRATION_2D* area_rule;

	if (this->area_rules.find(DUNAVANT)->first == DUNAVANT) {
		area_rule = this->area_rules.find(DUNAVANT)->second;
	}
	else {
		this->area_rules[DUNAVANT] = new INTEGRATION_2D(DUNAVANT, 2 * p);
		area_rule = this->area_rules.find(DUNAVANT)->second;
	}

	this->bases_2D[DUBINER] = new BASIS_2D(DUBINER, p, line_rule, area_rule);

	BASIS_2D* basis;

	if (this->bases_2D.find(DUBINER)->first == DUBINER) {
		basis = this->bases_2D.find(DUBINER)->second;
	}
	else {
		this->bases_2D[DUBINER] = new BASIS_2D(DUBINER, p, line_rule, area_rule);
		basis = this->bases_2D.find(DUBINER)->second;
	}

	unsigned int ID = 0;

	double x[3] = { -1,  1, -1 };
	double y[3] = { -1, -1,  1 };

	unsigned int neighbors[3] = { 1, DEFAULT_ID, DEFAULT_ID };
	unsigned char boundaries[3] = { INTERNAL, LAND, LAND };

	this->elements[ID] = new ELEMENT_2D(ID, neighbors, boundaries, x, y, basis);

	ID = 1;
	x[0] = 1;
	y[0] = 1;
	neighbors[0] = 0;

	this->elements[ID] = new ELEMENT_2D(ID, neighbors, boundaries, x, y, basis);
}

void MESH::InitializeInterfaces() {
	this->interfaces[INTERNAL] = std::vector<INTERFACE*>();
	this->interfaces[LAND] = std::vector<INTERFACE*>();

	for (auto it = this->elements.begin(); it != this->elements.end(); it++) {
		std::map<unsigned int, INTERFACE*> internal_interfaces(it->second->CreateInterfaces());

		for (auto itt = internal_interfaces.begin(); itt != internal_interfaces.end(); itt++) {
			this->elements.find(itt->first)->second->AppendInterface(it->first, itt->second);
		}

		std::vector<std::pair<unsigned char, INTERFACE*>> own_interfaces(it->second->GetOwnInterfaces());

		for (auto itt = own_interfaces.begin(); itt != own_interfaces.end(); itt++) {
			this->interfaces.find(itt->first)->second.push_back(itt->second);
		}
	}
}

//void MESH::solve() {
//}