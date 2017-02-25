#include "class_element.h"
#include <iostream>
using namespace std;

ELEMENT_TRI::ELEMENT_TRI(int ID, int p, int p_geom, double nodal_coordinates_x[],
	double nodal_coordinates_y[], BASIS_TRI* basis, BASIS_GEOM_TRI* basis_geom) {
	this->p = p;
	this->p_geom = p_geom;

	this->basis = basis;
	this->basis_geom = basis_geom;
}

ELEMENT_TRI::~ELEMENT_TRI(){
}

void ELEMENT_TRI::compute_geometry() {}

void ELEMENT_TRI::compute_integration_factors() {}