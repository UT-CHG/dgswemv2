#include "class_interface.h"

INTERFACE::INTERFACE(unsigned char dimension, int number_gp, double** u_boundary, double** normal, bool boundary) {
	this->dimension = dimension;
	this->boundary = boundary;
	this->number_gp = number_gp;
	
	this->normal = normal;
	
	this->u_boundary_in = u_boundary;

    if (this->boundary) {
		this->u_boundary_ex = new double*[SIZE_U_BOUNDARY];
        for (int i = 0; i < SIZE_U_BOUNDARY; i++) {
            this->u_boundary_ex[i] = new double[this->number_gp];
        }
    }
}

INTERFACE::~INTERFACE() {
	if (this->boundary) {
		for (int i = 0; i < SIZE_U_BOUNDARY; i++) {
			delete[] this->u_boundary_ex[i];
		}
		delete[] this->u_boundary_ex;
	}

	for (int i = 0; i < this->dimension; i++) {
		delete[] this->normal[i];
	}
	delete[] this->normal;
}