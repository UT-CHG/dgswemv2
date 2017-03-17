#include "class_interface_2D.h"

INTERFACE_2D::INTERFACE_2D(int number_gp, double** u_boundary, double* normal_x, double* normal_y, bool straight)
	: INTERFACE(number_gp, u_boundary) 
{
	this->normal_x = new double[this->number_gp];
	this->normal_y = new double[this->number_gp];

	if (straight) {
		for (int i = 0; i < this->number_gp; i++) {
			this->normal_x[i] = normal_x[0];
			this->normal_y[i] = normal_y[0];
		}
	}
	else if (!straight) {
		for (int i = 0; i < this->number_gp; i++) {
			this->normal_x[i] = normal_x[i];
			this->normal_y[i] = normal_y[i];
		}
	}
}

INTERFACE_2D::~INTERFACE_2D() {
	delete this->normal_x;
	delete this->normal_y;
}