#include "class_interface_2D.h"

INTERFACE_2D::INTERFACE_2D(double** U_edge) : INTERFACE() {
	this->U_edge_IN = U_edge;
}

void INTERFACE_2D::SetPointerEX(double** U_edge) {
	this->U_edge_EX = U_edge;
}