#include "class_interface_2D.h"

INTERFACE_2D::INTERFACE_2D(int number_gp, double** U_edge) : INTERFACE() {
	this->number_gp;
	this->U_edge_IN = U_edge;
}

void INTERFACE_2D::SetPointerEX(double** U_edge) {
	this->U_edge_EX = U_edge;
}

void INTERFACE_2D::ComputeAverageU(int u_flag_target, int u_flag_store) {
	for (int i = 0; i < this->number_gp; i++) {
		this->U_edge_IN[u_flag_store][i] = 0.5 * (this->U_edge_IN[u_flag_target][i] + this->U_edge_EX[u_flag_target][i]);
	}
}