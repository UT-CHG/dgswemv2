#include <iostream>

#include "class_interface.h"

INTERFACE_2D::INTERFACE_2D(double** U_edge) {
    this->U_edge_IN = U_edge;
}

INTERFACE_2D::~INTERFACE_2D() {

}