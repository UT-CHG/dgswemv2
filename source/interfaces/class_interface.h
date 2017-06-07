#ifndef CLASS_INTERFACE_2D_H
#define CLASS_INTERFACE_2D_H

#include "../general_definitions.h"

class INTERFACE{
    friend class PROBLEM;

private:
    unsigned char dimension;
    bool boundary;
    int number_gp;

    double** normal;

    double** u_boundary_in;
    double** u_boundary_ex;

public:
    INTERFACE(unsigned char, int, double**, double**, bool boundary = false);
    ~INTERFACE();

	void SetPointerEX(double** u_boundary) { this->u_boundary_ex = u_boundary; };
};

#endif