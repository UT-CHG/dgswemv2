#ifndef CLASS_INTERFACE_2D_H
#define CLASS_INTERFACE_2D_H

#include "../class_interface.h"

class INTERFACE_2D : public INTERFACE {
private:
	int number_gp;

	double** U_edge_IN;
	double** U_edge_EX;

public:
	INTERFACE_2D(int, double**);
	~INTERFACE_2D() = default;

	void SetPointerEX(double**);

	void ComputeAverageU(int, int);
};

#endif