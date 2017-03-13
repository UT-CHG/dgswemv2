#ifndef CLASS_ELEMENT_TRI_H
#define CLASS_ELEMENT_TRI_H

#include "../../general_definitions.h"
#include "../class_element_2D.h"

class ELEMENT_TRI : public ELEMENT_2D {
public:
	ELEMENT_TRI(unsigned int, unsigned int[], unsigned char[],
		double[], double[], BASIS_2D*, BASIS_GEOM_2D* basis_geom = nullptr);
	~ELEMENT_TRI();

	void CreateInterfaces();

protected:
	void ComputeGeometry();
	void ComputeIntegrationFactors();
};

#endif