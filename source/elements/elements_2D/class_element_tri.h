#ifndef CLASS_ELEMENT_TRI_H
#define CLASS_ELEMENT_TRI_H

#include "../../general_definitions.h"
#include "../class_element_2D.h"
#include "../../class_interface.h"

class ELEMENT_TRI : public ELEMENT_2D {
public:
	ELEMENT_TRI(unsigned int, unsigned int[], unsigned char[],
		double[], double[], BASIS_2D*, BASIS_GEOM_2D* basis_geom = nullptr);
	~ELEMENT_TRI();

	std::map<unsigned int, INTERFACE*> CreateInterfaces();
	void AppendInterface(unsigned int, INTERFACE*);

	std::vector<std::pair<unsigned char, INTERFACE*>> GetOwnInterfaces();

protected:
	void ComputeGeometry();
	void ComputeIntegrationFactors();
};

#endif