#ifndef CLASS_PROBLEM_H
#define CLASS_PROBLEM_H

#include "class_mesh_v2.h"

class PROBLEM {
private:
	MESH* mesh;

	std::vector <INTERFACE*> internal_interfaces;
	std::vector <INTERFACE*> land_interfaces;

public:
	PROBLEM();
	~PROBLEM();

	void Timestep();

private:
	void InterfaceFlowAverage(INTERFACE*);
	void ComputeUVA(ELEMENT*);
	void ComputeF(ELEMENT*);
	void LLFNumericalFlux(INTERFACE*);

	void LandInterfaceSetBC(INTERFACE*);
	void ComputeBoundaryInterfaceF(INTERFACE*);
};

#endif