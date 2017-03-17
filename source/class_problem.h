#ifndef CLASS_PROBLEM_H
#define CLASS_PROBLEM_H

#include "class_mesh_v2.h"

class PROBLEM {
private:
	MESH* mesh;

public:
	PROBLEM();
	~PROBLEM();

	void Timestep();

	void ComputeUVA(ELEMENT*);

	void ComputeF(ELEMENT*);

	void InternalInterfaceDiffusion(INTERFACE* interface);

	void LLFNumericalFlux(INTERFACE*);
};

#endif