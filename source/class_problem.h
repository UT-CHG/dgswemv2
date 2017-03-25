#ifndef CLASS_PROBLEM_H
#define CLASS_PROBLEM_H

#include "class_mesh.h"

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
        void ComputeUVA(std::unique_ptr<ELEMENT>&);
        void ComputeF(std::unique_ptr<ELEMENT>&);
	void LLFNumericalFlux(INTERFACE*);

	void LandInterfaceSetBC(INTERFACE*);
	void ComputeBoundaryInterfaceF(INTERFACE*);
};

#endif
