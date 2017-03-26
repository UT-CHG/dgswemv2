#ifndef CLASS_PROBLEM_H
#define CLASS_PROBLEM_H

#include "../class_mesh.h"

class PROBLEM {
private:
	MESH* mesh;
	
	double dt = 0.001;

	std::vector <INTERFACE*> internal_interfaces;
	std::vector <INTERFACE*> ocean_interfaces;
	std::vector <INTERFACE*> land_interfaces;
	std::vector <INTERFACE*> flow_interfaces;

public:
	PROBLEM();
	~PROBLEM();

	void EETimeStepper(int);
	void RK2TimeStepper(int);
	void RK3TimeStepper(int);
	void RK4TimeStepper(int);

private:
	void Timestep();

	void InterfaceFlowAverage(INTERFACE*);
	void ComputeUVA(ELEMENT*);
	void ComputeF(ELEMENT*);
	void LLFNumericalFlux(INTERFACE*);

	void OceanInterfaceSetBC(INTERFACE*);	
	void LandInterfaceSetBC(INTERFACE*);
	void FlowInterfaceSetBC(INTERFACE*);
//	void RadiationInterfaceSetBC(INTERFACE*);

//  void NodeInterfaceReceiveData(INTERFACE*);
//  void NodeInterfaceSendData(INTERFACE*);
	
	void ComputeBoundaryInterfaceF(INTERFACE*);
};

#endif