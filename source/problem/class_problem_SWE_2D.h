#ifndef CLASS_PROBLEM_H
#define CLASS_PROBLEM_H

#include "../class_mesh.h"

class PROBLEM {
private:
    MESH* mesh;
    
	double t = 0;
    double dt = 0.0001;

public:
    PROBLEM();
    ~PROBLEM();

    void EETimeStepper(int);
    void RK2TimeStepper(int);
    void RK3TimeStepper(int);
    void RK4TimeStepper(int);

    void WriteDataVTK();

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