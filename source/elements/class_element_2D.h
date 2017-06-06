#ifndef CLASS_ELEMENT_2D_H
#define CLASS_ELEMENT_2D_H

#include "../class_element.h"
#include "../interfaces/class_interface_2D.h"
#include "../class_basis.h"
#include "../class_basis_geometry.h"

class ELEMENT_2D : public ELEMENT {
protected:
    INTERFACE_2D** interfaces;
    bool* interface_owner;

    BASIS* basis;
    BASIS_GEOM* basis_geom = nullptr;
    
    bool orthogonal;
    double** m_inv;

    double*** J_inv_t_internal;
    double* det_J_internal;

    double** surface_J_boundary;
    double*** normal_boundary;

	double** phi_internal;
	double*** phi_boundary;

	double*** dphi_internal;
	
    double** internal_int_fac_phi;
    double*** internal_int_fac_dphi;
    
    double*** boundary_int_fac_phi;
    double**** boundary_int_fac_n;
    

public:
	ELEMENT_2D(int, unsigned int, unsigned int*, unsigned char*,
		double*, double*, BASIS*, BASIS_GEOM* basis_geom = nullptr);
    ~ELEMENT_2D();

	void allocate_memory();

    void Triangle(unsigned int*, unsigned char*, double*, double*);

	void ComputeDPhi();
	void ComputeIntegrationFactors();

    std::map<unsigned int, INTERFACE*> CreateInterfaces();
    void AppendInterface(unsigned int, INTERFACE*);

    std::vector<std::pair<unsigned char, INTERFACE*>> GetOwnInterfaces();

    void ComputeInternalU(int);
    void ComputeBoundaryU(int);

	void ComputeInternalDUDX(int,int);
	void ComputeInternalDUDY(int,int);

    
    double IntegrationInternalPhi(int, int);
    double IntegrationInternalDPhiDX(int, int);
    double IntegrationInternalDPhiDY(int, int);

    double IntegrationBoundaryPhi(int, int);
    double IntegrationBoundaryNX(int, int);
    double IntegrationBoundaryNY(int, int);

    void SolveLSE(int);

    void InitializeVTK(std::vector<double*>&, std::vector<unsigned int*>&);
    void WriteCellDataVTK(std::vector<double>&, int);
	void WritePointDataVTK(std::vector<double>&, int);

    void InitializeVTKTriangle(std::vector<double*>&, std::vector<unsigned int*>&);
	void WriteCellDataVTKTriangle(std::vector<double>&, int);
	void WritePointDataVTKTriangle(std::vector<double>&, int);

};

#endif