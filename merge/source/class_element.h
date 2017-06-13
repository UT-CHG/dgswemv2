#ifndef CLASS_ELEMENT_2D_H
#define CLASS_ELEMENT_2D_H

#include "general_definitions.h"

#include "class_interface.h"
#include "class_basis_geometry.h"
#include "class_master_element.h"

class ELEMENT{
    friend class PROBLEM;

private:
	int element_type;
    unsigned int ID;

	unsigned char dimension;
	unsigned char number_boundaries;

   	unsigned char* boundary_type;
	unsigned int* neighbor_ID;

    double** nodal_coordinates;

	INTERFACE** interfaces;
    bool* interface_owner;

	MasterElement& master;

	std::pair<bool, Array2D<double>>& m_inv;

	Array2D<double>& phi_internal;
	Array3D<double>& phi_boundary;
	Array2D<double>& phi_postprocessor_cell;
	Array2D<double>& phi_postprocessor_point;

    BASIS_GEOM* basis_geom = nullptr;
    
	int number_bf;
    int number_bf_geom;

    int number_gp_internal;
    int number_gp_boundary;

	std::vector<double> det_J_internal;
	Array3D<double> J_inv_internal;
	Array2D<double> surface_J_boundary;
	Array3D<double> normal_boundary;

	Array3D<double> internal_dphi_fac;	
    Array2D<double> internal_int_fac_phi;
	Array3D<double> internal_int_fac_dphi;
	Array3D<double> boundary_int_fac_phi;
 
	double** u;
    std::vector<double**> u_substep;

    double** u_internal;
    double*** u_boundary;

    double* RHS;

public:
	ELEMENT(MasterElement&,
		unsigned int, unsigned int*, unsigned char*, double**, BASIS_GEOM* basis_geom = nullptr);

    std::map<unsigned int, INTERFACE*> CreateInterfaces();
    void AppendInterface(unsigned int, INTERFACE*);
    std::vector<std::pair<unsigned char, INTERFACE*>> GetOwnInterfaces();

    void ComputeInternalU(int);
	void ComputeInternalDU(int, int, int);
    void ComputeBoundaryU(int);
   
    double IntegrationInternalPhi(int, int);
    double IntegrationInternalDPhi(int, int, int);

    double IntegrationBoundaryPhi(int, int);

    void SolveLSE(int);

    void InitializeVTK(std::vector<double*>&, std::vector<unsigned int*>&);
    void WriteCellDataVTK(std::vector<double>&, int);
	void WritePointDataVTK(std::vector<double>&, int);

private:
	void allocate_memory();
    
    void ComputeDifferentiationFactors();
	void ComputeIntegrationFactors();
    
    void Triangle();
    void InitializeVTKTriangle(std::vector<double*>&, std::vector<unsigned int*>&);
	void WriteCellDataVTKTriangle(std::vector<double>&, int);
	void WritePointDataVTKTriangle(std::vector<double>&, int);
};

#endif