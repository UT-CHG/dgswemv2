#ifndef CLASS_ELEMENT_2D_H
#define CLASS_ELEMENT_2D_H

#include <map>
#include <vector>
#include <fstream>
#include <string>

#include "../interfaces/class_interface.h"
#include "../class_basis.h"
#include "../class_basis_geometry.h"

class ELEMENT{
    friend class PROBLEM;

private:
    unsigned int ID;

	unsigned char dimension;
	unsigned char element_type;
	unsigned char number_boundaries;

   	unsigned char* boundary_type;
	unsigned int* neighbor_ID;

    double** nodal_coordinates;

    INTERFACE** interfaces;
    bool* interface_owner;

    BASIS* basis;
    BASIS_GEOM* basis_geom = nullptr;
    
	int number_bf;
    int number_bf_geom;

    int number_gp_internal;
    int number_gp_boundary;

    double*** J_inv_t_internal;
    double* det_J_internal;
    double** surface_J_boundary;
    double*** normal_boundary;

	double** phi_internal;
	double*** dphi_internal;
	double*** phi_boundary;
	
    double** internal_int_fac_phi;
    double*** internal_int_fac_dphi;
 
    double*** boundary_int_fac_phi;
    double**** boundary_int_fac_n;

    bool orthogonal;
    double** m_inv;

    double** u;
    std::vector<double**> u_substep;

    double** u_internal;
    double*** u_boundary;

    double* RHS;

public:
	ELEMENT(int, unsigned int, unsigned int*, unsigned char*,
		double**, BASIS*, BASIS_GEOM* basis_geom = nullptr);
    ~ELEMENT();

    std::map<unsigned int, INTERFACE*> CreateInterfaces();
    void AppendInterface(unsigned int, INTERFACE*);
    std::vector<std::pair<unsigned char, INTERFACE*>> GetOwnInterfaces();

    void ComputeInternalU(int);
	void ComputeInternalDU(int, int, int);
    void ComputeBoundaryU(int);
   
    double IntegrationInternalPhi(int, int);
    double IntegrationInternalDPhi(int, int, int);

    double IntegrationBoundaryPhi(int, int);
    double IntegrationBoundaryN(int, int, int);

    void SolveLSE(int);

    void InitializeVTK(std::vector<double*>&, std::vector<unsigned int*>&);
    void WriteCellDataVTK(std::vector<double>&, int);
	void WritePointDataVTK(std::vector<double>&, int);

private:
	void allocate_memory();
    
    void ComputeDPhi(); //rename ComputeDifferentiationFactors;
	void ComputeIntegrationFactors();
    
    void Triangle(unsigned int*, unsigned char*, double**);
    void InitializeVTKTriangle(std::vector<double*>&, std::vector<unsigned int*>&);
	void WriteCellDataVTKTriangle(std::vector<double>&, int);
	void WritePointDataVTKTriangle(std::vector<double>&, int);
};

#endif