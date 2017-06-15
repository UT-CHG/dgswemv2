#ifndef CLASS_ELEMENT_H
#define CLASS_ELEMENT_H

#include "general_definitions.h"

#include "class_interface.h"
#include "class_basis_geometry.h"
#include "class_master_element.h"

template<int dimension = 2, int element_type = TRIANGLE, 
class basis_type = Dubiner_2D, 
class integration_int_type = Dunavant_2D, 
class integration_bound_type = GaussLegendre_1D>
class ELEMENT{
    friend class PROBLEM;

private:
	MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>& master;

    unsigned int ID;

	unsigned char number_boundaries;

   	std::vector<unsigned char> boundary_type;
	std::vector<unsigned int> neighbor_ID;
	std::pair<std::vector<bool>, std::vector<INTERFACE*>> interfaces;

	Array2D<double> nodal_coordinates;
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
 
	Array2D<double> u;
    Array3D<double> u_substep;

    Array2D<double> u_internal;
    double*** u_boundary;

    std::vector<double> RHS;

public:
	ELEMENT(MasterElement<dimension, element_type, basis_type, integration_int_type, integration_bound_type>&,
		unsigned int, std::vector<unsigned int>&, std::vector<unsigned char>&, Array2D<double>&, BASIS_GEOM* basis_geom = nullptr);

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

    void InitializeVTK(std::vector<Point<3>>&, Array2D<unsigned int>&);
    void WriteCellDataVTK(std::vector<double>&, int);
	void WritePointDataVTK(std::vector<double>&, int);

private:
	void allocate_memory();
    
    void ComputeDifferentiationFactors();
	void ComputeIntegrationFactors();
    
    void Triangle();
    void InitializeVTKTriangle(std::vector<Point<3>>&, Array2D<unsigned int>&);
	void WriteCellDataVTKTriangle(std::vector<double>&, int);
	void WritePointDataVTKTriangle(std::vector<double>&, int);
};

#include "class_element.tpp"
#include "elements/elements_2D/element_tri.tpp"

#endif