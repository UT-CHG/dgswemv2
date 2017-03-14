#ifndef CLASS_ELEMENT_2D_H
#define CLASS_ELEMENT_2D_H

#include <vector>

#include "../class_element.h"
#include "../interfaces/class_interface_2D.h"
#include "../class_basis.h"
#include "../class_basis_geometry.h"

class ELEMENT_2D : public ELEMENT {
protected:
	INTERFACE_2D** interfaces;
	bool* interface_owner;

	BASIS_2D* basis;
    BASIS_GEOM_2D* basis_geom = nullptr;
    double** M_inv = nullptr;

    int number_bf;
    int number_bf_geom;

	int number_gp_area;
	int number_gp_edge;

    double*** J_inv_t_area;
    double* det_J_area;

    double** surface_J_edge;
    double** normal_edge_x;
    double** normal_edge_y;

	double** area_int_fac_phi;
    double** area_int_fac_dphidx;
    double** area_int_fac_dphidy;
    
	double*** edge_int_fac_nx;
    double*** edge_int_fac_ny;

    double** U; 
    std::vector<double**> U_substep;

    double** U_area;
    double*** U_edge;

public:
	ELEMENT_2D(unsigned int, BASIS_2D*, BASIS_GEOM_2D* basis_geom = nullptr);
	~ELEMENT_2D() = default;

    void ComputeInternalU(int);
	void ComputeBoundaryU(int, int);

    double IntegrationInternalPhi(int, int);
    double IntegrationInternalDPhiDX(int, int);
    double IntegrationInternalDPhiDY(int, int);

    double IntegrationBoundaryNX(int, int, int);
    double IntegrationBoundaryNY(int, int, int);
};

#endif