#ifndef CLASS_ELEMENT_H
#define CLASS_ELEMENT_H

#include "class_basis.h"
#include "class_basis_geometry.h"

class ELEMENT_TRI {
private:
    int ID;

    int p;
    int p_geom;

    double* nodal_coordinates_x;
    double* nodal_coordinates_y;

    BASIS_TRI* basis;
    
    BASIS_GEOM_TRI* basis_geom = nullptr;
    double** M_inv = nullptr;

    double*** J_inv_t_area;
    double* det_J_area;

    double** surface_J_edge;
    double** normal_edge_x;
    double** normal_edge_y;

    double** area_int_fac_phi;
    double** area_int_fac_dphidx;
    double** area_int_fac_dphidy;
    double*** edge_int_fac_x;
    double*** edge_int_fac_y;

public:
    ELEMENT_TRI(int, int, int, double[], double[], BASIS_TRI*, BASIS_GEOM_TRI* basis_geom = nullptr);
    ~ELEMENT_TRI();

private:
    void ComputeGeometry();
    void ComputeIntegrationFactors();
};

#endif
