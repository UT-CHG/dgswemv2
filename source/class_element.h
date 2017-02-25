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
	BASIS_GEOM_TRI* basis_geom;

	double*** J_2D_inv_area;
	double* det_J_2D_area;

	double** surface_J_edge;
	double** normal_edge;

	double** area_int_fac_phi;
	double** area_int_fac_dphidx;
	double** area_int_fac_dphidy;
	double*** edge_int_fac_x;
	double*** edge_int_fac_y;

public:
	ELEMENT_TRI(int, int, int, double[], double[], BASIS_TRI*, BASIS_GEOM_TRI*);
	~ELEMENT_TRI();

private:
	void compute_geometry();
	void compute_integration_factors();
};