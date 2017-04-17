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

    BASIS_2D* basis;
    BASIS_GEOM_2D* basis_geom = nullptr;
    
    bool orthogonal;
    double** m_inv;

    double*** J_inv_t_area;
    double* det_J_area;

    double** surface_J_edge;
    double** normal_edge_x;
    double** normal_edge_y;

	double** phi_area;
	double*** phi_edge;

	double** dphi_dx_area;
	double** dphi_dy_area;

    double** area_int_fac_phi;
    double** area_int_fac_dphidx;
    double** area_int_fac_dphidy;
    
    double*** edge_int_fac_phi;
    double*** edge_int_fac_nx;
    double*** edge_int_fac_ny;


public:
	ELEMENT_2D(int, unsigned int, unsigned int*, unsigned char*,
		double*, double*, BASIS_2D*, BASIS_GEOM_2D* basis_geom = nullptr);
    ~ELEMENT_2D();

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