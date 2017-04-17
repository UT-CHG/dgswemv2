#ifndef CLASS_ELEMENT_H
#define CLASS_ELEMENT_H

#include <map>
#include <vector>
#include <fstream>
#include <string>

#include "class_interface.h"

class ELEMENT {
    friend class PROBLEM;

protected:
    unsigned int ID;

	unsigned char element_type;
    unsigned char number_interfaces;

	unsigned int* neighbor_ID;
	unsigned char* boundary_type;

	double* nodal_coordinates_x;
	double* nodal_coordinates_y;

	int number_bf;
    int number_bf_geom;

    int number_gp_internal;
    int number_gp_boundary;

    double** u;

    std::vector<double**> u_substep;

    double** u_internal;
    double*** u_boundary;

    double* RHS;
    
public:
    ELEMENT(int ID) : ID(ID) {}
    virtual ~ELEMENT() = default;

    virtual std::map<unsigned int, INTERFACE*> CreateInterfaces() = 0;
    virtual void AppendInterface(unsigned int, INTERFACE*) = 0;

    virtual std::vector<std::pair<unsigned char, INTERFACE*>> GetOwnInterfaces() = 0;

    virtual void ComputeInternalU(int) = 0;
    virtual void ComputeBoundaryU(int) = 0;

	virtual void ComputeInternalDUDX(int, int) = 0;
	virtual void ComputeInternalDUDY(int, int) = 0;

    virtual double IntegrationInternalPhi(int, int) = 0;
    virtual double IntegrationInternalDPhiDX(int, int) = 0;
    virtual double IntegrationInternalDPhiDY(int, int) = 0;

    virtual double IntegrationBoundaryPhi(int, int) = 0;
    virtual double IntegrationBoundaryNX(int, int) = 0;
    virtual double IntegrationBoundaryNY(int, int) = 0;

    virtual void SolveLSE(int) = 0;

    virtual void InitializeVTK(std::vector<double*>&, std::vector<unsigned int*>&) = 0;
    virtual void WriteCellDataVTK(std::vector<double>&, int) = 0;
	virtual void WritePointDataVTK(std::vector<double>&, int) = 0;

    //virtual double test_against_phi(double f_at_gp[]) = 0;
    //virtual double test_against_dphidx(double f_at_gp[]) = 0;
    //virtual double test_against_dphidy(double f_at_gp[]) = 0;

    //virtual double* invert_mass_matrix(double f_bf[]) = 0;
};

#endif