#ifndef CLASS_MESH_H
#define CLASS_MESH_H

#include "general_definitions.h"

#include "class_element.h"
#include "boundary/class_boundary.h"

class MESH {
    friend class PROBLEM;

protected:
    int p;
    int p_geom;
	
	MasterElement<2, TRIANGLE, Basis::Dubiner_2D, Integration::Dunavant_2D, Integration::GaussLegendre_1D>* triangle; 
	Shape::StraightTriangle* shape;

    std::map<unsigned int, Element<2, TRIANGLE, Basis::Dubiner_2D, 
		Integration::Dunavant_2D, Integration::GaussLegendre_1D, Shape::StraightTriangle>*> elements;
	std::map<unsigned char, std::vector<INTERFACE*>> interfaces;

	std::map<unsigned char, std::vector<Boundary<>*>> boundaries;

public:
    MESH(int, int);
	~MESH();

	void RectangularDomainTest(double, double, int, int, int);

private:
	void InitializeBoundaries();

    void InitializeInterfaces();
    void InitializeVTK();
};

#endif