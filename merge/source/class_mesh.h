#ifndef CLASS_MESH_H
#define CLASS_MESH_H

#include "problem/SWE/swe_data.hpp"
#include "problem/SWE/swe_kernels.hpp"
#include "problem/SWE/swe_boundary_conditions.hpp"
#include "problem/SWE/swe_definitions.hpp"

#include "general_definitions.h"
#include "stepper.hpp"

#include "class_element.h"
#include "class_boundary.h"

class MESH {
private:
    int p;
	
	Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>* triangle; 
	Shape::StraightTriangle* shape;

    std::map<unsigned int, Element<2, Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>, 
		Shape::StraightTriangle>*> elements;

	std::map<unsigned char, std::vector<Boundary<>*>> boundaries;
	std::vector<Interface<>*> interfaces;

public:
	MESH(int p) : p(p), snapshot(0) {}
	~MESH();

	void RectangularDomainTest(double, double, int, int, int);
	void Solve();
	void WriteDataVTK();
	int snapshot;

private:
	void InitializeBoundariesInterfaces();
    void InitializeVTK();
};

#endif