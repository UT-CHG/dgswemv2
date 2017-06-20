#ifndef CLASS_MESH_H
#define CLASS_MESH_H

#include "general_definitions.h"

#include "class_element.h"
#include "class_boundary.h"

class MESH {
protected:
    int p;
	
	Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>* triangle; 
	Shape::StraightTriangle* shape;

    std::map<unsigned int, Element<2, Master::Triangle<Basis::Dubiner_2D, Integration::Dunavant_2D>, 
		Shape::StraightTriangle>*> elements;

	std::map<unsigned char, std::vector<Boundary<>*>> boundaries;
	std::vector<Interface<>*> interfaces_;
public:
	MESH(int p) : p(p) {}
	~MESH();

	void RectangularDomainTest(double, double, int, int, int);

private:
	void InitializeBoundariesInterfaces();
    void InitializeVTK();
};

#endif