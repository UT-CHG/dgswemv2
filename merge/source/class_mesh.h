#ifndef CLASS_MESH_H
#define CLASS_MESH_H

#include "general_definitions.h"

#include "class_basis_geometry.h"
#include "class_element.h"

class MESH {
    friend class PROBLEM;

protected:
    int p;
    int p_geom;

	MasterElement triangle = MasterElement(1);

    std::map<unsigned int, ELEMENT*> elements;
	std::map<unsigned char, std::vector<INTERFACE*>> interfaces;
public:
    MESH(int, int);
    ~MESH();

	void RectangularDomainTest(double, double, int, int, int);

private:
    void InitializeInterfaces();
    void InitializeVTK();
};

#endif