#ifndef CLASS_MESH_H
#define CLASS_MESH_H

#include <map>
#include <vector>

#include "general_definitions.h"

#include "class_integration.h"
#include "class_basis.h"
#include "class_basis_geometry.h"

#include "class_element.h"

class MESH {
    friend class PROBLEM;

protected:
    int p;
    int p_geom;

    std::map<unsigned int, ELEMENT*> elements;
    std::map<unsigned char, std::vector<INTERFACE*>> interfaces;

    std::map<unsigned char, INTEGRATION*> boundary_rules;
    std::map<unsigned char, INTEGRATION*> internal_rules;
    
    std::map<unsigned char, BASIS*> bases;
    std::map<unsigned char, BASIS_GEOM*> geometric_bases;

public:
    MESH(int, int);
    ~MESH();

	void RectangularDomainTest(double, double, int, int, int);

private:
    void InitializeInterfaces();
    void InitializeVTK();
};

#endif