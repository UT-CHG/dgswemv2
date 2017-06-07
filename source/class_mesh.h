#ifndef CLASS_MESH_H
#define CLASS_MESH_H

#include <map>
#include <vector>

#include "general_definitions.h"

#include "class_integration.h"
#include "class_basis.h"
#include "class_basis_geometry.h"

#include "elements\class_element.h"

class MESH {
    friend class PROBLEM;

protected:
    int p;
    int p_geom;

    std::map<unsigned int, ELEMENT*> elements;
    std::map<unsigned char, std::vector<INTERFACE*>> interfaces;

    std::map<unsigned char, INTEGRATION*> line_rules;
    std::map<unsigned char, INTEGRATION*> area_rules;
    
    std::map<unsigned char, BASIS*> bases_2D;
    std::map<unsigned char, BASIS_GEOM*> geometric_bases_2D;

public:
    MESH(int, int);
    ~MESH();

	void RectangularDomainTest(double, double, int, int, int);

private:
    void InitializeInterfaces();
    void InitializeVTK();
};

#endif