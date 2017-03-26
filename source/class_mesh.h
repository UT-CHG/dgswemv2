#ifndef CLASS_MESH_H
#define CLASS_MESH_H

#include <map>
#include <vector>

#include "general_definitions.h"

#include "class_integration.h"
#include "class_basis.h"
#include "class_basis_geometry.h"

#include "class_element.h"
#include "elements\class_element_2D.h"

class MESH {
    friend class PROBLEM;

protected:
    int p;
    int p_geom;

    std::map<unsigned int, ELEMENT*> elements;
    std::map<unsigned char, std::vector<INTERFACE*>> interfaces;
    
    std::map<unsigned char, INTEGRATION_1D*> line_rules;
    std::map<unsigned char, INTEGRATION_2D*> area_rules;
    
    std::map<unsigned char, BASIS_2D*> bases_2D;
    std::map<unsigned char, BASIS_GEOM_2D*> geometric_bases_2D;

public:
    MESH(int, int);
    ~MESH();

    //void solve();

protected:
    void InitializeElements();
    void InitializeInterfaces();
};

#endif