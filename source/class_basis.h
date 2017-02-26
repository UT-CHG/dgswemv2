#ifndef CLASS_BASIS_H
#define CLASS_BASIS_H

#include "class_integration.h"

class BASIS_TRI {
private: 
    int p;
    bool orthogonal;

    AREA_INTEGRATION* integration_rule_area;
    LINE_INTEGRATION* integration_rule_line;

    int number_bf;

    double** phi_area;
    double** dphi_dz1_area;
    double** dphi_dz2_area;
    double*** phi_edge;
    double** m_inv;

public:
    BASIS_TRI(int, AREA_INTEGRATION*, LINE_INTEGRATION*);
    ~BASIS_TRI();

    int GetPolynomial();
    int GetNumberBasisFunctions();

    AREA_INTEGRATION* GetIntegrationRuleArea();
    LINE_INTEGRATION* GetIntegrationRuleLine();
    
    double** GetPhiArea();
    double** GetDPhiDZ1Area();
    double** GetDPhiDZ2Area();
    double*** GetPhiEdge();

private:
    void Dubiner();
};

#endif