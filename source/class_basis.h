#ifndef CLASS_BASIS_H
#define CLASS_BASIS_H

#include "general_definitions.h"
#include "basis_functions/basis_functions.h"

#include "class_integration.h"

class BASIS_2D {
private: 
    int p;
    bool orthogonal;

    INTEGRATION_1D* integration_rule_line;
    INTEGRATION_2D* integration_rule_area;

    int number_bf;

    double** phi_area;
    double** dphi_dz1_area;
    double** dphi_dz2_area;
    double*** phi_edge;
    double** m_inv;

public:
    BASIS_2D(int, int, INTEGRATION_1D*, INTEGRATION_2D*);
    ~BASIS_2D();

    int GetPolynomial();
    int GetNumberBasisFunctions();

    INTEGRATION_1D* GetIntegrationRuleLine();
    INTEGRATION_2D* GetIntegrationRuleArea();
    
    double** GetPhiArea();
    double** GetDPhiDZ1Area();
    double** GetDPhiDZ2Area();
    double*** GetPhiEdge();
    bool GetOrthogonal();
    double** GetMInv();

private:
    void Dubiner();
};

#endif