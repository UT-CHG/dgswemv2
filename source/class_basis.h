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

    double** phi_postprocessor_cell;
	double** phi_postprocessor_point;

public:
    BASIS_2D(int, int, INTEGRATION_1D*, INTEGRATION_2D*);
    ~BASIS_2D();

	int GetPolynomial() { return this->p; }
	int GetNumberBasisFunctions() { return this->number_bf; }

	INTEGRATION_1D* GetIntegrationRuleLine() { return this->integration_rule_line; }
	INTEGRATION_2D* GetIntegrationRuleArea() { return this->integration_rule_area; }

	double** GetPhiArea() { return this->phi_area; };
	double** GetDPhiDZ1Area() { return this->dphi_dz1_area; };
	double** GetDPhiDZ2Area() { return this->dphi_dz2_area; };
	double*** GetPhiEdge() { return this->phi_edge; };

	bool GetOrthogonal() { return this->orthogonal; };
	double** GetMInv() { return this->m_inv; };

	double** GetPhiPostProcessorCell() { return this->phi_postprocessor_cell; };
	double** GetPhiPostProcessorPoint() { return this->phi_postprocessor_point; };

private:
    void Dubiner();
};

#endif