#ifndef CLASS_BASIS_H
#define CLASS_BASIS_H

#include "general_definitions.h"
#include "basis_functions/basis_functions.h"

#include "class_integration.h"

class BASIS {
private: 
	bool orthogonal;
	int dimension;
	int number_boundaries;
	int p;

    INTEGRATION* integration_rule_boundary;
    INTEGRATION* integration_rule_internal;

    int number_bf;

    double** phi_internal;
    double*** dphi_dz_internal;
	double*** phi_boundary;
    double** m_inv;

    double** phi_postprocessor_cell;
	double** phi_postprocessor_point;

public:
    BASIS(int, int, INTEGRATION*, INTEGRATION*);
    ~BASIS();

	int GetPolynomial() { return this->p; }
	int GetNumberBasisFunctions() { return this->number_bf; }

	INTEGRATION* GetIntegrationRuleBoundary() { return this->integration_rule_boundary; }
	INTEGRATION* GetIntegrationRuleInternal() { return this->integration_rule_internal; }

	double** GetPhiInternal() { return this->phi_internal; };
	double*** GetDPhiDZInternal() { return this->dphi_dz_internal; };
	double*** GetPhiBoundary() { return this->phi_boundary; };

	bool GetOrthogonal() { return this->orthogonal; };
	double** GetMInv() { return this->m_inv; };

	double** GetPhiPostProcessorCell() { return this->phi_postprocessor_cell; };
	double** GetPhiPostProcessorPoint() { return this->phi_postprocessor_point; };

private:
	void allocate_memory(int, int);

    void Dubiner2D();
};

#endif