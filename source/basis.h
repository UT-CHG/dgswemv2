#include "integration.h"

class BASIS_TRI {
private: 
	int p;
	bool orthogonal;

	AREA_INTEGRATION* integration_rule_area;
	LINE_INTEGRATION* integration_rule_line;	

	int number_bf;
	int number_gp;

	double** phi_area;
	double*** phi_edge;
	double*** dphi_area;
	double** m_inv;

public:
	BASIS_TRI(int);
	~BASIS_TRI();

	int GetPolynomial();

private:
	void Dubiner();
};