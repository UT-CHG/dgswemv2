class AREA_INTEGRATION {
private: 
	int p;
	int number_gp;
	double* w;
	double* z1;
	double* z2;

public:
	AREA_INTEGRATION(int);
	~AREA_INTEGRATION();

	int GetPolynomial();
	int GetNumberGP();
	double GetWeight(int);
	double GetZ1(int);
	double GetZ2(int);

private:	
	void Dunavant();
};

class LINE_INTEGRATION {

private: 
	int p;
	int number_gp;
	double* w;
	double* z;

public:
	LINE_INTEGRATION(int);
	~LINE_INTEGRATION();

	int GetPolynomial();
	int GetNumberGP();
	double GetWeight(int);
	double GetZ(int);

private:
	void GaussLegendre();
};