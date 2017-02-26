#ifndef CLASS_INTEGRATION_H
#define CLASS_INTEGRATION_H

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
    double* GetWeight();
    double* GetZ1();
    double* GetZ2();

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
    double* GetWeight();
    double* GetZ();

private:
    void GaussLegendre();
};

#endif