#ifndef CLASS_INTEGRATION_H
#define CLASS_INTEGRATION_H

class INTEGRATION_1D {

private:
    int p;
    int number_gp;
    double* w;
    double* z;

public:
    INTEGRATION_1D(int);
    ~INTEGRATION_1D();

    int GetPolynomial();
    int GetNumberGP();
    double* GetWeight();
    double* GetZ();

private:
    void GaussLegendre();
};

class INTEGRATION_2D {
private: 
    int p;
    int number_gp;
    double* w;
    double* z1;
    double* z2;

public:
    INTEGRATION_2D(int);
    ~INTEGRATION_2D();

    int GetPolynomial();
    int GetNumberGP();
    double* GetWeight();
    double* GetZ1();
    double* GetZ2();

private:	
    void Dunavant();
};

#endif