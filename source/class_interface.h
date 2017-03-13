#ifndef CLASS_INTERFACE_2D_H
#define CLASS_INTERFACE_2D_H

class INTERFACE_2D {	
private:
    double** U_edge_IN;
    double** U_edge_EX;

public:
    INTERFACE_2D(double**);
    ~INTERFACE_2D();
};

#endif