#ifndef CLASS_EDGE_H
#define CLASS_EDGE_H

class EDGE {
protected:
  int ID;

public:
 EDGE(int ID) : ID(ID) {}

  virtual double* get_f_at_gp(double f_bf_coeffs[])=0;
  virtual double* test_against_phi(double f_at_gp[])=0;
};

#endif
