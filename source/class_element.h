#ifndef CLASS_ELEMENT_H
#define CLASS_ELEMENT_H

class ELEMENT {
protected:
  int ID;
  int p;

public:
 ELEMENT(int ID, int p) : ID(ID), p(p) {}

  virtual double* get_f_at_gp(double[] f_bf_coeffs);

  virtual void test_against_phi(double[] f_at_gp);
  virtual void test_against_dphidx(double[] f_at_gp);
  virtual void test_against_dphidy(double[] f_at_gp);

  virtual void invert_mass_matrix(double[] f_bf);
}

#endif
