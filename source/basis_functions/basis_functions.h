void jacobi_polynomial(int n, int a, int b, int number_gp, double x[], double P[]);
void jacobi_polynomial_derivative(int n, int a, int b, int number_gp, double x[], double dP[]);

void dubiner_phi(int p, int q, int number_gp, double n1[], double n2[], double phi_pq[]);
void dubiner_dphi(int p, int q, int number_gp, double n1[], double n2[], double dphi_dz1[], double dphi_dz2[]);
void dubiner_test(int p, int number_gp, double** phi_area, double* w);

TEST
