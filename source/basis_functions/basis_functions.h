void jacobi_polynomial(int n, int a, int b, int number_gp, const double* x, double* P);
void jacobi_polynomial_derivative(int n, int a, int b, int number_gp, const double* x, double* dP);

void dubiner_2d_phi(int p, int q, int m, int number_gp, const double* n1, const double* n2, double** phi_pq);
void dubiner_2d_dphi(int p, int q, int m, int number_gp, const double* n1, const double* n2, double*** dphi_dz);
void dubiner_2d_test(int p, int number_gp, const double* const* phi_area, const double* w);