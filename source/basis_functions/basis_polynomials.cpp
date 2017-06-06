#include "basis_functions.h"

void jacobi_polynomial(int n, int a, int b, int number_gp, const double* x, double* P) {
    if (n == 0) {
        for (int i = 0; i < number_gp; i++) { P[i] = 1; }
    }
    else if (n == 1) {
        for (int i = 0; i < number_gp; i++) { P[i] = (a - b + (a + b + 2)*x[i]) / 2; }
    }
    else {
        double** J = new double*[3];

        J[0] = new double[number_gp];
        J[1] = new double[number_gp];
        J[2] = new double[number_gp];

        for (int i = 0; i < number_gp; i++) {
            J[0][i] = 1;
            J[1][i] = (a - b + (a + b + 2)*x[i]) / 2;
        }

        double a1;
        double a2;
        double a3;
        double a4;

        for (int i = 1; i < n; i++) {
            a1 = 2 * (i + 1)*(i + a + b + 1)*(2 * i + a + b);
            a2 = (2 * i + a + b + 1)*(a*a - b*b);
            a3 = (2 * i + a + b)*(2 * i + a + b + 1)*(2 * i + a + b + 2);
            a4 = 2 * (i + a)*(i + b)*(2 * i + a + b + 2);

            for (int j = 0; j < number_gp; j++) {
                J[(i + 1) % 3][j] = ((a2 + a3*x[j])*J[i % 3][j] - a4*J[(i - 1) % 3][j]) / a1;
            }
        }

        for (int i = 0; i < number_gp; i++) {
            P[i] = J[n % 3][i];
        }

        delete[] J[0];
        delete[] J[1];
        delete[] J[2];

        delete[] J;
    }
}

void jacobi_polynomial_derivative(int n, int a, int b, int number_gp, const double* x, double* dP) {
    if (n == 0) {
        for (int i = 0; i < number_gp; i++) { dP[i] = 0; }
    }
    else if (n == 1) {
        for (int i = 0; i < number_gp; i++) { dP[i] = (a + b + n + 1) / 2.0; }
    }
    else {
        jacobi_polynomial(n - 1, a + 1, b + 1, number_gp, x, dP);

        for (int i = 0; i < number_gp; i++) { dP[i] = dP[i] * (a + b + n + 1) / 2.0; }
    }
}