#include "../basis_polynomials.hpp"

namespace Basis {
DynVector<double> jacobi_polynomial(const uint n, const uint a, const uint b, const std::vector<double>& x) {
    uint npt = x.size();

    DynVector<double> P(npt);

    if (n == 0) {
        for (uint pt = 0; pt < npt; ++pt) {
            P[pt] = 1;
        }
    } else if (n == 1) {
        for (uint pt = 0; pt < npt; ++pt) {
            P[pt] = (a - b + (a + b + 2) * x[pt]) / 2;
        }
    } else {
        DynMatrix<double> J(3, npt);

        for (uint pt = 0; pt < npt; ++pt) {
            J(0, pt) = 1.0;
            J(1, pt) = (a - b + (a + b + 2) * x[pt]) / 2;
            J(2, pt) = 0.0;
        }

        for (uint i = 1; i < n; ++i) {
            double a1 = 2 * (i + 1) * (i + a + b + 1) * (2 * i + a + b);
            double a2 = (2 * i + a + b + 1) * (a * a - b * b);
            double a3 = (2 * i + a + b) * (2 * i + a + b + 1) * (2 * i + a + b + 2);
            double a4 = 2 * (i + a) * (i + b) * (2 * i + a + b + 2);

            for (uint pt = 0; pt < npt; ++pt) {
                J((i + 1) % 3, pt) = ((a2 + a3 * x[pt]) * J(i % 3, pt) - a4 * J((i - 1) % 3, pt)) / a1;
            }
        }

        for (uint pt = 0; pt < npt; ++pt) {
            P[pt] = J(n % 3, pt);
        }
    }

    return P;
}

DynVector<double> jacobi_polynomial_derivative(const uint n, const uint a, const uint b, const std::vector<double>& x) {
    uint npt = x.size();

    DynVector<double> dP(npt);

    if (n == 0) {
        for (uint pt = 0; pt < npt; ++pt) {
            dP[pt] = 0.0;
        }
    } else if (n == 1) {
        for (uint pt = 0; pt < npt; ++pt) {
            dP[pt] = (a + b + n + 1) / 2.0;
        }
    } else {
        dP = jacobi_polynomial(n - 1, a + 1, b + 1, x);

        for (uint pt = 0; pt < npt; ++pt) {
            dP[pt] *= (a + b + n + 1) / 2.0;
        }
    }

    return dP;
}
}