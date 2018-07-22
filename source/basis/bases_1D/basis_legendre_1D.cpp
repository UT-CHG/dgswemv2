#include "../bases_1D.hpp"

namespace Basis {
DynMatrix<double> Legendre_1D::GetPhi(const uint p, const DynVector<Point<1>>& points) {
    uint ndof = p + 1;
    uint npt  = points.size();

    DynMatrix<double> phi(ndof, npt);

    DynVector<double> l1(npt);

    for (uint pt = 0; pt < npt; pt++) {
        l1[pt] = points[pt][LocalCoordLin::l1];
    }

    for (uint dof = 0; dof < ndof; dof++) {
        row(phi, dof) = transpose(jacobi_polynomial(dof, 0, 0, l1));
    }

    return phi;
}

StatVector<DynMatrix<double>, 1> Legendre_1D::GetDPhi(const uint p, const DynVector<Point<1>>& points) {
    uint ndof = p + 1;
    uint npt  = points.size();

    StatVector<DynMatrix<double>, 1> dphi;

    DynMatrix<double> dphi_dl1(ndof, npt);

    DynVector<double> l1(npt);

    for (uint pt = 0; pt < npt; pt++) {
        l1[pt] = points[pt][LocalCoordLin::l1];
    }

    for (uint dof = 0; dof < ndof; dof++) {
        row(dphi_dl1, dof) = transpose(jacobi_polynomial_derivative(dof, 0, 0, l1));
    }

    dphi[LocalCoordLin::l1] = dphi_dl1;

    return dphi;
}

DynMatrix<double> Legendre_1D::GetMinv(const uint p) {
    uint ndof = p + 1;

    DynMatrix<double> m_inv(ndof, ndof, 0.0);

    for (uint dof = 0; dof < ndof; dof++) {
        m_inv(dof, dof) = (2 * dof + 1) / 2.0;
    }

    return m_inv;
}

template <typename InputArrayType>
inline decltype(auto) Legendre_1D::ProjectBasisToLinear(const InputArrayType& u) {
    /*u_lin[0] = 0.5 * u[0] - 0.5 * u[1];
    u_lin[1] = 0.5 * u[0] + 0.5 * u[1];*/
}

template <typename InputArrayType>
inline decltype(auto) Legendre_1D::ProjectLinearToBasis(const uint ndof, const InputArrayType& u_lin) {
    /*u[0] = u_lin[0] + u_lin[1];
    u[1] = -u_lin[0] + u_lin[1];*/
}
}