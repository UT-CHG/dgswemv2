#include "../bases_1D.hpp"

namespace Basis {
DynMatrix<double> Legendre_1D::GetPhi(const uint p, const std::vector<Point<1>>& points) {
    uint ndof = p + 1;
    uint npt  = points.size();

    DynMatrix<double> phi(ndof, npt);

    std::vector<double> l1(npt);

    for (uint pt = 0; pt < npt; ++pt) {
        l1[pt] = points[pt][LocalCoordLin::l1];
    }

    for (uint dof = 0; dof < ndof; ++dof) {
        row(phi, dof) = transpose(jacobi_polynomial(dof, 0, 0, l1));
    }

    return phi;
}

std::array<DynMatrix<double>, 1> Legendre_1D::GetDPhi(const uint p, const std::vector<Point<1>>& points) {
    uint ndof = p + 1;
    uint npt  = points.size();

    std::array<DynMatrix<double>, 1> dphi;

    DynMatrix<double> dphi_dl1(ndof, npt);

    std::vector<double> l1(npt);

    for (uint pt = 0; pt < npt; ++pt) {
        l1[pt] = points[pt][LocalCoordLin::l1];
    }

    for (uint dof = 0; dof < ndof; ++dof) {
        row(dphi_dl1, dof) = transpose(jacobi_polynomial_derivative(dof, 0, 0, l1));
    }

    dphi[LocalCoordLin::l1] = dphi_dl1;

    return dphi;
}

DynMatrix<double> Legendre_1D::GetMinv(const uint p) {
    uint ndof = p + 1;

    DynMatrix<double> m_inv(ndof, ndof);

    set_constant(m_inv, 0.0);

    for (uint dof = 0; dof < ndof; ++dof) {
        m_inv(dof, dof) = (2 * dof + 1) / 2.0;
    }

    return m_inv;
}

DynMatrix<double> Legendre_1D::ProjectBasisToLinear(const DynMatrix<double>& u) {
    std::cout << "Legendre_1D::ProjectBasisToLinear not implemented!" << std::endl;
    abort();
    /*u_lin[0] = 0.5 * u[0] - 0.5 * u[1];
    u_lin[1] = 0.5 * u[0] + 0.5 * u[1];*/
}

DynMatrix<double> Legendre_1D::ProjectLinearToBasis(const uint ndof, const DynMatrix<double>& u_lin) {
    std::cout << "Legendre_1D::ProjectLinearToBasis not implemented!" << std::endl;
    abort();
    /*u[0] = u_lin[0] + u_lin[1];
    u[1] = -u_lin[0] + u_lin[1];*/
}
}