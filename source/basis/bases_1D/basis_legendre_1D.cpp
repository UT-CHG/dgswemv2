#include "../bases_1D.hpp"

namespace Basis {
Array2D<double> Legendre_1D::GetPhi(const uint p, const std::vector<Point<1>>& points) {
    uint n_pts = points.size();

    Array2D<double> phi(p + 1);

    std::vector<double> l1(n_pts);

    for (uint pt = 0; pt < n_pts; pt++) {
        l1[pt] = points[pt][LocalCoordLin::l1];
    }

    for (uint dof = 0; dof < phi.size(); dof++) {
        phi[dof] = jacobi_polynomial(dof, 0, 0, l1);
    }

    return phi;
}

Array3D<double> Legendre_1D::GetDPhi(const uint p, const std::vector<Point<1>>& points) {
    uint n_pts = points.size();

    Array3D<double> dphi(p + 1);

    std::vector<double> l1(n_pts);

    for (uint pt = 0; pt < n_pts; pt++) {
        l1[pt] = points[pt][LocalCoordLin::l1];
    }

    for (uint dof = 0; dof < dphi.size(); dof++) {
        dphi[dof].push_back(jacobi_polynomial_derivative(dof, 0, 0, l1));
    }

    return dphi;
}

std::pair<bool, Array2D<double>> Legendre_1D::GetMinv(const uint p) {
    std::pair<bool, Array2D<double>> m_inv(true, Array2D<double>(1));  // diagonal

    m_inv.second[0].resize(p + 1);

    for (uint dof = 0; dof < m_inv.second[0].size(); dof++) {
        m_inv.second[0][dof] = (2 * dof + 1) / 2.0;
    }

    return m_inv;
}

void Legendre_1D::ProjectBasisToLinear(const std::vector<double>& u, std::vector<double>& u_lin) {
    u_lin[0] = 0.5 * u[0] - 0.5 * u[1];
    u_lin[1] = 0.5 * u[0] + 0.5 * u[1];
}

void Legendre_1D::ProjectLinearToBasis(const std::vector<double>& u_lin, std::vector<double>& u) {
    u[0] = u_lin[0] + u_lin[1];
    u[1] = -u_lin[0] + u_lin[1];
}
}