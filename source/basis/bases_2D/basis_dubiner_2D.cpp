#include "../bases_2D.hpp"

namespace Basis {
DynMatrix<double> Dubiner_2D::GetPhi(const uint p, const std::vector<Point<2>>& points) {
    uint ndof = (p + 1) * (p + 2) / 2;
    uint npt  = points.size();

    DynMatrix<double> phi(ndof, npt);

    std::vector<double> n1(npt);
    std::vector<double> n2(npt);

    for (uint pt = 0; pt < npt; ++pt) {
        n1[pt] = 2 * (1 + points[pt][LocalCoordTri::z1]) / (1 - points[pt][LocalCoordTri::z2]) - 1;
        n2[pt] = points[pt][LocalCoordTri::z2];
    }

    for (uint dof = 0; dof < ndof; ++dof) {
        uint tri_num_indx  = (uint)std::ceil((-3. + std::sqrt(1. + 8. * (dof + 1))) / 2.);
        uint lower_tri_num = (tri_num_indx + 1) * tri_num_indx / 2;

        uint p = dof - lower_tri_num;
        uint q = tri_num_indx - p;

        row(phi, dof) = transpose(this->ComputePhi(p, q, n1, n2));
    }

    return phi;
}

std::array<DynMatrix<double>, 2> Dubiner_2D::GetDPhi(const uint p, const std::vector<Point<2>>& points) {
    uint ndof = (p + 1) * (p + 2) / 2;
    uint npt  = points.size();

    std::array<DynMatrix<double>, 2> dphi;

    DynMatrix<double> dphi_dz1(ndof, npt);
    DynMatrix<double> dphi_dz2(ndof, npt);

    std::vector<double> n1(npt);
    std::vector<double> n2(npt);

    for (uint pt = 0; pt < npt; ++pt) {
        n1[pt] = 2 * (1 + points[pt][LocalCoordTri::z1]) / (1 - points[pt][LocalCoordTri::z2]) - 1;
        n2[pt] = points[pt][LocalCoordTri::z2];
    }

    for (uint dof = 0; dof < ndof; ++dof) {
        uint tri_num_indx  = (uint)std::ceil((-3. + std::sqrt(1. + 8. * (dof + 1))) / 2.);
        uint lower_tri_num = (tri_num_indx + 1) * tri_num_indx / 2;

        uint p = dof - lower_tri_num;
        uint q = tri_num_indx - p;

        row(dphi_dz1, dof) = transpose(this->ComputeDPhiDZ1(p, q, n1, n2));
        row(dphi_dz2, dof) = transpose(this->ComputeDPhiDZ2(p, q, n1, n2));
    }

    dphi[LocalCoordTri::z1] = dphi_dz1;
    dphi[LocalCoordTri::z2] = dphi_dz2;

    return dphi;
}

DynMatrix<double> Dubiner_2D::GetMinv(const uint p) {
    uint ndof = (p + 1) * (p + 2) / 2;

    DynMatrix<double> m_inv(ndof, ndof);

    set_constant(m_inv, 0.0);

    for (uint dof = 0; dof < ndof; ++dof) {
        uint tri_num_indx  = (uint)std::ceil((-3. + std::sqrt(1. + 8. * (dof + 1))) / 2.);
        uint lower_tri_num = (tri_num_indx + 1) * tri_num_indx / 2;

        uint p = dof - lower_tri_num;
        uint q = tri_num_indx - p;

        m_inv(dof, dof) = ((2 * p + 1) * (p + q + 1) / 2.0);
    }

    return m_inv;
}

DynVector<double> Dubiner_2D::ComputePhi(const uint p,
                                         const uint q,
                                         const std::vector<double>& n1,
                                         const std::vector<double>& n2) {
    assert(n1.size() == n2.size());

    uint npt = n1.size();

    DynVector<double> phi(npt);

    DynVector<double> psi_p  = jacobi_polynomial(p, 0, 0, n1);
    DynVector<double> psi_pq = jacobi_polynomial(q, 2 * p + 1, 0, n2);

    for (uint pt = 0; pt < npt; ++pt) {
        phi[pt] = psi_p[pt] * std::pow((1 - n2[pt]) / 2, (int)p) * psi_pq[pt];

        if (std::isnan(n1[pt])) {  // value of Dubiner polynomial at singular point (-1,1)
            if (p == 0) {
                phi[pt] = q + 1;
            } else {
                phi[pt] = 0;
            }
        }
    }

    return phi;
}

DynVector<double> Dubiner_2D::ComputeDPhiDZ1(const uint p,
                                             const uint q,
                                             const std::vector<double>& n1,
                                             const std::vector<double>& n2) {
    assert(n1.size() == n2.size());

    uint npt = n1.size();

    DynVector<double> dphi_dz1(npt);

    DynVector<double> psi_pq     = jacobi_polynomial(q, 2 * p + 1, 0, n2);
    DynVector<double> dpsi_p_dn1 = jacobi_polynomial_derivative(p, 0, 0, n1);

    for (uint pt = 0; pt < npt; ++pt) {
        dphi_dz1[pt] = (2.0 / (1.0 - n2[pt])) * dpsi_p_dn1[pt] * pow((1.0 - n2[pt]) / 2.0, (int)p) * psi_pq[pt];

        if (std::isnan(n1[pt])) {  // value of Dubiner polynomial derivatives at singular point (-1,1)
            if (p == 1) {
                dphi_dz1[pt] = this->ComputeSingularDPhiDZ1(q)[1];
            } else {
                dphi_dz1[pt] = 0.0;
            }
        }
    }

    return dphi_dz1;
}

DynVector<double> Dubiner_2D::ComputeDPhiDZ2(const uint p,
                                             const uint q,
                                             const std::vector<double>& n1,
                                             const std::vector<double>& n2) {
    assert(n1.size() == n2.size());

    uint npt = n1.size();

    DynVector<double> dphi_dz2(npt);

    DynVector<double> psi_p  = jacobi_polynomial(p, 0, 0, n1);
    DynVector<double> psi_pq = jacobi_polynomial(q, 2 * p + 1, 0, n2);

    DynVector<double> dpsi_p_dn1  = jacobi_polynomial_derivative(p, 0, 0, n1);
    DynVector<double> dpsi_pq_dn2 = jacobi_polynomial_derivative(q, 2 * p + 1, 0, n2);

    for (uint pt = 0; pt < npt; ++pt) {
        dphi_dz2[pt] =
            ((1 + n1[pt]) / (1.0 - n2[pt])) * dpsi_p_dn1[pt] * pow((1.0 - n2[pt]) / 2.0, (int)p) * psi_pq[pt] +
            psi_p[pt] * (pow((1.0 - n2[pt]) / 2.0, (int)p) * dpsi_pq_dn2[pt] -
                         (p / 2.0) * pow((1.0 - n2[pt]) / 2.0, (int)(p - 1)) * psi_pq[pt]);

        if (std::isnan(n1[pt])) {  // value of Dubiner polynomial derivatives at singular point (-1,1)
            if (p == 0) {
                dphi_dz2[pt] = 3 * this->ComputeSingularDPhiDZ2(q)[1];
            } else if (p == 1) {
                dphi_dz2[pt] = this->ComputeSingularDPhiDZ2(q + 1)[1];
            } else if (p > 1) {
                dphi_dz2[pt] = 0.0;
            }
        }
    }

    return dphi_dz2;
}

std::array<double, 2> Dubiner_2D::ComputeSingularDPhiDZ1(const uint q) {
    assert(q >= 0);

    std::array<double, 2> dphi_data;

    if (q == 0) {
        dphi_data[0] = 3.;
        dphi_data[1] = 1.;
    } else {
        std::array<double, 2> temp_dphi_data = this->ComputeSingularDPhiDZ1(q - 1);

        dphi_data[0] = temp_dphi_data[0] + q + 2;
        dphi_data[1] = temp_dphi_data[0] + temp_dphi_data[1];
    }

    return dphi_data;
}

std::array<double, 2> Dubiner_2D::ComputeSingularDPhiDZ2(const uint q) {
    assert(q >= 0);

    std::array<double, 2> dphi_data;

    if (q == 0) {
        dphi_data[0] = 0.5;
        dphi_data[1] = 0;
    } else {
        std::array<double, 2> temp_dphi_data = this->ComputeSingularDPhiDZ2(q - 1);

        dphi_data[0] = temp_dphi_data[0] + 0.5 * (q + 1);
        dphi_data[1] = temp_dphi_data[0] + temp_dphi_data[1];
    }

    return dphi_data;
}

DynMatrix<double> Dubiner_2D::GetBasisLinearT(const uint p) {
    uint ndof     = (p + 1) * (p + 2) / 2;
    uint ndof_lin = 3;

    DynMatrix<double> T_basis_linear(ndof, ndof_lin);

    set_constant(T_basis_linear, 0.0);

    T_basis_linear(0, 0) = 1.0;
    T_basis_linear(1, 0) = -1.0;
    T_basis_linear(2, 0) = -1.0;

    T_basis_linear(0, 1) = 1.0;
    T_basis_linear(1, 1) = -1.0;
    T_basis_linear(2, 1) = 1.0;

    T_basis_linear(0, 2) = 1.0;
    T_basis_linear(1, 2) = 2.0;
    T_basis_linear(2, 2) = 0.0;

    return T_basis_linear;
}

DynMatrix<double> Dubiner_2D::GetLinearBasisT(const uint p) {
    uint ndof     = (p + 1) * (p + 2) / 2;
    uint ndof_lin = 3;

    DynMatrix<double> T_linear_basis(ndof_lin, ndof);

    set_constant(T_linear_basis, 0.0);

    T_linear_basis(0, 0) = 1.0 / 3.0;
    T_linear_basis(1, 0) = 1.0 / 3.0;
    T_linear_basis(2, 0) = 1.0 / 3.0;

    T_linear_basis(0, 1) = -1.0 / 6.0;
    T_linear_basis(1, 1) = -1.0 / 6.0;
    T_linear_basis(2, 1) = 1.0 / 3.0;

    T_linear_basis(0, 2) = -1.0 / 2.0;
    T_linear_basis(1, 2) = 1.0 / 2.0;
    T_linear_basis(2, 2) = 0.0;

    return T_linear_basis;
}

DynMatrix<double> Dubiner_2D::ProjectBasisToLinear(const DynMatrix<double>& u) {
    uint nvar = rows(u);

    DynMatrix<double> u_lin(nvar, 3);

    column(u_lin, 0) = column(u, 0) - column(u, 1) - column(u, 2);
    column(u_lin, 1) = column(u, 0) - column(u, 1) + column(u, 2);
    column(u_lin, 2) = column(u, 0) + 2.0 * column(u, 1);

    return u_lin;
}

DynMatrix<double> Dubiner_2D::ProjectLinearToBasis(const uint ndof, const DynMatrix<double>& u_lin) {
    uint nvar = rows(u_lin);

    DynMatrix<double> u(nvar, ndof);

    set_constant(u, 0.0);

    column(u, 0) = (column(u_lin, 0) + column(u_lin, 1) + column(u_lin, 2)) / 3.0;
    column(u, 1) = (-column(u_lin, 0) - column(u_lin, 1) + 2.0 * column(u_lin, 2)) / 6.0;
    column(u, 2) = (-column(u_lin, 0) + column(u_lin, 1)) / 2.0;

    return u;
}
}