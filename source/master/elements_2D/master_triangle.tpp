#include "../master_elements_2D.hpp"

namespace Master {
template <typename BasisType, typename IntegrationType>
Triangle<BasisType, IntegrationType>::Triangle(const uint p) : Master<2>(p) {
    this->nvrtx  = 3;
    this->nbound = 3;

    this->integration_rule = this->integration.GetRule(2 * this->p);

    this->ndof = (p + 1) * (p + 2) / 2;
    this->ngp  = this->integration_rule.first.size();

    this->T_basis_linear = this->basis.GetBasisLinearT(p);
    this->T_linear_basis = this->basis.GetLinearBasisT(p);

    this->chi_baryctr.resize(this->nvrtx);

    this->chi_baryctr[0] = 1.0 / 3.0;
    this->chi_baryctr[1] = 1.0 / 3.0;
    this->chi_baryctr[2] = 1.0 / 3.0;

    this->chi_midpts.resize(this->nvrtx, this->nbound);

    this->chi_midpts(0, 0) = 0.0;
    this->chi_midpts(0, 1) = 1.0 / 2.0;
    this->chi_midpts(0, 2) = 1.0 / 2.0;
    this->chi_midpts(1, 0) = 1.0 / 2.0;
    this->chi_midpts(1, 1) = 0.0;
    this->chi_midpts(1, 2) = 1.0 / 2.0;
    this->chi_midpts(2, 0) = 1.0 / 2.0;
    this->chi_midpts(2, 1) = 1.0 / 2.0;
    this->chi_midpts(2, 2) = 0.0;

    this->chi_gp.resize(this->nvrtx, this->ngp);
    this->dchi_gp[LocalCoordTri::z1].resize(this->nvrtx, this->ngp);
    this->dchi_gp[LocalCoordTri::z2].resize(this->nvrtx, this->ngp);

    for (uint gp = 0; gp < this->ngp; ++gp) {
        this->chi_gp(0, gp) = -(this->integration_rule.second[gp][LocalCoordTri::z1] +
                                this->integration_rule.second[gp][LocalCoordTri::z2]) /
                              2.0;
        this->chi_gp(1, gp) = (1 + this->integration_rule.second[gp][LocalCoordTri::z1]) / 2.0;
        this->chi_gp(2, gp) = (1 + this->integration_rule.second[gp][LocalCoordTri::z2]) / 2.0;

        this->dchi_gp[LocalCoordTri::z1](0, gp) = -0.5;
        this->dchi_gp[LocalCoordTri::z2](0, gp) = -0.5;

        this->dchi_gp[LocalCoordTri::z1](1, gp) = 0.5;
        this->dchi_gp[LocalCoordTri::z2](1, gp) = 0.0;

        this->dchi_gp[LocalCoordTri::z1](2, gp) = 0.0;
        this->dchi_gp[LocalCoordTri::z2](2, gp) = 0.5;
    }

    this->phi_gp  = this->basis.GetPhi(this->p, this->integration_rule.second);
    this->dphi_gp = this->basis.GetDPhi(this->p, this->integration_rule.second);

    std::vector<Point<2>> z_postprocessor_cell = this->VTKPostCell();
    this->phi_postprocessor_cell               = this->basis.GetPhi(this->p, z_postprocessor_cell);

    std::vector<Point<2>> z_postprocessor_point = this->VTKPostPoint();
    this->phi_postprocessor_point               = this->basis.GetPhi(this->p, z_postprocessor_point);

    this->int_phi_fact = transpose(this->phi_gp);
    for (uint dof = 0; dof < this->ndof; ++dof) {
        for (uint gp = 0; gp < this->ngp; ++gp) {
            this->int_phi_fact(gp, dof) *= this->integration_rule.first[gp];
        }
    }

    for (uint dir = 0; dir < 2; ++dir) {
        this->int_dphi_fact[dir] = transpose(this->dphi_gp[dir]);
        for (uint dof = 0; dof < this->ndof; ++dof) {
            for (uint gp = 0; gp < this->ngp; ++gp) {
                this->int_dphi_fact[dir](gp, dof) *= this->integration_rule.first[gp];
            }
        }
    }

    this->m_inv = this->basis.GetMinv(this->p);
}

template <typename BasisType, typename IntegrationType>
std::vector<Point<2>> Triangle<BasisType, IntegrationType>::BoundaryToMasterCoordinates(
    const uint bound_id,
    const std::vector<Point<1>>& z_boundary) const {
    // *** //
    uint ngp = z_boundary.size();

    std::vector<Point<2>> z_master(ngp);

    if (bound_id == 0) {
        for (uint gp = 0; gp < ngp; ++gp) {
            assert(std::abs(z_boundary[gp][LocalCoordTri::z1]) < 1 + 100 * std::numeric_limits<double>::epsilon());

            z_master[gp][LocalCoordTri::z1] = -z_boundary[gp][LocalCoordTri::z1];
            z_master[gp][LocalCoordTri::z2] = z_boundary[gp][LocalCoordTri::z1];
        }
    } else if (bound_id == 1) {
        for (uint gp = 0; gp < ngp; ++gp) {
            assert(std::abs(z_boundary[gp][LocalCoordTri::z1]) < 1 + 100 * std::numeric_limits<double>::epsilon());

            z_master[gp][LocalCoordTri::z1] = -1.0;
            z_master[gp][LocalCoordTri::z2] = -z_boundary[gp][LocalCoordTri::z1];
        }
    } else if (bound_id == 2) {
        for (uint gp = 0; gp < ngp; ++gp) {
            assert(std::abs(z_boundary[gp][LocalCoordTri::z1]) < 1 + 100 * std::numeric_limits<double>::epsilon());

            z_master[gp][LocalCoordTri::z1] = z_boundary[gp][LocalCoordTri::z1];
            z_master[gp][LocalCoordTri::z2] = -1.0;
        }
    }

    return z_master;
}

template <typename BasisType, typename IntegrationType>
template <typename InputArrayType>
decltype(auto) Triangle<BasisType, IntegrationType>::ProjectBasisToLinear(const InputArrayType& u) const {
    return u * this->T_basis_linear;
}

template <typename BasisType, typename IntegrationType>
template <typename InputArrayType>
decltype(auto) Triangle<BasisType, IntegrationType>::ProjectLinearToBasis(const InputArrayType& u_lin) const {
    return u_lin * this->T_linear_basis;
}

template <typename BasisType, typename IntegrationType>
template <typename InputArrayType>
inline decltype(auto) Triangle<BasisType, IntegrationType>::ComputeLinearUbaryctr(const InputArrayType& u_lin) const {
    return u_lin * this->chi_baryctr;
}

template <typename BasisType, typename IntegrationType>
template <typename InputArrayType>
inline decltype(auto) Triangle<BasisType, IntegrationType>::ComputeLinearUmidpts(const InputArrayType& u_lin) const {
    return u_lin * this->chi_midpts;
}

template <typename BasisType, typename IntegrationType>
template <typename InputArrayType>
inline decltype(auto) Triangle<BasisType, IntegrationType>::ComputeLinearUvrtx(const InputArrayType& u_lin) const {
    return u_lin;
}

template <typename BasisType, typename IntegrationType>
std::vector<Point<2>> Triangle<BasisType, IntegrationType>::VTKPostCell() const {
    std::vector<Point<2>> z_postprocessor_cell(N_DIV * N_DIV);

    double dz = 2.0 / N_DIV;

    uint n_pt = 0;
    for (uint i = 0; i < N_DIV; ++i) {
        for (uint j = 0; j < N_DIV - i; ++j) {
            z_postprocessor_cell[n_pt][LocalCoordTri::z1] = -1.0 + dz * j + dz / 3.0;  // CENTROID
            z_postprocessor_cell[n_pt][LocalCoordTri::z2] = -1.0 + dz * i + dz / 3.0;
            n_pt++;
        }
    }
    for (uint i = 1; i < N_DIV; ++i) {
        for (uint j = 0; j < N_DIV - i; ++j) {
            z_postprocessor_cell[n_pt][LocalCoordTri::z1] = -1.0 + dz * j + 2 * dz / 3.0;  // CENTROID
            z_postprocessor_cell[n_pt][LocalCoordTri::z2] = -1.0 + dz * i - dz / 3.0;
            n_pt++;
        }
    }

    return z_postprocessor_cell;
}

template <typename BasisType, typename IntegrationType>
std::vector<Point<2>> Triangle<BasisType, IntegrationType>::VTKPostPoint() const {
    std::vector<Point<2>> z_postprocessor_point((N_DIV + 1) * (N_DIV + 2) / 2);

    double dz = 2.0 / N_DIV;

    uint n_pt = 0;
    for (uint i = 0; i <= N_DIV; ++i) {
        for (uint j = 0; j <= N_DIV - i; ++j) {
            z_postprocessor_point[n_pt][LocalCoordTri::z1] = -1.0 + dz * j;
            z_postprocessor_point[n_pt][LocalCoordTri::z2] = -1.0 + dz * i;
            n_pt++;
        }
    }

    return z_postprocessor_point;
}
}