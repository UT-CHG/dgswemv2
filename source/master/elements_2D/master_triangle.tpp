#include "../master_elements_2D.hpp"

namespace Master {
template <typename BasisType, typename IntegrationType>
Triangle<BasisType, IntegrationType>::Triangle(const uint p) : Master<2>(p) {
    this->nvrtx  = 3;
    this->nbound = 3;

    this->integration_rule = this->integration.GetRule(2 * this->p);

    this->ndof = (p + 1) * (p + 2) / 2;
    this->ngp  = this->integration_rule.first.size();

    this->chi_gp.resize(this->ngp, 3);
    this->dchi_gp[LocalCoordTri::z1].resize(this->ngp, 3);
    this->dchi_gp[LocalCoordTri::z2].resize(this->ngp, 3);

    for (uint gp = 0; gp < this->ngp; gp++) {
        this->chi_gp(gp, 0) = -(this->integration_rule.second[gp][LocalCoordTri::z1] +
                                this->integration_rule.second[gp][LocalCoordTri::z2]) /
                              2.0;
        this->chi_gp(gp, 1) = (1 + this->integration_rule.second[gp][LocalCoordTri::z1]) / 2.0;
        this->chi_gp(gp, 2) = (1 + this->integration_rule.second[gp][LocalCoordTri::z2]) / 2.0;

        this->dchi_gp[LocalCoordTri::z1](gp, 0) = -0.5;
        this->dchi_gp[LocalCoordTri::z2](gp, 0) = -0.5;

        this->dchi_gp[LocalCoordTri::z1](gp, 1) = 0.5;
        this->dchi_gp[LocalCoordTri::z2](gp, 1) = 0.0;

        this->dchi_gp[LocalCoordTri::z1](gp, 2) = 0.0;
        this->dchi_gp[LocalCoordTri::z2](gp, 2) = 0.5;
    }

    this->phi_gp  = this->basis.GetPhi(this->p, this->integration_rule.second);
    this->dphi_gp = this->basis.GetDPhi(this->p, this->integration_rule.second);

    DynVector<Point<2>> z_postprocessor_cell = this->VTKPostCell();
    this->phi_postprocessor_cell             = this->basis.GetPhi(this->p, z_postprocessor_cell);

    DynVector<Point<2>> z_postprocessor_point = this->VTKPostPoint();
    this->phi_postprocessor_point             = this->basis.GetPhi(this->p, z_postprocessor_point);

    this->int_phi_fact = transpose(this->phi_gp);
    for (uint dof = 0; dof < this->ndof; dof++) {
        for (uint gp = 0; gp < this->ngp; gp++) {
            this->int_phi_fact(dof, gp) *= this->integration_rule.first[gp];
        }
    }

    for (uint dir = 0; dir < 2; dir++) {
        this->int_dphi_fact[dir] = transpose(this->dphi_gp[dir]);
        for (uint dof = 0; dof < this->ndof; dof++) {
            for (uint gp = 0; gp < this->ngp; gp++) {
                this->int_dphi_fact[dir](dof, gp) *= this->integration_rule.first[gp];
            }
        }
    }

    this->m_inv = this->basis.GetMinv(this->p);
}

template <typename BasisType, typename IntegrationType>
DynVector<Point<2>> Triangle<BasisType, IntegrationType>::BoundaryToMasterCoordinates(
    const uint bound_id,
    const DynVector<Point<1>>& z_boundary) {
    // *** //
    uint ngp = z_boundary.size();

    DynVector<Point<2>> z_master(ngp);

    if (bound_id == 0) {
        for (uint gp = 0; gp < ngp; gp++) {
            assert(std::abs(z_boundary[gp][LocalCoordTri::z1]) < 1 + 100 * std::numeric_limits<double>::epsilon());

            z_master[gp][LocalCoordTri::z1] = -z_boundary[gp][LocalCoordTri::z1];
            z_master[gp][LocalCoordTri::z2] = z_boundary[gp][LocalCoordTri::z1];
        }
    } else if (bound_id == 1) {
        for (uint gp = 0; gp < ngp; gp++) {
            assert(std::abs(z_boundary[gp][LocalCoordTri::z1]) < 1 + 100 * std::numeric_limits<double>::epsilon());

            z_master[gp][LocalCoordTri::z1] = -1.0;
            z_master[gp][LocalCoordTri::z2] = -z_boundary[gp][LocalCoordTri::z1];
        }
    } else if (bound_id == 2) {
        for (uint gp = 0; gp < ngp; gp++) {
            assert(std::abs(z_boundary[gp][LocalCoordTri::z1]) < 1 + 100 * std::numeric_limits<double>::epsilon());

            z_master[gp][LocalCoordTri::z1] = z_boundary[gp][LocalCoordTri::z1];
            z_master[gp][LocalCoordTri::z2] = -1.0;
        }
    }

    return z_master;
}

template <typename BasisType, typename IntegrationType>
template <typename T>
inline void Triangle<BasisType, IntegrationType>::ComputeLinearUbaryctr(const std::vector<T>& u_lin, T& u_lin_baryctr) {
    u_lin_baryctr = (u_lin[0] + u_lin[1] + u_lin[2]) / 3.0;
}

template <typename BasisType, typename IntegrationType>
template <typename T>
inline void Triangle<BasisType, IntegrationType>::ComputeLinearUmidpts(const std::vector<T>& u_lin,
                                                                       std::vector<T>& u_lin_midpts) {
    u_lin_midpts[0] = (u_lin[1] + u_lin[2]) / 2.0;
    u_lin_midpts[1] = (u_lin[2] + u_lin[0]) / 2.0;
    u_lin_midpts[2] = (u_lin[0] + u_lin[1]) / 2.0;
}

template <typename BasisType, typename IntegrationType>
template <typename T>
inline void Triangle<BasisType, IntegrationType>::ComputeLinearUvrtx(const std::vector<T>& u_lin,
                                                                     std::vector<T>& u_lin_vrtx) {
    u_lin_vrtx = u_lin;
}

template <typename BasisType, typename IntegrationType>
DynVector<Point<2>> Triangle<BasisType, IntegrationType>::VTKPostCell() {
    DynVector<Point<2>> z_postprocessor_cell(N_DIV * N_DIV);

    double dz = 2.0 / N_DIV;

    uint n_pt = 0;
    for (uint i = 0; i < N_DIV; i++) {
        for (uint j = 0; j < N_DIV - i; j++) {
            z_postprocessor_cell[n_pt][LocalCoordTri::z1] = -1.0 + dz * j + dz / 3.0;  // CENTROID
            z_postprocessor_cell[n_pt][LocalCoordTri::z2] = -1.0 + dz * i + dz / 3.0;
            n_pt++;
        }
    }
    for (uint i = 1; i < N_DIV; i++) {
        for (uint j = 0; j < N_DIV - i; j++) {
            z_postprocessor_cell[n_pt][LocalCoordTri::z1] = -1.0 + dz * j + 2 * dz / 3.0;  // CENTROID
            z_postprocessor_cell[n_pt][LocalCoordTri::z2] = -1.0 + dz * i - dz / 3.0;
            n_pt++;
        }
    }

    return z_postprocessor_cell;
}

template <typename BasisType, typename IntegrationType>
DynVector<Point<2>> Triangle<BasisType, IntegrationType>::VTKPostPoint() {
    DynVector<Point<2>> z_postprocessor_point((N_DIV + 1) * (N_DIV + 2) / 2);

    double dz = 2.0 / N_DIV;

    uint n_pt = 0;
    for (uint i = 0; i <= N_DIV; i++) {
        for (uint j = 0; j <= N_DIV - i; j++) {
            z_postprocessor_point[n_pt][LocalCoordTri::z1] = -1.0 + dz * j;
            z_postprocessor_point[n_pt][LocalCoordTri::z2] = -1.0 + dz * i;
            n_pt++;
        }
    }

    return z_postprocessor_point;
}
}