#include "../master_elements_2D.hpp"

namespace Master {
template <typename BasisType, typename IntegrationType>
Triangle<BasisType, IntegrationType>::Triangle(const uint p)
    : Master<2>(p) {
    this->nvrtx = 3;
    this->nbound = 3;

    this->integration_rule = this->integration.GetRule(2 * this->p);

    this->psi_gp = Array2D<double>(3, std::vector<double>(this->integration_rule.first.size()));

    for (uint gp = 0; gp < this->integration_rule.first.size(); gp++) {
        this->psi_gp[0][gp] = -(this->integration_rule.second[gp][LocalCoordTri::z1] +
                                this->integration_rule.second[gp][LocalCoordTri::z2]) /
                              2.0;
        this->psi_gp[1][gp] = (1 + this->integration_rule.second[gp][LocalCoordTri::z1]) / 2.0;
        this->psi_gp[2][gp] = (1 + this->integration_rule.second[gp][LocalCoordTri::z2]) / 2.0;
    }

    this->dpsi = Array2D<double>{{-0.5, -0.5}, {0.5, 0.0}, {0.0, 0.5}};

    this->phi_gp = this->basis.GetPhi(this->p, this->integration_rule.second);
    this->dphi_gp = this->basis.GetDPhi(this->p, this->integration_rule.second);

    std::vector<Point<2>> z_postprocessor_cell = this->VTKPostCell();
    this->phi_postprocessor_cell = this->basis.GetPhi(this->p, z_postprocessor_cell);

    std::vector<Point<2>> z_postprocessor_point = this->VTKPostPoint();
    this->phi_postprocessor_point = this->basis.GetPhi(this->p, z_postprocessor_point);

    this->int_fact_phi = this->phi_gp;
    for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {
        for (uint gp = 0; gp < this->int_fact_phi[dof].size(); gp++) {
            this->int_fact_phi[dof][gp] *= this->integration_rule.first[gp];
        }
    }

    this->int_fact_dphi = this->dphi_gp;
    for (uint dof = 0; dof < this->int_fact_dphi.size(); dof++) {
        for (uint dir = 0; dir < this->int_fact_dphi[dof].size(); dir++) {
            for (uint gp = 0; gp < this->int_fact_dphi[dof][dir].size(); gp++) {
                this->int_fact_dphi[dof][dir][gp] *= this->integration_rule.first[gp];
            }
        }
    }

    this->m_inv = this->basis.GetMinv(this->p);
}

template <typename BasisType, typename IntegrationType>
std::vector<Point<2>> Triangle<BasisType, IntegrationType>::BoundaryToMasterCoordinates(
    const uint bound_id,
    const std::vector<Point<1>>& z_boundary) {
    std::vector<Point<2>> z_master(z_boundary.size());

    if (bound_id == 0) {
        for (uint gp = 0; gp < z_master.size(); gp++) {
            assert(std::abs(z_boundary[gp][LocalCoordTri::z1]) < 1 + 100 * std::numeric_limits<double>::epsilon());

            z_master[gp][LocalCoordTri::z1] = -z_boundary[gp][LocalCoordTri::z1];
            z_master[gp][LocalCoordTri::z2] = z_boundary[gp][LocalCoordTri::z1];
        }
    } else if (bound_id == 1) {
        for (uint gp = 0; gp < z_master.size(); gp++) {
            assert(std::abs(z_boundary[gp][LocalCoordTri::z1]) < 1 + 100 * std::numeric_limits<double>::epsilon());

            z_master[gp][LocalCoordTri::z1] = -1;
            z_master[gp][LocalCoordTri::z2] = -z_boundary[gp][LocalCoordTri::z1];
        }
    } else if (bound_id == 2) {
        for (uint gp = 0; gp < z_master.size(); gp++) {
            assert(std::abs(z_boundary[gp][LocalCoordTri::z1]) < 1 + 100 * std::numeric_limits<double>::epsilon());

            z_master[gp][LocalCoordTri::z1] = z_boundary[gp][LocalCoordTri::z1];
            z_master[gp][LocalCoordTri::z2] = -1;
        }
    }

    return z_master;
}

template <typename BasisType, typename IntegrationType>
inline void Triangle<BasisType, IntegrationType>::ComputeLinearUbaryctr(const std::vector<double>& u_lin,
                                                                        double& u_lin_baryctr) {
    u_lin_baryctr = (u_lin[0] + u_lin[1] + u_lin[2]) / 3.0;
}

template <typename BasisType, typename IntegrationType>
inline void Triangle<BasisType, IntegrationType>::ComputeLinearUmidpts(const std::vector<double>& u_lin,
                                                                       std::vector<double>& u_lin_midpts) {
    u_lin_midpts[0] = (u_lin[1] + u_lin[2]) / 2.0;
    u_lin_midpts[1] = (u_lin[2] + u_lin[0]) / 2.0;
    u_lin_midpts[2] = (u_lin[0] + u_lin[1]) / 2.0;
}

template <typename BasisType, typename IntegrationType>
inline void Triangle<BasisType, IntegrationType>::ComputeLinearUvrtx(const std::vector<double>& u_lin,
                                                                     std::vector<double>& u_lin_vrtx) {
    u_lin_vrtx = u_lin;
}

template <typename BasisType, typename IntegrationType>
std::vector<Point<2>> Triangle<BasisType, IntegrationType>::VTKPostCell() {
    std::vector<Point<2>> z_postprocessor_cell(N_DIV * N_DIV);

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
std::vector<Point<2>> Triangle<BasisType, IntegrationType>::VTKPostPoint() {
    std::vector<Point<2>> z_postprocessor_point((N_DIV + 1) * (N_DIV + 2) / 2);

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