#include "../master_elements_2D.hpp"

namespace Master {
template <class basis_type, class integration_type>
Triangle<basis_type, integration_type>::Triangle(uint p)
    : Master<2>(p) {
    this->integration_rule = this->integration.GetRule(2 * this->p);

    this->phi_gp = this->basis.GetPhi(this->p, this->integration_rule.second);
    this->dphi_gp = this->basis.GetDPhi(this->p, this->integration_rule.second);

    std::vector<Point<2>> z_postprocessor_cell = this->VTKPostCell();
    this->phi_postprocessor_cell = this->basis.GetPhi(this->p, z_postprocessor_cell);

    std::vector<Point<2>> z_postprocessor_point = this->VTKPostPoint();
    this->phi_postprocessor_point = this->basis.GetPhi(this->p, z_postprocessor_point);

    this->int_fact_phi = this->phi_gp;
    for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {  // iterate through basis functions
        for (uint gp = 0; gp < this->int_fact_phi[dof].size(); gp++) {  // iterate through internal GPs
            this->int_fact_phi[dof][gp] *= this->integration_rule.first[gp];  // apply weight
        }
    }

    this->int_fact_dphi = this->dphi_gp;
    for (uint dof = 0; dof < this->int_fact_dphi.size(); dof++) {  // iterate through basis functions
        for (uint dir = 0; dir < this->int_fact_dphi[dof].size(); dir++) {  // iterate through differentiation
                                                                            // directions
            for (uint gp = 0; gp < this->int_fact_dphi[dof][dir].size(); gp++) {  // iterate through internal GPs
                this->int_fact_dphi[dof][dir][gp] *= this->integration_rule.first[gp];  // apply weight
            }
        }
    }

    this->m_inv = this->basis.GetMinv(this->p);
}

template <class basis_type, class integration_type>
std::vector<Point<2>> Triangle<basis_type, integration_type>::BoundaryToMasterCoordinates(
    uint boundary,
    const std::vector<Point<1>>& z_boundary) {
    std::vector<Point<2>> z_master(z_boundary.size());

    if (boundary == 0) {
        for (uint gp = 0; gp < z_master.size(); gp++) {
            z_master[gp][LocalCoordTri::z1] = -z_boundary[gp][LocalCoordTri::z1];
            z_master[gp][LocalCoordTri::z2] = z_boundary[gp][LocalCoordTri::z1];
        }
    } else if (boundary == 1) {
        for (uint gp = 0; gp < z_master.size(); gp++) {
            z_master[gp][LocalCoordTri::z1] = -1;
            z_master[gp][LocalCoordTri::z2] = -z_boundary[gp][LocalCoordTri::z1];
        }
    } else if (boundary == 2) {
        for (uint gp = 0; gp < z_master.size(); gp++) {
            z_master[gp][LocalCoordTri::z1] = z_boundary[gp][LocalCoordTri::z1];
            z_master[gp][LocalCoordTri::z2] = -1;
        }
    }

    return z_master;
}

template <class basis_type, class integration_type>
std::vector<Point<2>> Triangle<basis_type, integration_type>::VTKPostCell() {
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

template <class basis_type, class integration_type>
std::vector<Point<2>> Triangle<basis_type, integration_type>::VTKPostPoint() {
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