#ifndef CLASS_INTERFACE_HPP
#define CLASS_INTERFACE_HPP

namespace Geometry {
template <uint dimension, class IntegrationType, class DataType>
class Interface {
  public:
    uint bound_id_in;
    uint bound_id_ex;
    DataType& data_in;
    DataType& data_ex;

    Array2D<double> surface_normal;

  private:
    Array2D<double> phi_gp_in;
    Array2D<double> phi_gp_ex;
    Array2D<double> int_fact_phi_in;
    Array2D<double> int_fact_phi_ex;

  public:
    Interface(const RawBoundary<dimension, DataType>& raw_boundary_in,
              const RawBoundary<dimension, DataType>& raw_boundary_ex);

    void ComputeUgpIN(const std::vector<double>& u, std::vector<double>& u_gp);
    void ComputeUgpEX(const std::vector<double>& u, std::vector<double>& u_gp);
    double IntegrationPhiIN(const uint dof, const std::vector<double>& u_gp);
    double IntegrationPhiEX(const uint dof, const std::vector<double>& u_gp);
};

template <uint dimension, class IntegrationType, class DataType>
Interface<dimension, IntegrationType, DataType>::Interface(const RawBoundary<dimension, DataType>& raw_boundary_in,
                                                           const RawBoundary<dimension, DataType>& raw_boundary_ex)
    : data_in(raw_boundary_in.data), data_ex(raw_boundary_ex.data) {
    uint p = std::max(raw_boundary_in.p, raw_boundary_ex.p);

    IntegrationType integration;
    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * p);

    std::vector<Point<dimension + 1>> z_master =
        raw_boundary_ex.master.BoundaryToMasterCoordinates(raw_boundary_ex.bound_id, integration_rule.second);
    this->phi_gp_ex = raw_boundary_ex.basis.GetPhi(raw_boundary_ex.p, z_master);

    z_master = raw_boundary_in.master.BoundaryToMasterCoordinates(raw_boundary_in.bound_id, integration_rule.second);
    this->phi_gp_in = raw_boundary_in.basis.GetPhi(raw_boundary_in.p, z_master);

    std::vector<double> surface_J = raw_boundary_in.shape.GetSurfaceJ(raw_boundary_in.bound_id, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_fact_phi_in = this->phi_gp_in;
        for (uint dof = 0; dof < this->int_fact_phi_in.size(); dof++) {
            for (uint gp = 0; gp < this->int_fact_phi_in[dof].size(); gp++) {
                this->int_fact_phi_in[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_fact_phi_ex = this->phi_gp_ex;
        for (uint dof = 0; dof < this->int_fact_phi_ex.size(); dof++) {
            for (uint gp = 0; gp < this->int_fact_phi_ex[dof].size(); gp++) {
                this->int_fact_phi_ex[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->surface_normal =
            Array2D<double>(integration_rule.first.size(),
                            *raw_boundary_in.shape.GetSurfaceNormal(raw_boundary_in.bound_id, z_master).begin());
    }

    this->bound_id_in = raw_boundary_in.bound_id;
    this->bound_id_ex = raw_boundary_ex.bound_id;
    this->data_in.set_ngp_boundary(raw_boundary_in.bound_id, integration_rule.first.size());
    this->data_ex.set_ngp_boundary(raw_boundary_ex.bound_id, integration_rule.first.size());
}

template <uint dimension, class IntegrationType, class DataType>
inline void Interface<dimension, IntegrationType, DataType>::ComputeUgpIN(const std::vector<double>& u,
                                                                          std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->phi_gp_in[dof][gp];
        }
    }
}

template <uint dimension, class IntegrationType, class DataType>
inline double Interface<dimension, IntegrationType, DataType>::IntegrationPhiIN(const uint dof,
                                                                                const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_phi_in[dof][gp];
    }

    return integral;
}

template <uint dimension, class IntegrationType, class DataType>
inline void Interface<dimension, IntegrationType, DataType>::ComputeUgpEX(const std::vector<double>& u,
                                                                          std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->phi_gp_ex[dof][gp];
        }
    }
}

template <uint dimension, class IntegrationType, class DataType>
inline double Interface<dimension, IntegrationType, DataType>::IntegrationPhiEX(const uint dof,
                                                                                const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_phi_ex[dof][gp];
    }

    return integral;
}
}

#endif