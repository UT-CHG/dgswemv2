#ifndef CLASS_INTERFACE_HPP
#define CLASS_INTERFACE_HPP

namespace Geometry {
template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
class Interface {
  public:
    SpecializationType specialization;

    uint      bound_id_in;
    uint      bound_id_ex;
    DataType& data_in;
    DataType& data_ex;

    Array2D<double> surface_normal_in;
    Array2D<double> surface_normal_ex;

  private:
    Array2D<double> psi_gp_in;
    Array2D<double> psi_gp_ex;
    Array2D<double> phi_gp_in;
    Array2D<double> phi_gp_ex;

    std::vector<double> int_fact_in;
    std::vector<double> int_fact_ex;
    Array2D<double>     int_fact_phi_in;
    Array2D<double>     int_fact_phi_ex;

  public:
    Interface(const RawBoundary<dimension, DataType>& raw_boundary_in,
              const RawBoundary<dimension, DataType>& raw_boundary_ex,
              const SpecializationType&               specialization = SpecializationType());

    void ComputeUgpIN(const std::vector<double>& u, std::vector<double>& u_gp);
    void ComputeUgpEX(const std::vector<double>& u, std::vector<double>& u_gp);

    void ComputeNodalUgpIN(const std::vector<double>& u_nodal, std::vector<double>& u_nodal_gp);
    void ComputeNodalUgpEX(const std::vector<double>& u_nodal, std::vector<double>& u_nodal_gp);

    double IntegrationIN(const std::vector<double>& u_gp);
    double IntegrationEX(const std::vector<double>& u_gp);
    double IntegrationPhiIN(const uint dof, const std::vector<double>& u_gp);
    double IntegrationPhiEX(const uint dof, const std::vector<double>& u_gp);
};

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
Interface<dimension, IntegrationType, DataType, SpecializationType>::Interface(
    const RawBoundary<dimension, DataType>& raw_boundary_in,
    const RawBoundary<dimension, DataType>& raw_boundary_ex,
    const SpecializationType&               specialization)
    : specialization(std::move(specialization)),
      bound_id_in(raw_boundary_in.bound_id),
      bound_id_ex(raw_boundary_ex.bound_id),
      data_in(raw_boundary_in.data),
      data_ex(raw_boundary_ex.data) {
    // *** //
    uint p = std::max(raw_boundary_in.p, raw_boundary_ex.p);

    IntegrationType integration;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * p + 1);

    // transfrom gp to master coord in
    std::vector<Point<dimension + 1>> z_master =
        raw_boundary_in.master.BoundaryToMasterCoordinates(this->bound_id_in, integration_rule.second);

    // Compute factors to expand nodal values in
    this->psi_gp_in.resize(raw_boundary_in.shape.nodal_coordinates.size());

    std::vector<double> u_temp_in(raw_boundary_in.shape.nodal_coordinates.size());
    for (uint dof = 0; dof < raw_boundary_in.shape.nodal_coordinates.size(); dof++) {
        std::fill(u_temp_in.begin(), u_temp_in.end(), 0.0);
        u_temp_in[dof] = 1.0;

        this->psi_gp_in[dof] = raw_boundary_in.shape.InterpolateNodalValues(u_temp_in, z_master);
    }

    // Compute factors to expand modal values in
    this->phi_gp_in = raw_boundary_in.basis.GetPhi(raw_boundary_in.p, z_master);

    // transfrom gp to master coord ex
    z_master = raw_boundary_ex.master.BoundaryToMasterCoordinates(this->bound_id_ex, integration_rule.second);

    // Compute factors to expand nodal values ex
    this->psi_gp_ex.resize(raw_boundary_ex.shape.nodal_coordinates.size());

    std::vector<double> u_temp_ex(raw_boundary_ex.shape.nodal_coordinates.size());
    for (uint dof = 0; dof < raw_boundary_ex.shape.nodal_coordinates.size(); dof++) {
        std::fill(u_temp_ex.begin(), u_temp_ex.end(), 0.0);
        u_temp_ex[dof] = 1.0;

        this->psi_gp_ex[dof] = raw_boundary_ex.shape.InterpolateNodalValues(u_temp_ex, z_master);
    }

    // Compute factors to expand modal values ex
    this->phi_gp_ex = raw_boundary_ex.basis.GetPhi(raw_boundary_ex.p, z_master);

    std::vector<double> surface_J = raw_boundary_in.shape.GetSurfaceJ(this->bound_id_in, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_fact_in = integration_rule.first;
        for (uint gp = 0; gp < this->int_fact_in.size(); gp++) {
            this->int_fact_in[gp] *= surface_J[0];
        }

        this->int_fact_ex = integration_rule.first;
        for (uint gp = 0; gp < this->int_fact_ex.size(); gp++) {
            this->int_fact_ex[gp] *= surface_J[0];
        }

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

        this->surface_normal_in =
            Array2D<double>(integration_rule.first.size(),
                            *raw_boundary_in.shape.GetSurfaceNormal(this->bound_id_in, z_master).begin());

        this->surface_normal_ex = this->surface_normal_in;  // same dimensions

        uint ngp   = this->surface_normal_ex.size();
        uint gp_ex = 0;
        for (uint gp = 0; gp < this->surface_normal_ex.size(); gp++) {
            gp_ex = ngp - gp - 1;
            for (uint dir = 0; dir < this->surface_normal_ex[gp].size(); dir++) {
                this->surface_normal_ex[gp_ex][dir] = -this->surface_normal_in[gp][dir];
            }
        }
    }

    this->data_in.set_ngp_boundary(this->bound_id_in, integration_rule.first.size());
    this->data_ex.set_ngp_boundary(this->bound_id_ex, integration_rule.first.size());
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline void Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeUgpIN(
    const std::vector<double>& u,
    std::vector<double>&       u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->phi_gp_in[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline void Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeNodalUgpIN(
    const std::vector<double>& u_nodal,
    std::vector<double>&       u_nodal_gp) {
    std::fill(u_nodal_gp.begin(), u_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_nodal_gp.size(); gp++) {
            u_nodal_gp[gp] += u_nodal[dof] * this->psi_gp_in[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline double Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationIN(
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_in[gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline double Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiIN(
    const uint                 dof,
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_phi_in[dof][gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline void Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeUgpEX(
    const std::vector<double>& u,
    std::vector<double>&       u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->phi_gp_ex[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline void Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeNodalUgpEX(
    const std::vector<double>& u_nodal,
    std::vector<double>&       u_nodal_gp) {
    std::fill(u_nodal_gp.begin(), u_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_nodal_gp.size(); gp++) {
            u_nodal_gp[gp] += u_nodal[dof] * this->psi_gp_ex[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline double Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationEX(
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_ex[gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline double Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiEX(
    const uint                 dof,
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_phi_ex[dof][gp];
    }

    return integral;
}
}

#endif
