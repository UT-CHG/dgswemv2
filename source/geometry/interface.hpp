#ifndef CLASS_INTERFACE_HPP
#define CLASS_INTERFACE_HPP

namespace Geometry {
template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
class Interface {
  public:
    SpecializationType specialization;

    DataType& data_in;
    DataType& data_ex;

    uint bound_id_in;
    uint bound_id_ex;

    Array2D<double> surface_normal_in;
    Array2D<double> surface_normal_ex;

  private:
    Master::Master<dimension + 1>& master_in;
    Master::Master<dimension + 1>& master_ex;

    Shape::Shape<dimension + 1>& shape_in;
    Shape::Shape<dimension + 1>& shape_ex;

    std::vector<uint> node_ID_in;
    std::vector<uint> node_ID_ex;

    Array2D<double> psi_gp_in;
    Array2D<double> psi_gp_ex;
    Array2D<double> psi_bound_gp_in;
    Array2D<double> psi_bound_gp_ex;
    Array2D<double> phi_gp_in;
    Array2D<double> phi_gp_ex;

    std::vector<double> int_fact_in;
    std::vector<double> int_fact_ex;
    Array2D<double> int_phi_fact_in;
    Array2D<double> int_phi_fact_ex;
    Array3D<double> int_phi_phi_fact_in;
    Array3D<double> int_phi_phi_fact_ex;

  public:
    Interface(RawBoundary<dimension, DataType>&& raw_boundary_in,
              RawBoundary<dimension, DataType>&& raw_boundary_ex,
              SpecializationType&& specialization = SpecializationType());

    Master::Master<dimension + 1>& GetMasterIN() { return this->master_in; }
    Master::Master<dimension + 1>& GetMasterEX() { return this->master_ex; }

    Shape::Shape<dimension + 1>& GetShapeIN() { return this->shape_in; }
    Shape::Shape<dimension + 1>& GetShapeEX() { return this->shape_ex; }

    std::vector<uint>& GetNodeIDIN() { return this->node_ID_in; }
    std::vector<uint>& GetNodeIDEX() { return this->node_ID_ex; }

    template <typename T>
    void ComputeUgpIN(const std::vector<T>& u, std::vector<T>& u_gp);
    template <typename T>
    void ComputeUgpEX(const std::vector<T>& u, std::vector<T>& u_gp);

    void ComputeNodalUgpIN(const std::vector<double>& u_nodal, std::vector<double>& u_nodal_gp);
    void ComputeNodalUgpEX(const std::vector<double>& u_nodal, std::vector<double>& u_nodal_gp);
    void ComputeBoundaryNodalUgpIN(const std::vector<double>& u_bound_nodal, std::vector<double>& u_bound_nodal_gp);
    void ComputeBoundaryNodalUgpEX(const std::vector<double>& u_bound_nodal, std::vector<double>& u_bound_nodal_gp);

    template <typename T>
    T IntegrationIN(const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationEX(const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhiIN(const uint dof, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhiEX(const uint dof, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhiPhiIN(const uint dof_i, const uint dof_j, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhiPhiEX(const uint dof_i, const uint dof_j, const std::vector<T>& u_gp);

  public:
    using InterfaceIntegrationType = IntegrationType;

  private:
    template <uint d, typename BT, typename EDT, typename InT>
    friend class EdgeInterface;
};

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
Interface<dimension, IntegrationType, DataType, SpecializationType>::Interface(
    RawBoundary<dimension, DataType>&& raw_boundary_in,
    RawBoundary<dimension, DataType>&& raw_boundary_ex,
    SpecializationType&& specialization)
    : specialization(specialization),
      data_in(raw_boundary_in.data),
      data_ex(raw_boundary_ex.data),
      bound_id_in(raw_boundary_in.bound_id),
      bound_id_ex(raw_boundary_ex.bound_id),
      master_in(raw_boundary_in.master),
      master_ex(raw_boundary_ex.master),
      shape_in(raw_boundary_in.shape),
      shape_ex(raw_boundary_ex.shape),
      node_ID_in(std::move(raw_boundary_in.node_ID)),
      node_ID_ex(std::move(raw_boundary_ex.node_ID)) {
    // *** //
    uint p = std::max(raw_boundary_in.p, raw_boundary_ex.p);

    IntegrationType integration;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * p + 1);

    // transfrom gp to master coord in
    std::vector<Point<dimension + 1>> z_master_in =
        this->master_in.BoundaryToMasterCoordinates(this->bound_id_in, integration_rule.second);

    // Compute factors to expand nodal values in
    this->psi_gp_in = this->shape_in.GetPsi(z_master_in);

    // Compute factors to expand boundary nodal values in
    this->psi_bound_gp_in = this->shape_in.GetBoundaryPsi(this->bound_id_in, integration_rule.second);

    // Compute factors to expand modal values in
    this->phi_gp_in = raw_boundary_in.basis.GetPhi(raw_boundary_in.p, z_master_in);

    // transfrom gp to master coord ex
    std::vector<Point<dimension + 1>> z_master_ex =
        this->master_ex.BoundaryToMasterCoordinates(this->bound_id_ex, integration_rule.second);

    // Compute factors to expand nodal values ex
    this->psi_gp_ex = this->shape_ex.GetPsi(z_master_ex);

    // Compute factors to expand boundary nodal values ex
    this->psi_bound_gp_ex = this->shape_ex.GetBoundaryPsi(this->bound_id_ex, integration_rule.second);

    // Compute factors to expand modal values ex
    this->phi_gp_ex = raw_boundary_ex.basis.GetPhi(raw_boundary_ex.p, z_master_ex);

    std::vector<double> surface_J = this->shape_in.GetSurfaceJ(this->bound_id_in, z_master_in);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_fact_in = integration_rule.first;
        for (uint gp = 0; gp < this->int_fact_in.size(); gp++) {
            this->int_fact_in[gp] *= surface_J[0];
        }

        this->int_fact_ex = integration_rule.first;
        for (uint gp = 0; gp < this->int_fact_ex.size(); gp++) {
            this->int_fact_ex[gp] *= surface_J[0];
        }

        this->int_phi_fact_in = this->phi_gp_in;
        for (uint dof = 0; dof < this->int_phi_fact_in.size(); dof++) {
            for (uint gp = 0; gp < this->int_phi_fact_in[dof].size(); gp++) {
                this->int_phi_fact_in[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_phi_fact_ex = this->phi_gp_ex;
        for (uint dof = 0; dof < this->int_phi_fact_ex.size(); dof++) {
            for (uint gp = 0; gp < this->int_phi_fact_ex[dof].size(); gp++) {
                this->int_phi_fact_ex[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_phi_phi_fact_in.resize(this->phi_gp_in.size());
        for (uint dof_i = 0; dof_i < this->phi_gp_in.size(); dof_i++) {
            this->int_phi_phi_fact_in[dof_i] = this->int_phi_fact_in;
            for (uint dof_j = 0; dof_j < this->phi_gp_in.size(); dof_j++) {
                for (uint gp = 0; gp < this->int_phi_phi_fact_in[dof_i][dof_j].size(); gp++) {
                    this->int_phi_phi_fact_in[dof_i][dof_j][gp] *= this->phi_gp_in[dof_i][gp];
                }
            }
        }

        this->int_phi_phi_fact_ex.resize(this->phi_gp_ex.size());
        for (uint dof_i = 0; dof_i < this->phi_gp_ex.size(); dof_i++) {
            this->int_phi_phi_fact_ex[dof_i] = this->int_phi_fact_ex;
            for (uint dof_j = 0; dof_j < this->phi_gp_ex.size(); dof_j++) {
                for (uint gp = 0; gp < this->int_phi_phi_fact_ex[dof_i][dof_j].size(); gp++) {
                    this->int_phi_phi_fact_ex[dof_i][dof_j][gp] *= this->phi_gp_ex[dof_i][gp];
                }
            }
        }

        this->surface_normal_in = Array2D<double>(
            integration_rule.first.size(), *this->shape_in.GetSurfaceNormal(this->bound_id_in, z_master_in).begin());

        this->surface_normal_ex = this->surface_normal_in;  // same dimensions

        uint ngp = this->surface_normal_ex.size();
        uint gp_ex;
        for (uint gp = 0; gp < this->surface_normal_ex.size(); gp++) {
            gp_ex = ngp - gp - 1;
            for (uint dir = 0; dir < this->surface_normal_ex[gp].size(); dir++) {
                this->surface_normal_ex[gp_ex][dir] = -this->surface_normal_in[gp][dir];  // but opposite direction
            }
        }
    }

    this->data_in.set_ngp_boundary(this->bound_id_in, integration_rule.first.size());
    this->data_ex.set_ngp_boundary(this->bound_id_ex, integration_rule.first.size());
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename T>
inline void Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeUgpIN(const std::vector<T>& u,
                                                                                              std::vector<T>& u_gp) {
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
    std::vector<double>& u_nodal_gp) {
    std::fill(u_nodal_gp.begin(), u_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_nodal_gp.size(); gp++) {
            u_nodal_gp[gp] += u_nodal[dof] * this->psi_gp_in[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline void Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeBoundaryNodalUgpIN(
    const std::vector<double>& u_bound_nodal,
    std::vector<double>& u_bound_nodal_gp) {
    std::fill(u_bound_nodal_gp.begin(), u_bound_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_bound_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_bound_nodal_gp.size(); gp++) {
            u_bound_nodal_gp[gp] += u_bound_nodal[dof] * this->psi_bound_gp_in[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename T>
inline T Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationIN(
    const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_in[gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename T>
inline T Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiIN(
    const uint dof,
    const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_fact_in[dof][gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename T>
inline T Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiPhiIN(
    const uint dof_i,
    const uint dof_j,
    const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_phi_fact_in[dof_i][dof_j][gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename T>
inline void Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeUgpEX(const std::vector<T>& u,
                                                                                              std::vector<T>& u_gp) {
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
    std::vector<double>& u_nodal_gp) {
    std::fill(u_nodal_gp.begin(), u_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_nodal_gp.size(); gp++) {
            u_nodal_gp[gp] += u_nodal[dof] * this->psi_gp_ex[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
inline void Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeBoundaryNodalUgpEX(
    const std::vector<double>& u_bound_nodal,
    std::vector<double>& u_bound_nodal_gp) {
    std::fill(u_bound_nodal_gp.begin(), u_bound_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_bound_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_bound_nodal_gp.size(); gp++) {
            u_bound_nodal_gp[gp] += u_bound_nodal[dof] * this->psi_bound_gp_ex[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename T>
inline T Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationEX(
    const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_ex[gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename T>
inline T Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiEX(
    const uint dof,
    const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_fact_ex[dof][gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename T>
inline T Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiPhiEX(
    const uint dof_i,
    const uint dof_j,
    const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_phi_fact_ex[dof_i][dof_j][gp];
    }

    return integral;
}
}

#endif
