#ifndef CLASS_BOUNDARY_HPP
#define CLASS_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
class Boundary {
  public:
    ConditonType boundary_condition;

    DataType& data;

    uint bound_id;

    Array2D<double> surface_normal;

  private:
    Master::Master<dimension + 1>& master;
    Shape::Shape<dimension + 1>& shape;

    std::vector<uint> node_ID;

    Array2D<double> psi_gp;
    Array2D<double> psi_bound_gp;
    Array2D<double> phi_gp;

    std::vector<double> int_fact;
    Array2D<double> int_phi_fact;
    Array3D<double> int_phi_phi_fact;

  public:
    Boundary(RawBoundary<dimension, DataType>&& raw_boundary, ConditonType&& boundary_condition = ConditonType());

    Master::Master<dimension + 1>& GetMaster() { return this->master; }
    Shape::Shape<dimension + 1>& GetShape() { return this->shape; }

    std::vector<uint>& GetNodeID() { return this->node_ID; }

    template <typename T>
    void ComputeUgp(const std::vector<T>& u, std::vector<T>& u_gp);

    void ComputeNodalUgp(const std::vector<double>& u_nodal, std::vector<double>& u_nodal_gp);
    void ComputeBoundaryNodalUgp(const std::vector<double>& u_bound_nodal, std::vector<double>& u_bound_nodal_gp);

    template <typename T>
    T Integration(const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhi(const uint dof, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhiPhi(const uint dof_i, const uint dof_j, const std::vector<T>& u_gp);

  public:
    using BoundaryIntegrationType = IntegrationType;
};

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
Boundary<dimension, IntegrationType, DataType, ConditonType>::Boundary(RawBoundary<dimension, DataType>&& raw_boundary,
                                                                       ConditonType&& boundary_condition)
    : boundary_condition(boundary_condition),
      data(raw_boundary.data),
      bound_id(raw_boundary.bound_id),
      master(raw_boundary.master),
      shape(raw_boundary.shape),
      node_ID(std::move(raw_boundary.node_ID)) {
    // *** //
    IntegrationType integration;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule =
        integration.GetRule(2 * raw_boundary.p + 1);

    std::vector<Point<dimension + 1>> z_master =
        this->master.BoundaryToMasterCoordinates(this->bound_id, integration_rule.second);

    // Compute factors to expand nodal values
    this->psi_gp = this->shape.GetPsi(z_master);

    // Compute factors to expand boundary nodal values
    this->psi_bound_gp = this->shape.GetBoundaryPsi(this->bound_id, integration_rule.second);

    // Compute factors to expand modal values
    this->phi_gp = raw_boundary.basis.GetPhi(raw_boundary.p, z_master);

    std::vector<double> surface_J = this->shape.GetSurfaceJ(this->bound_id, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_fact = integration_rule.first;
        for (uint gp = 0; gp < this->int_fact.size(); gp++) {
            this->int_fact[gp] *= surface_J[0];
        }

        this->int_phi_fact = this->phi_gp;
        for (uint dof = 0; dof < this->int_phi_fact.size(); dof++) {
            for (uint gp = 0; gp < this->int_phi_fact[dof].size(); gp++) {
                this->int_phi_fact[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_phi_phi_fact.resize(this->phi_gp.size());
        for (uint dof_i = 0; dof_i < this->phi_gp.size(); dof_i++) {
            this->int_phi_phi_fact[dof_i] = this->int_phi_fact;
            for (uint dof_j = 0; dof_j < this->phi_gp.size(); dof_j++) {
                for (uint gp = 0; gp < this->int_phi_phi_fact[dof_i][dof_j].size(); gp++) {
                    this->int_phi_phi_fact[dof_i][dof_j][gp] *= this->phi_gp[dof_i][gp];
                }
            }
        }

        this->surface_normal = Array2D<double>(integration_rule.first.size(),
                                               *this->shape.GetSurfaceNormal(this->bound_id, z_master).begin());
    }

    this->data.set_ngp_boundary(this->bound_id, integration_rule.first.size());
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename T>
inline void Boundary<dimension, IntegrationType, DataType, ConditonType>::ComputeUgp(const std::vector<T>& u,
                                                                                     std::vector<T>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->phi_gp[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
inline void Boundary<dimension, IntegrationType, DataType, ConditonType>::ComputeNodalUgp(
    const std::vector<double>& u_nodal,
    std::vector<double>& u_nodal_gp) {
    std::fill(u_nodal_gp.begin(), u_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_nodal_gp.size(); gp++) {
            u_nodal_gp[gp] += u_nodal[dof] * this->psi_gp[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
inline void Boundary<dimension, IntegrationType, DataType, ConditonType>::ComputeBoundaryNodalUgp(
    const std::vector<double>& u_bound_nodal,
    std::vector<double>& u_bound_nodal_gp) {
    std::fill(u_bound_nodal_gp.begin(), u_bound_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_bound_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_bound_nodal_gp.size(); gp++) {
            u_bound_nodal_gp[gp] += u_bound_nodal[dof] * this->psi_bound_gp[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename T>
inline T Boundary<dimension, IntegrationType, DataType, ConditonType>::Integration(const std::vector<T>& u_gp) {
    T integral;

    integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact[gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename T>
inline T Boundary<dimension, IntegrationType, DataType, ConditonType>::IntegrationPhi(const uint dof,
                                                                                      const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_fact[dof][gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename T>
inline T Boundary<dimension, IntegrationType, DataType, ConditonType>::IntegrationPhiPhi(const uint dof_i,
                                                                                         const uint dof_j,
                                                                                         const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_phi_fact[dof_i][dof_j][gp];
    }

    return integral;
}
}

#endif
