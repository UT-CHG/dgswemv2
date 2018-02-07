#ifndef CLASS_BOUNDARY_HPP
#define CLASS_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, typename IntegrationType, typename DataType, typename BoundaryType>
class Boundary {
  public:
    using BoundaryIntegrationType = IntegrationType;

  public:
    BoundaryType boundary_condition;

    DataType* data = nullptr;

    uint bound_id;
    std::vector<uint> node_ID;

    Master::Master<dimension + 1>& master;
    Shape::Shape<dimension + 1>& shape;

    Array2D<double> surface_normal;

  private:
    Array2D<double> psi_gp;
    Array2D<double> psi_bound_gp;
    Array2D<double> phi_gp;

    std::vector<double> int_fact;
    Array2D<double> int_phi_fact;

  public:
    Boundary(const RawBoundary<dimension, DataType>& raw_boundary,
             const BoundaryType& boundary_condition = BoundaryType());

    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);

    void ComputeNodalUgp(const std::vector<double>& u_nodal, std::vector<double>& u_nodal_gp);
    void ComputeBoundaryNodalUgp(const std::vector<double>& u_bound_nodal, std::vector<double>& u_bound_nodal_gp);

    double Integration(const std::vector<double>& u_gp);
    double IntegrationPhi(const uint dof, const std::vector<double>& u_gp);

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

template <uint dimension, typename IntegrationType, typename DataType, typename BoundaryType>
Boundary<dimension, IntegrationType, DataType, BoundaryType>::Boundary(
    const RawBoundary<dimension, DataType>& raw_boundary,
    const BoundaryType& boundary_condition)
    : boundary_condition(boundary_condition),
      data(&raw_boundary.data),
      bound_id(raw_boundary.bound_id),
      node_ID(raw_boundary.node_ID),
      master(raw_boundary.master),
      shape(raw_boundary.shape) {
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

        this->surface_normal = Array2D<double>(integration_rule.first.size(),
                                               *this->shape.GetSurfaceNormal(this->bound_id, z_master).begin());
    }

    this->data.set_ngp_boundary(this->bound_id, integration_rule.first.size());
}

template <uint dimension, typename IntegrationType, typename DataType, typename BoundaryType>
inline void Boundary<dimension, IntegrationType, DataType, BoundaryType>::ComputeUgp(const std::vector<double>& u,
                                                                                     std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->phi_gp[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename BoundaryType>
inline void Boundary<dimension, IntegrationType, DataType, BoundaryType>::ComputeNodalUgp(
    const std::vector<double>& u_nodal,
    std::vector<double>& u_nodal_gp) {
    std::fill(u_nodal_gp.begin(), u_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_nodal_gp.size(); gp++) {
            u_nodal_gp[gp] += u_nodal[dof] * this->psi_gp[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename BoundaryType>
inline void Boundary<dimension, IntegrationType, DataType, BoundaryType>::ComputeBoundaryNodalUgp(
    const std::vector<double>& u_bound_nodal,
    std::vector<double>& u_bound_nodal_gp) {
    std::fill(u_bound_nodal_gp.begin(), u_bound_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_bound_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_bound_nodal_gp.size(); gp++) {
            u_bound_nodal_gp[gp] += u_bound_nodal[dof] * this->psi_bound_gp[dof][gp];
        }
    }
}

template <uint dimension, typename IntegrationType, typename DataType, typename BoundaryType>
inline double Boundary<dimension, IntegrationType, DataType, BoundaryType>::Integration(
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact[gp];
    }

    return integral;
}

template <uint dimension, typename IntegrationType, typename DataType, typename BoundaryType>
inline double Boundary<dimension, IntegrationType, DataType, BoundaryType>::IntegrationPhi(
    const uint dof,
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_fact[dof][gp];
    }

    return integral;
}

#ifdef HAS_HPX
template <uint dimension, class IntegrationType, class DataType, class BoundaryType>
template<typename Archive>
void Boundary<dimension, IntegrationType, DataType, BoundaryType>::serialize(Archive& ar, unsigned) {
    ar & boundary_condition
       & bound_id
       & surface_normal
       & phi_gp
       & int_fact
       & int_fact_phi;
}
#endif
}
#endif
