#ifndef CLASS_BOUNDARY_HPP
#define CLASS_BOUNDARY_HPP

#include "general_definitions.hpp"

namespace Geometry {
template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
class Boundary {
  public:
    ConditonType boundary_condition;

    DataType& data;

    uint bound_id;

    std::array<DynRowVector<double>, dimension + 1> surface_normal;

  private:
    Master::Master<dimension + 1>& master;
    Shape::Shape<dimension + 1>& shape;

    std::vector<uint> node_ID;

    std::vector<Point<dimension + 1>> gp_global_coordinates;

    DynMatrix<double> psi_gp;
    DynMatrix<double> psi_bound_gp;
    DynMatrix<double> phi_gp;

    DynVector<double> int_fact;
    DynMatrix<double> int_phi_fact;
    DynMatrix<double> int_phi_phi_fact;

  public:
    template <typename... Args>
    Boundary(RawBoundary<dimension, DataType>&& raw_boundary, Args&&... args);

    Master::Master<dimension + 1>& GetMaster() { return this->master; }
    Shape::Shape<dimension + 1>& GetShape() { return this->shape; }

    std::vector<uint>& GetNodeID() { return this->node_ID; }

    template <typename F>
    DynMatrix<double> ComputeFgp(const F& f);
    template <typename InputArrayType>
    decltype(auto) ComputeUgp(const InputArrayType& u);

    template <typename InputArrayType>
    decltype(auto) ComputeNodalUgp(const InputArrayType& u_nodal);
    template <typename InputArrayType>
    decltype(auto) ComputeBoundaryNodalUgp(const InputArrayType& u_bound_nodal);

    template <typename InputArrayType>
    decltype(auto) Integration(const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhi(const uint dof, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhi(const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiPhi(const uint dof_i, const uint dof_j, const InputArrayType& u_gp);

  public:
    using BoundaryIntegrationType = IntegrationType;

  private:
    template <uint d, typename BT, typename EDT, typename BdT>
    friend class EdgeBoundary;
};

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename... Args>
Boundary<dimension, IntegrationType, DataType, ConditonType>::Boundary(RawBoundary<dimension, DataType>&& raw_boundary,
                                                                       Args&&... args)
    : boundary_condition(std::forward<Args>(args)...),
      data(raw_boundary.data),
      bound_id(raw_boundary.bound_id),
      master(raw_boundary.master),
      shape(raw_boundary.shape),
      node_ID(std::move(raw_boundary.node_ID)) {
    // *** //
    IntegrationType integration;

    std::pair<DynVector<double>, std::vector<Point<dimension>>> integration_rule =
        integration.GetRule(2 * raw_boundary.p + 1);

    uint ngp = integration_rule.first.size();

    std::vector<Point<dimension + 1>> z_master =
        this->master.BoundaryToMasterCoordinates(this->bound_id, integration_rule.second);

    // Global coordinates of gps
    this->gp_global_coordinates = this->shape.LocalToGlobalCoordinates(z_master);

    // Compute factors to expand nodal values
    this->psi_gp = this->shape.GetPsi(z_master);

    // Compute factors to expand boundary nodal values
    this->psi_bound_gp = this->shape.GetBoundaryPsi(this->bound_id, integration_rule.second);

    // Compute factors to expand modal values
    this->phi_gp = raw_boundary.basis.GetPhi(raw_boundary.p, z_master);

    DynVector<double> surface_J = this->shape.GetSurfaceJ(this->bound_id, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_fact = integration_rule.first * surface_J[0];

        this->int_phi_fact = transpose(this->phi_gp);
        for (uint dof = 0; dof < this->master.ndof; ++dof) {
            for (uint gp = 0; gp < ngp; ++gp) {
                this->int_phi_fact(gp, dof) *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_phi_phi_fact.resize(ngp, std::pow(this->master.ndof, 2));
        for (uint dof_i = 0; dof_i < this->master.ndof; ++dof_i) {
            for (uint dof_j = 0; dof_j < this->master.ndof; ++dof_j) {
                uint lookup = this->master.ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; ++gp) {
                    this->int_phi_phi_fact(gp, lookup) = this->phi_gp(dof_i, gp) * this->int_phi_fact(gp, dof_j);
                }
            }
        }

        StatVector<double, dimension + 1> normal = this->shape.GetSurfaceNormal(this->bound_id, z_master)[0];

        this->surface_normal.fill(DynRowVector<double>(ngp));

        for (uint gp = 0; gp < ngp; ++gp) {
            this->surface_normal[GlobalCoord::x][gp] = normal[GlobalCoord::x];
            this->surface_normal[GlobalCoord::y][gp] = normal[GlobalCoord::y];
        }
    }

    this->data.set_ngp_boundary(this->bound_id, ngp);
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename F>
inline DynMatrix<double> Boundary<dimension, IntegrationType, DataType, ConditonType>::ComputeFgp(const F& f) {
    uint nvar = f(this->gp_global_coordinates[0]).size();
    uint ngp  = this->gp_global_coordinates.size();

    DynMatrix<double> f_vals(nvar, ngp);

    for (uint gp = 0; gp < this->gp_global_coordinates.size(); ++gp) {
        column(f_vals, gp) = f(this->gp_global_coordinates[gp]);
    }

    return f_vals;
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename InputArrayType>
inline decltype(auto) Boundary<dimension, IntegrationType, DataType, ConditonType>::ComputeUgp(
    const InputArrayType& u) {
    // u_gp(q, gp) = u(q, dof) * phi_gp(dof, gp)
    return u * this->phi_gp;
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename InputArrayType>
inline decltype(auto) Boundary<dimension, IntegrationType, DataType, ConditonType>::ComputeNodalUgp(
    const InputArrayType& u_nodal) {
    // u_nodal_gp(q, gp) = u_nodal(q, dof) * psi_gp(dof, gp)
    return u_nodal * this->psi_gp;
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename InputArrayType>
inline decltype(auto) Boundary<dimension, IntegrationType, DataType, ConditonType>::ComputeBoundaryNodalUgp(
    const InputArrayType& u_bound_nodal) {
    // u_nodal_gp(q, gp) = u_nodal(q, dof) * psi_gp(dof, gp)
    return u_bound_nodal * this->psi_bound_gp;
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename InputArrayType>
inline decltype(auto) Boundary<dimension, IntegrationType, DataType, ConditonType>::Integration(
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_fact[gp]
    return u_gp * this->int_fact;
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename InputArrayType>
inline decltype(auto) Boundary<dimension, IntegrationType, DataType, ConditonType>::IntegrationPhi(
    const uint dof,
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_phi_fact(gp, dof)
    return u_gp * column(this->int_phi_fact, dof);
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename InputArrayType>
inline decltype(auto) Boundary<dimension, IntegrationType, DataType, ConditonType>::IntegrationPhi(
    const InputArrayType& u_gp) {
    // integral(q, dof) =  u_gp(q, gp) * int_phi_fact(gp, dof)
    return u_gp * this->int_phi_fact;
}

template <uint dimension, typename IntegrationType, typename DataType, typename ConditonType>
template <typename InputArrayType>
inline decltype(auto) Boundary<dimension, IntegrationType, DataType, ConditonType>::IntegrationPhiPhi(
    const uint dof_i,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_phi_phi_fact(gp, lookup)
    uint lookup = this->master.ndof * dof_i + dof_j;

    return u_gp * column(this->int_phi_phi_fact, lookup);
}
}

#endif
