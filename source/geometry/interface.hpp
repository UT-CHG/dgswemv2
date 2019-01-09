#ifndef CLASS_INTERFACE_HPP
#define CLASS_INTERFACE_HPP

#include "general_definitions.hpp"

namespace Geometry {
template <uint dim, typename IntegrationType, typename DataType, typename SpecializationType>
class Interface {
  public:
    constexpr static uint dimension = dim;
    using specialization_t = SpecializationType;
    using InterfaceIntegrationType = IntegrationType;

  public:
    SpecializationType specialization;

    DataType& data_in;
    DataType& data_ex;

    uint bound_id_in;
    uint bound_id_ex;

    std::array<DynRowVector<double>, dimension +1> surface_normal_in;
    std::array<DynRowVector<double>, dimension +1> surface_normal_ex;

  private:
    Master::Master<dimension + 1>& master_in;
    Master::Master<dimension + 1>& master_ex;

    Shape::Shape<dimension + 1>& shape_in;
    Shape::Shape<dimension + 1>& shape_ex;

    std::vector<uint> node_ID_in;
    std::vector<uint> node_ID_ex;

    DynMatrix<double> psi_gp_in;
    DynMatrix<double> psi_gp_ex;
    DynMatrix<double> psi_bound_gp_in;
    DynMatrix<double> psi_bound_gp_ex;
    DynMatrix<double> phi_gp_in;
    DynMatrix<double> phi_gp_ex;

    DynVector<double> int_fact_in;
    DynVector<double> int_fact_ex;
    DynMatrix<double> int_phi_fact_in;
    DynMatrix<double> int_phi_fact_ex;
    DynMatrix<double> int_phi_phi_fact_in;
    DynMatrix<double> int_phi_phi_fact_ex;

  public:
    template <typename... Args>
    Interface(RawBoundary<dimension, DataType>&& raw_boundary_in,
              RawBoundary<dimension, DataType>&& raw_boundary_ex,
              Args&&... args);

    Master::Master<dimension + 1>& GetMasterIN() { return this->master_in; }
    Master::Master<dimension + 1>& GetMasterEX() { return this->master_ex; }

    Shape::Shape<dimension + 1>& GetShapeIN() { return this->shape_in; }
    Shape::Shape<dimension + 1>& GetShapeEX() { return this->shape_ex; }

    std::vector<uint>& GetNodeIDIN() { return this->node_ID_in; }
    std::vector<uint>& GetNodeIDEX() { return this->node_ID_ex; }

    template <typename InputArrayType>
    decltype(auto) ComputeUgpIN(const InputArrayType& u);
    template <typename InputArrayType>
    decltype(auto) ComputeUgpEX(const InputArrayType& u);

    template <typename InputArrayType>
    decltype(auto) ComputeNodalUgpIN(const InputArrayType& u_nodal);
    template <typename InputArrayType>
    decltype(auto) ComputeNodalUgpEX(const InputArrayType& u_nodal);
    template <typename InputArrayType>
    decltype(auto) ComputeBoundaryNodalUgpIN(const InputArrayType& u_bound_nodal);
    template <typename InputArrayType>
    decltype(auto) ComputeBoundaryNodalUgpEX(const InputArrayType& u_bound_nodal);

    template <typename InputArrayType>
    decltype(auto) IntegrationIN(const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationEX(const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiIN(const uint dof, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiEX(const uint dof, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiIN(const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiEX(const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiPhiIN(const uint dof_i, const uint dof_j, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiPhiEX(const uint dof_i, const uint dof_j, const InputArrayType& u_gp);

  private:
    template <uint d, typename BT, typename EDT, typename InT>
    friend class EdgeInterface;
};

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename... Args>
Interface<dimension, IntegrationType, DataType, SpecializationType>::Interface(
    RawBoundary<dimension, DataType>&& raw_boundary_in,
    RawBoundary<dimension, DataType>&& raw_boundary_ex,
    Args&&... args)
    : specialization(std::forward<Args>(args)...),
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

    std::pair<DynVector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * p + 1);

    uint ngp = integration_rule.first.size();

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
    this->psi_gp_ex = reverse_columns(this->shape_ex.GetPsi(z_master_ex));

    // Compute factors to expand boundary nodal values ex
    this->psi_bound_gp_ex = reverse_columns(this->shape_ex.GetBoundaryPsi(this->bound_id_ex, integration_rule.second));

    // Compute factors to expand modal values ex
    this->phi_gp_ex = reverse_columns(raw_boundary_ex.basis.GetPhi(raw_boundary_ex.p, z_master_ex));

    DynVector<double> surface_J = this->shape_in.GetSurfaceJ(this->bound_id_in, z_master_in);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_fact_in = integration_rule.first * surface_J[0];

        this->int_fact_ex = reverse(integration_rule.first * surface_J[0]);

        this->int_phi_fact_in = transpose(this->phi_gp_in);
        for (uint dof = 0; dof < this->master_in.ndof; ++dof) {
            for (uint gp = 0; gp < ngp; ++gp) {
                this->int_phi_fact_in(gp, dof) *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_phi_fact_ex = transpose(this->phi_gp_ex);
        for (uint dof = 0; dof < this->master_ex.ndof; ++dof) {
            for (uint gp = 0; gp < ngp; ++gp) {
                this->int_phi_fact_ex(ngp-1 - gp, dof) *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_phi_phi_fact_in.resize(ngp, std::pow(this->master_in.ndof, 2));
        for (uint dof_i = 0; dof_i < this->master_in.ndof; ++dof_i) {
            for (uint dof_j = 0; dof_j < this->master_in.ndof; ++dof_j) {
                uint lookup = this->master_in.ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; ++gp) {
                    this->int_phi_phi_fact_in(gp, lookup) =
                        this->phi_gp_in(dof_i, gp) * this->int_phi_fact_in(gp, dof_j);
                }
            }
        }

        this->int_phi_phi_fact_ex.resize(ngp, std::pow(this->master_ex.ndof, 2));
        for (uint dof_i = 0; dof_i < this->master_ex.ndof; ++dof_i) {
            for (uint dof_j = 0; dof_j < this->master_ex.ndof; ++dof_j) {
                uint lookup = this->master_ex.ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; ++gp) {
                    this->int_phi_phi_fact_ex(gp, lookup) =
                        this->phi_gp_ex(dof_i, gp) * this->int_phi_fact_ex(gp, dof_j);
                }
            }
        }

        StatVector<double, dimension + 1> normal = this->shape_in.GetSurfaceNormal(this->bound_id_in, z_master_in)[0];

        this->surface_normal_in.fill(DynRowVector<double>(ngp));
        this->surface_normal_ex.fill(DynRowVector<double>(ngp));

        for (uint gp = 0; gp < ngp; ++gp) {
            this->surface_normal_in[GlobalCoord::x][gp]    = normal[GlobalCoord::x];
            this->surface_normal_in[GlobalCoord::y][gp]    = normal[GlobalCoord::y];

            this->surface_normal_ex[GlobalCoord::x][gp] = -normal[GlobalCoord::x];
            this->surface_normal_ex[GlobalCoord::y][gp] = -normal[GlobalCoord::y];
        }
    }

    this->data_in.set_ngp_boundary(this->bound_id_in, ngp);
    this->data_ex.set_ngp_boundary(this->bound_id_ex, ngp);
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeUgpIN(
    const InputArrayType& u) {
    // u_gp(q, gp) = u(q, dof) * phi_gp(dof, gp)
    return u * this->phi_gp_in;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeNodalUgpIN(
    const InputArrayType& u_nodal) {
    // u_nodal_gp(q, gp) = u_nodal(q, dof) * psi_gp(dof, gp)
    return u_nodal * this->psi_gp_in;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeBoundaryNodalUgpIN(
    const InputArrayType& u_bound_nodal) {
    // u_nodal_gp(q, gp) = u_nodal(q, dof) * psi_gp(dof, gp)
    return u_bound_nodal * this->psi_bound_gp_in;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationIN(
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_fact[gp]
    return u_gp * this->int_fact_in;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiIN(
    const uint dof,
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_phi_fact(gp, dof)
    return u_gp * column(this->int_phi_fact_in, dof);
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiIN(
    const InputArrayType& u_gp) {
    // integral(q, dof) =  u_gp(q, gp) * int_phi_fact(gp, dof)
    return u_gp * this->int_phi_fact_in;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiPhiIN(
    const uint dof_i,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_phi_phi_fact(gp, lookup)
    uint lookup = this->master_in.ndof * dof_i + dof_j;

    return u_gp * column(this->int_phi_phi_fact_in, lookup);
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeUgpEX(
    const InputArrayType& u) {
    // u_gp(q, gp) = u(q, dof) * phi_gp(dof, gp)
    return u * this->phi_gp_ex;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeNodalUgpEX(
    const InputArrayType& u_nodal) {
    // u_nodal_gp(q, gp) = u_nodal(q, dof) * psi_gp(dof, gp)
    return u_nodal * this->psi_gp_ex;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::ComputeBoundaryNodalUgpEX(
    const InputArrayType& u_bound_nodal) {
    // u_nodal_gp(q, gp) = u_nodal(q, dof) * psi_gp(dof, gp)
    return u_bound_nodal * this->psi_bound_gp_ex;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationEX(
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_fact[gp]
    return u_gp * this->int_fact_ex;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiEX(
    const uint dof,
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_phi_fact(gp, dof)
    return u_gp * column(this->int_phi_fact_ex, dof);
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiEX(
    const InputArrayType& u_gp) {
    // integral(q, dof) =  u_gp(q, gp) * int_phi_fact(gp, dof)
    return u_gp * this->int_phi_fact_ex;
}

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType>
template <typename InputArrayType>
inline decltype(auto) Interface<dimension, IntegrationType, DataType, SpecializationType>::IntegrationPhiPhiEX(
    const uint dof_i,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_phi_phi_fact(gp, lookup)
    uint lookup = this->master_ex.ndof * dof_i + dof_j;

    return u_gp * column(this->int_phi_phi_fact_ex, lookup);
}
}

#endif
