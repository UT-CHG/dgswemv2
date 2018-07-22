#ifndef CLASS_EDGE_INTERFACE_HPP
#define CLASS_EDGE_INTERFACE_HPP

namespace Geometry {
template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
class EdgeInterface {
  public:
    EdgeDataType edge_data;

    InterfaceType& interface;

  private:
    uint ID;

    DynMatrix<double> lambda_gp;
    DynMatrix<double> int_lambda_fact;
    DynMatrix<double> int_lambda_lambda_fact;
    DynMatrix<double> int_phi_lambda_fact_in;
    DynMatrix<double> int_phi_lambda_fact_ex;

    DynMatrix<double> m_inv;

  public:
    EdgeInterface(InterfaceType& interface);

    uint GetID() { return this->ID; }
    void SetID(uint ID) { this->ID = ID; }

    template <typename InputArrayType>
    decltype(auto) L2Projection(const InputArrayType& u_gp);

    template <typename InputArrayType>
    decltype(auto) ComputeUgp(const InputArrayType& u);

    template <typename InputArrayType>
    decltype(auto) IntegrationLambda(const uint dof, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationLambdaLambda(const uint dof_i, const uint dof_j, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiLambdaIN(const uint dof_i, const uint dof_j, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiLambdaEX(const uint dof_i, const uint dof_j, const InputArrayType& u_gp);

    template <typename InputArrayType>
    decltype(auto) ApplyMinv(const InputArrayType& rhs);
};

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::EdgeInterface(InterfaceType& interface)
    : interface(interface) {
    // *** //
    uint p = std::max(this->interface.GetMasterIN().p, this->interface.GetMasterEX().p);

    typename InterfaceType::InterfaceIntegrationType integration;

    std::pair<DynVector<double>, DynVector<Point<dimension>>> integration_rule = integration.GetRule(2 * p + 1);

    uint ngp = integration_rule.first.size();

    BasisType basis;

    this->lambda_gp = basis.GetPhi(p, integration_rule.second);

    uint ndof = row(this->lambda_gp, 0).size();

    // transfrom gp to master coord in
    DynVector<Point<dimension + 1>> z_master_in =
        this->interface.GetMasterIN().BoundaryToMasterCoordinates(this->interface.bound_id_in, integration_rule.second);

    DynVector<double> surface_J = this->interface.GetShapeIN().GetSurfaceJ(this->interface.bound_id_in, z_master_in);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_lambda_fact = transpose(this->lambda_gp);
        for (uint dof = 0; dof < ndof; dof++) {
            for (uint gp = 0; gp < ngp; gp++) {
                this->int_lambda_fact(gp, dof) *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_lambda_lambda_fact.resize(ngp, std::pow(ndof, 2));
        for (uint dof_i = 0; dof_i < ndof; dof_i++) {
            for (uint dof_j = 0; dof_j < ndof; dof_j++) {
                uint lookup = ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; gp++) {
                    this->int_lambda_lambda_fact(gp, lookup) =
                        this->lambda_gp(dof_i, gp) * this->int_lambda_fact(gp, dof_j);
                }
            }
        }

        uint ndof_loc_in = this->interface.data_in.get_ndof();
        this->int_phi_lambda_fact_in.resize(ngp, ndof * ndof_loc_in);
        for (uint dof_i = 0; dof_i < ndof_loc_in; dof_i++) {
            for (uint dof_j = 0; dof_j < ndof; dof_j++) {
                uint lookup = ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; gp++) {
                    this->int_phi_lambda_fact_in(gp, lookup) =
                        this->interface.phi_gp_in(dof_i, gp) * this->int_lambda_fact(gp, dof_j);
                }
            }
        }

        uint ndof_loc_ex = this->interface.data_ex.get_ndof();
        this->int_phi_lambda_fact_ex.resize(ngp, ndof * ndof_loc_ex);
        for (uint dof_i = 0; dof_i < ndof_loc_ex; dof_i++) {
            for (uint dof_j = 0; dof_j < ndof; dof_j++) {
                uint lookup = ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; gp++) {
                    uint gp_ex = ngp - gp - 1;
                    this->int_phi_lambda_fact_ex(gp, lookup) =
                        this->interface.phi_gp_ex(dof_i, gp) * this->int_lambda_fact(gp_ex, dof_j);
                }
            }
        }

        this->m_inv = basis.GetMinv(p) / surface_J[0];
    }

    this->edge_data.set_ndof(ndof);
    this->edge_data.set_ngp(ngp);

    this->edge_data.initialize();
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename InputArrayType>
inline decltype(auto) EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::L2Projection(
    const InputArrayType& u_gp) {
    // *** //
    return this->ApplyMinv(u_gp * this->int_lambda_fact);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename InputArrayType>
inline decltype(auto) EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::ComputeUgp(
    const InputArrayType& u) {
    // u_gp(q, gp) = u(q, dof) * lambda_gp(dof, gp)
    return u * this->lambda_gp;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename InputArrayType>
inline decltype(auto) EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationLambda(
    const uint dof,
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_lambda_fact(gp, dof)
    return u_gp * column(this->int_lambda_fact, dof);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename InputArrayType>
inline decltype(auto) EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationLambdaLambda(
    const uint dof_i,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] = u_gp(q, gp) * int_lambda_lambda_fact(gp, lookup)
    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    return u_gp * column(this->int_lambda_lambda_fact, lookup);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename InputArrayType>
inline decltype(auto) EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationPhiLambdaIN(
    const uint dof_i,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] = u_gp(q, gp) * int_phi_lambda_fact(gp, lookup)
    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    return u_gp * column(this->int_phi_lambda_fact_in, lookup);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename InputArrayType>
inline decltype(auto) EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationPhiLambdaEX(
    const uint dof_i,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] = u_gp(q, gp) * int_phi_lambda_fact(gp, lookup)
    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    return u_gp * column(this->int_phi_lambda_fact_ex, lookup);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename InputArrayType>
inline decltype(auto) EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::ApplyMinv(
    const InputArrayType& rhs) {
    // solution(q, dof) = rhs(q, dof) * this->m_inv(dof, dof)
    return rhs * this->m_inv;
}
}

#endif