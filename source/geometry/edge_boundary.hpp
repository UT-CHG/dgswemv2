#ifndef CLASS_EDGE_BOUNDARY_HPP
#define CLASS_EDGE_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
class EdgeBoundary {
  public:
    EdgeDataType edge_data;

    BoundaryType& boundary;

  private:
    uint ID;

    DynMatrix<double> lambda_gp;
    DynMatrix<double> int_lambda_fact;
    DynMatrix<double> int_lambda_lambda_fact;
    DynMatrix<double> int_phi_lambda_fact;

    DynMatrix<double> m_inv;

  public:
    EdgeBoundary(BoundaryType& boundary, const bool ccw = true);

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
    decltype(auto) IntegrationPhiLambda(const uint dof_i, const uint dof_j, const InputArrayType& u_gp);

    template <typename InputArrayType>
    decltype(auto) ApplyMinv(const InputArrayType& rhs);
};

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::EdgeBoundary(BoundaryType& boundary, const bool ccw)
    : boundary(boundary) {
    // *** //
    typename BoundaryType::BoundaryIntegrationType integration;

    std::pair<DynVector<double>, std::vector<Point<dimension>>> integration_rule =
        integration.GetRule(2 * this->boundary.GetMaster().p + 1);

    uint ngp = integration_rule.first.size();

    BasisType basis;

    if (ccw) {  // if trace basis is defined counterclockwise around element
        this->lambda_gp = basis.GetPhi(this->boundary.GetMaster().p, integration_rule.second);
    } else {
        std::vector<Point<dimension>> gauss_points = integration_rule.second;

        for (uint gp = 0; gp < (ngp / 2); ++gp) {
            std::swap(gauss_points[gp], gauss_points[ngp - gp - 1]);
        }

        this->lambda_gp = basis.GetPhi(this->boundary.GetMaster().p, gauss_points);
    }

    uint ndof = row(this->lambda_gp, 0).size();

    std::vector<Point<dimension + 1>> z_master =
        this->boundary.GetMaster().BoundaryToMasterCoordinates(this->boundary.bound_id, integration_rule.second);

    DynVector<double> surface_J = this->boundary.GetShape().GetSurfaceJ(this->boundary.bound_id, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_lambda_fact = transpose(this->lambda_gp);
        for (uint dof = 0; dof < ndof; ++dof) {
            for (uint gp = 0; gp < ngp; ++gp) {
                this->int_lambda_fact(gp, dof) *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_lambda_lambda_fact.resize(ngp, std::pow(ndof, 2));
        for (uint dof_i = 0; dof_i < ndof; ++dof_i) {
            for (uint dof_j = 0; dof_j < ndof; ++dof_j) {
                uint lookup = ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; ++gp) {
                    this->int_lambda_lambda_fact(gp, lookup) =
                        this->lambda_gp(dof_i, gp) * this->int_lambda_fact(gp, dof_j);
                }
            }
        }

        uint ndof_loc = this->boundary.data.get_ndof();
        this->int_phi_lambda_fact.resize(ngp, ndof * ndof_loc);
        for (uint dof_i = 0; dof_i < ndof_loc; ++dof_i) {
            for (uint dof_j = 0; dof_j < ndof; ++dof_j) {
                uint lookup = ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; ++gp) {
                    this->int_phi_lambda_fact(gp, lookup) =
                        this->boundary.phi_gp(dof_i, gp) * this->int_lambda_fact(gp, dof_j);
                }
            }
        }

        this->m_inv = basis.GetMinv(this->boundary.GetMaster().p) / surface_J[0];
    }

    this->edge_data.set_ndof(ndof);
    this->edge_data.set_ngp(ngp);

    this->edge_data.initialize();
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename InputArrayType>
inline decltype(auto) EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::L2Projection(
    const InputArrayType& u_gp) {
    // projection(q, dof) = gp_values(q, gp) * int_lambda_fact(gp, dof) * m_inv(dof, dof)
    InputArrayType projection = u_gp * this->int_lambda_fact * this->m_inv;

    return projection;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename InputArrayType>
inline decltype(auto) EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::ComputeUgp(
    const InputArrayType& u) {
    // u_gp(q, gp) = u(q, dof) * lambda_gp(dof, gp)
    return u * this->lambda_gp;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename InputArrayType>
inline decltype(auto) EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::IntegrationLambda(
    const uint dof,
    const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * int_lambda_fact(gp, dof)
    return u_gp * column(this->int_lambda_fact, dof);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename InputArrayType>
inline decltype(auto) EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::IntegrationLambdaLambda(
    const uint dof_i,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] = u_gp(q, gp) * int_lambda_lambda_fact(gp, lookup)
    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    return u_gp * column(this->int_lambda_lambda_fact, lookup);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename InputArrayType>
inline decltype(auto) EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::IntegrationPhiLambda(
    const uint dof_i,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] = u_gp(q, gp) * int_phi_lambda_fact(gp, lookup)
    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    return u_gp * column(this->int_phi_lambda_fact, lookup);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename InputArrayType>
inline decltype(auto) EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::ApplyMinv(
    const InputArrayType& rhs) {
    // solution(q, dof) = rhs(q, dof) * this->m_inv(dof, dof)
    return rhs * this->m_inv;
}
}

#endif