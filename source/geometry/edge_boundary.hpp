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

    std::pair<bool, DynMatrix<double>> m_inv;

  public:
    EdgeBoundary(BoundaryType& boundary);

    uint GetID() { return this->ID; }
    void SetID(uint ID) { this->ID = ID; }

    template <typename T>
    void L2Projection(const std::vector<T>& u_gp, std::vector<T>& projection);

    template <typename T>
    void ComputeUgp(const std::vector<T>& u, std::vector<T>& u_gp);

    template <typename T>
    T IntegrationLambda(const uint dof, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationLambdaLambda(const uint dof_i, const uint dof_j, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhiLambda(const uint dof_i, const uint dof_j, const std::vector<T>& u_gp);

    template <typename T>
    void ApplyMinv(const std::vector<T>& rhs, std::vector<T>& solution);
};

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::EdgeBoundary(BoundaryType& boundary)
    : boundary(boundary) {
    // *** //
    typename BoundaryType::BoundaryIntegrationType integration;

    std::pair<DynVector<double>, DynVector<Point<dimension>>> integration_rule =
        integration.GetRule(2 * this->boundary.GetMaster().p + 1);

    uint ngp = integration_rule.first.size();

    BasisType basis;

    this->lambda_gp = basis.GetPhi(this->boundary.GetMaster().p, integration_rule.second);

    uint ndof = row(this->lambda_gp, 0).size();

    DynVector<Point<dimension + 1>> z_master =
        this->boundary.GetMaster().BoundaryToMasterCoordinates(this->boundary.bound_id, integration_rule.second);

    DynVector<double> surface_J = this->boundary.GetShape().GetSurfaceJ(this->boundary.bound_id, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_lambda_fact = transpose(this->lambda_gp);
        for (uint dof = 0; dof < ndof; dof++) {
            for (uint gp = 0; gp < ngp; gp++) {
                this->int_lambda_fact(dof, gp) *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_lambda_lambda_fact.resize(std::pow(ndof, 2), ngp);
        for (uint dof_i = 0; dof_i < ndof; dof_i++) {
            for (uint dof_j = 0; dof_j < ndof; dof_j++) {
                uint lookup = ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; gp++) {
                    this->int_lambda_lambda_fact(lookup, gp) =
                        this->lambda_gp(gp, dof_i) * this->int_lambda_fact(dof_j, gp);
                }
            }
        }

        uint ndof_loc = this->boundary.data.get_ndof();
        this->int_phi_lambda_fact.resize(ndof * ndof_loc, ngp);
        for (uint dof_i = 0; dof_i < ndof_loc; dof_i++) {
            for (uint dof_j = 0; dof_j < ndof; dof_j++) {
                uint lookup = ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; gp++) {
                    this->int_phi_lambda_fact(lookup, gp) =
                        this->boundary.phi_gp(gp, dof_i) * this->int_lambda_fact(dof_j, gp);
                }
            }
        }

        this->m_inv = basis.GetMinv(this->boundary.GetMaster().p);
        this->m_inv.second /= surface_J[0];
    }

    this->edge_data.set_ndof(ndof);
    this->edge_data.set_ngp(ngp);

    this->edge_data.initialize();
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename T>
void EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::L2Projection(const std::vector<T>& u_gp,
                                                                                  std::vector<T>& projection) {
    std::vector<T> rhs;

    for (uint dof = 0; dof < this->edge_data.get_ndof(); dof++) {
        rhs.push_back(this->IntegrationLambda(dof, u_gp));
    }

    this->ApplyMinv(rhs, projection);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename T>
void EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::ComputeUgp(const std::vector<T>& u,
                                                                                std::vector<T>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->lambda_gp(gp, dof);
        }
    }
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename T>
T EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::IntegrationLambda(const uint dof,
                                                                                    const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_fact(dof, gp);
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename T>
T EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::IntegrationLambdaLambda(const uint dof_i,
                                                                                          const uint dof_j,
                                                                                          const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_lambda_fact(lookup, gp);
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename T>
T EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::IntegrationPhiLambda(const uint dof_i,
                                                                                       const uint dof_j,
                                                                                       const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_lambda_fact(lookup, gp);
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
template <typename T>
void EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::ApplyMinv(const std::vector<T>& rhs,
                                                                               std::vector<T>& solution) {
    if (this->m_inv.first) {  // diagonal
        for (uint i = 0; i < rhs.size(); i++) {
            solution[i] = this->m_inv.second(i, i) * rhs[i];
        }
    } else if (!(this->m_inv.first)) {  // not diagonal
        for (uint i = 0; i < rhs.size(); i++) {
            solution[i] = 0.0;
            for (uint j = 0; j < rhs.size(); j++) {
                solution[i] += this->m_inv.second(i, j) * rhs[j];
            }
        }
    }
}
}

#endif