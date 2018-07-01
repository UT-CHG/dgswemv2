#ifndef CLASS_EDGE_BOUNDARY_HPP
#define CLASS_EDGE_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
class EdgeBoundary {
  public:
    EdgeDataType edge_data;

    BoundaryType& boundary;

  private:
    Array2D<double> lambda_gp;
    Array2D<double> int_lambda_fact;
    Array3D<double> int_lambda_lambda_fact;

    std::pair<bool, Array2D<double>> m_inv;

  public:
    EdgeBoundary(BoundaryType& boundary);

    void L2Projection(const std::vector<double>& u_gp, std::vector<double>& projection);

    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);
    double IntegrationLambda(const uint dof, const std::vector<double>& u_gp);
    double IntegrationLambdaLambda(const uint dof_i, const uint dof_j, const std::vector<double>& u_gp);

    void ApplyMinv(const std::vector<double>& rhs, std::vector<double>& solution);
};

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::EdgeBoundary(BoundaryType& boundary)
    : boundary(boundary) {
    // *** //
    typename BoundaryType::BoundaryIntegrationType integration;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule =
        integration.GetRule(2 * this->boundary.GetMaster().p + 1);

    BasisType basis;

    this->lambda_gp = basis.GetPhi(this->boundary.GetMaster().p, integration_rule.second);

    std::vector<Point<dimension + 1>> z_master =
        this->boundary.GetMaster().BoundaryToMasterCoordinates(this->boundary.bound_id, integration_rule.second);

    std::vector<double> surface_J = this->boundary.GetShape().GetSurfaceJ(this->boundary.bound_id, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_lambda_fact = this->lambda_gp;
        for (uint dof = 0; dof < this->int_lambda_fact.size(); dof++) {
            for (uint gp = 0; gp < this->int_lambda_fact[dof].size(); gp++) {
                this->int_lambda_fact[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_lambda_lambda_fact.resize(this->lambda_gp.size());
        for (uint dof_i = 0; dof_i < this->lambda_gp.size(); dof_i++) {
            this->int_lambda_lambda_fact[dof_i] = this->int_lambda_fact;
            for (uint dof_j = 0; dof_j < this->lambda_gp.size(); dof_j++) {
                for (uint gp = 0; gp < this->int_lambda_lambda_fact[dof_i][dof_j].size(); gp++) {
                    this->int_lambda_lambda_fact[dof_i][dof_j][gp] *= this->lambda_gp[dof_i][gp];
                }
            }
        }

        this->m_inv = basis.GetMinv(this->boundary.GetMaster().p);
        for (uint i = 0; i < this->m_inv.second.size(); i++) {
            for (uint j = 0; j < this->m_inv.second[i].size(); j++) {
                this->m_inv.second[i][j] /= surface_J[0];
            }
        }
    }

    uint ndof = this->lambda_gp.size();
    uint ngp  = integration_rule.first.size();

    this->edge_data = EdgeDataType(ndof, ngp);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
void EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::L2Projection(const std::vector<double>& u_gp,
                                                                                  std::vector<double>& projection) {
    std::vector<double> rhs;

    for (uint dof = 0; dof < this->lambda_gp.size(); dof++) {
        rhs.push_back(this->IntegrationLambda(dof, u_gp));
    }

    this->ApplyMinv(rhs, projection);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
void EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::ComputeUgp(const std::vector<double>& u,
                                                                                std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->lambda_gp[dof][gp];
        }
    }
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
double EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::IntegrationLambda(
    const uint dof,
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_fact[dof][gp];
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
double EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::IntegrationLambdaLambda(
    const uint dof_i,
    const uint dof_j,
    const std::vector<double>& u_gp) {
    // *** //
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_lambda_fact[dof_i][dof_j][gp];
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename BoundaryType>
void EdgeBoundary<dimension, BasisType, EdgeDataType, BoundaryType>::ApplyMinv(const std::vector<double>& rhs,
                                                                               std::vector<double>& solution) {
    if (this->m_inv.first) {  // diagonal
        for (uint i = 0; i < rhs.size(); i++) {
            solution[i] = this->m_inv.second[0][i] * rhs[i];
        }
    } else if (!(this->m_inv.first)) {  // not diagonal
        for (uint i = 0; i < this->m_inv.second.size(); i++) {
            solution[i] = 0.0;
            for (uint j = 0; j < rhs.size(); j++) {
                solution[i] += this->m_inv.second[i][j] * rhs[j];
            }
        }
    }
}
}

#endif