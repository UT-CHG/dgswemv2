#ifndef CLASS_EDGE_INTERFACE_HPP
#define CLASS_EDGE_INTERFACE_HPP

namespace Geometry {
template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
class EdgeInterface {
  public:
    EdgeDataType edge_data;

    InterfaceType& interface;

  private:
    Array2D<double> lambda_gp;
    Array2D<double> int_lambda_fact;
    Array3D<double> int_lambda_lambda_fact;

    std::pair<bool, Array2D<double>> m_inv;

  public:
    EdgeInterface(InterfaceType& interface);

    void L2Projection(const std::vector<double>& u_gp, std::vector<double>& projection);

    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);
    double IntegrationLambda(const uint dof, const std::vector<double>& u_gp);
    double IntegrationLambdaLambda(const uint dof_i, const uint dof_j, const std::vector<double>& u_gp);

    void ApplyMinv(const std::vector<double>& rhs, std::vector<double>& solution);
};

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::EdgeInterface(InterfaceType& interface)
    : interface(interface) {
    // *** //
    uint p = std::max(this->interface.GetMasterIN().p, this->interface.GetMasterEX().p);

    typename InterfaceType::InterfaceIntegrationType integration;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * p + 1);

    BasisType basis;

    this->lambda_gp = basis.GetPhi(p, integration_rule.second);

    // transfrom gp to master coord in
    std::vector<Point<dimension + 1>> z_master_in =
        this->interface.GetMasterIN().BoundaryToMasterCoordinates(this->interface.bound_id_in, integration_rule.second);

    std::vector<double> surface_J = this->interface.GetShapeIN().GetSurfaceJ(this->interface.bound_id_in, z_master_in);

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

        this->m_inv = basis.GetMinv(p);
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

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
void EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::L2Projection(const std::vector<double>& u_gp,
                                                                                    std::vector<double>& projection) {
    std::vector<double> rhs;

    for (uint dof = 0; dof < this->lambda_gp.size(); dof++) {
        rhs.push_back(this->IntegrationLambda(dof, u_gp));
    }

    this->ApplyMinv(rhs, projection);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
void EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::ComputeUgp(const std::vector<double>& u,
                                                                                  std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->lambda_gp[dof][gp];
        }
    }
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
double EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationLambda(
    const uint dof,
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_fact[dof][gp];
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
double EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationLambdaLambda(
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

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
void EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::ApplyMinv(const std::vector<double>& rhs,
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