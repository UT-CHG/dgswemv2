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

    std::pair<bool, DynMatrix<double>> m_inv;

  public:
    EdgeInterface(InterfaceType& interface);

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
    T IntegrationPhiLambdaIN(const uint dof_i, const uint dof_j, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhiLambdaEX(const uint dof_i, const uint dof_j, const std::vector<T>& u_gp);

    template <typename T>
    void ApplyMinv(const std::vector<T>& rhs, std::vector<T>& solution);
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

        uint ndof_loc_in = this->interface.data_in.get_ndof();
        this->int_phi_lambda_fact_in.resize(ndof * ndof_loc_in, ngp);
        for (uint dof_i = 0; dof_i < ndof_loc_in; dof_i++) {
            for (uint dof_j = 0; dof_j < ndof; dof_j++) {
                uint lookup = ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; gp++) {
                    this->int_phi_lambda_fact_in(lookup, gp) =
                        this->interface.phi_gp_in(gp, dof_i) * this->int_lambda_fact(dof_j, gp);
                }
            }
        }

        uint ndof_loc_ex = this->interface.data_ex.get_ndof();
        this->int_phi_lambda_fact_ex.resize(ndof * ndof_loc_ex, ngp);
        for (uint dof_i = 0; dof_i < ndof_loc_ex; dof_i++) {
            for (uint dof_j = 0; dof_j < ndof; dof_j++) {
                uint lookup = ndof * dof_i + dof_j;
                for (uint gp = 0; gp < ngp; gp++) {
                    uint gp_ex = ngp - gp - 1;
                    this->int_phi_lambda_fact_ex(lookup, gp) =
                        this->interface.phi_gp_ex(gp, dof_i) * this->int_lambda_fact(dof_j, gp_ex);
                }
            }
        }

        this->m_inv = basis.GetMinv(p);
        this->m_inv.second /= surface_J[0];
    }

    this->edge_data.set_ndof(ndof);
    this->edge_data.set_ngp(ngp);

    this->edge_data.initialize();
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename T>
void EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::L2Projection(const std::vector<T>& u_gp,
                                                                                    std::vector<T>& projection) {
    std::vector<T> rhs;

    for (uint dof = 0; dof < this->edge_data.get_ndof(); dof++) {
        rhs.push_back(this->IntegrationLambda(dof, u_gp));
    }

    this->ApplyMinv(rhs, projection);
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename T>
void EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::ComputeUgp(const std::vector<T>& u,
                                                                                  std::vector<T>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->lambda_gp(gp, dof);
        }
    }
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename T>
T EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationLambda(const uint dof,
                                                                                      const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_fact(dof, gp);
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename T>
T EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationLambdaLambda(
    const uint dof_i,
    const uint dof_j,
    const std::vector<T>& u_gp) {
    // *** //
    T integral;

    integral = 0.0;

    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_lambda_fact(lookup, gp);
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename T>
T EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationPhiLambdaIN(const uint dof_i,
                                                                                           const uint dof_j,
                                                                                           const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_lambda_fact_in(lookup, gp);
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename T>
T EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::IntegrationPhiLambdaEX(const uint dof_i,
                                                                                           const uint dof_j,
                                                                                           const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    uint lookup = this->edge_data.get_ndof() * dof_i + dof_j;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_lambda_fact_ex(lookup, gp);
    }

    return integral;
}

template <uint dimension, typename BasisType, typename EdgeDataType, typename InterfaceType>
template <typename T>
void EdgeInterface<dimension, BasisType, EdgeDataType, InterfaceType>::ApplyMinv(const std::vector<T>& rhs,
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