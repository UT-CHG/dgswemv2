#ifndef CLASS_EDGE_INTERNAL_HPP
#define CLASS_EDGE_INTERNAL_HPP

namespace Geometry {
template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
class EdgeInternal {
  public:
    EdgeDataType edge_data;

    DataType& data_in;
    DataType& data_ex;

    uint bound_id_in;
    uint bound_id_ex;

    Array2D<double> surface_normal_in;
    Array2D<double> surface_normal_ex;

  private:
    Master::Master<dimension + 1>& master_in;
    Master::Master<dimension + 1>& master_ex;

    Shape::Shape<dimension + 1>& shape_in;
    Shape::Shape<dimension + 1>& shape_ex;

    Array2D<double> lambda_gp;
    Array2D<double> int_lambda_fact;
    Array3D<double> int_lambda_lambda_fact;

    std::pair<bool, Array2D<double>> m_inv;

  public:
    template <typename InterfaceType>
    EdgeInternal(InterfaceType& intface);

    Master::Master<dimension + 1>& GetMasterIN() { return this->master_in; }
    Master::Master<dimension + 1>& GetMasterEX() { return this->master_ex; }

    Shape::Shape<dimension + 1>& GetShapeIN() { return this->shape_in; }
    Shape::Shape<dimension + 1>& GetShapeEX() { return this->shape_ex; }

    void L2Projection(const std::vector<double>& u_gp, std::vector<double>& projection);

    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);
    double IntegrationLambda(const uint dof, const std::vector<double>& u_gp);
    double IntegrationLambdaLambda(const uint dof_i, const uint dof_j, const std::vector<double>& u_gp);

    void ApplyMinv(const std::vector<double>& rhs, std::vector<double>& solution);
};

template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
template <typename InterfaceType>
EdgeInternal<dimension, BasisType, DataType, EdgeDataType>::EdgeInternal(InterfaceType& intface)
    : data_in(intface.data_in),
      data_ex(intface.data_ex),
      bound_id_in(intface.bound_id_in),
      bound_id_ex(intface.bound_id_ex),
      master_in(intface.GetMasterIN()),
      master_ex(intface.GetMasterEX()),
      shape_in(intface.GetShapeIN()),
      shape_ex(intface.GetShapeEX()) {
    // *** //
    uint p = std::max(this->master_in.p, this->master_ex.p);

    typename InterfaceType::InterfaceIntegrationType integration;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * p + 1);

    BasisType basis;

    this->lambda_gp = basis.GetPhi(p, integration_rule.second);

    // transfrom gp to master coord in
    std::vector<Point<dimension + 1>> z_master_in =
        this->master_in.BoundaryToMasterCoordinates(this->bound_id_in, integration_rule.second);

    std::vector<double> surface_J = this->shape_in.GetSurfaceJ(this->bound_id_in, z_master_in);

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

        this->surface_normal_in = Array2D<double>(
            integration_rule.first.size(), *this->shape_in.GetSurfaceNormal(this->bound_id_in, z_master_in).begin());

        this->surface_normal_ex = this->surface_normal_in;  // same dimensions

        uint ngp   = this->surface_normal_ex.size();
        uint gp_ex = 0;
        for (uint gp = 0; gp < this->surface_normal_ex.size(); gp++) {
            gp_ex = ngp - gp - 1;
            for (uint dir = 0; dir < this->surface_normal_ex[gp].size(); dir++) {
                this->surface_normal_ex[gp_ex][dir] = -this->surface_normal_in[gp][dir];  // but opposite direction
            }
        }
    }

    uint ndof = this->lambda_gp.size();
    uint ngp  = integration_rule.first.size();

    this->edge_data = EdgeDataType(ndof, ngp);
}

template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
void EdgeInternal<dimension, BasisType, DataType, EdgeDataType>::L2Projection(const std::vector<double>& u_gp,
                                                                              std::vector<double>& projection) {
    std::vector<double> rhs;

    for (uint dof = 0; dof < this->lambda_gp.size(); dof++) {
        rhs.push_back(this->IntegrationLambda(dof, u_gp));
    }

    this->ApplyMinv(rhs, projection);
}

template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
void EdgeInternal<dimension, BasisType, DataType, EdgeDataType>::ComputeUgp(const std::vector<double>& u,
                                                                            std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->lambda_gp[dof][gp];
        }
    }
}

template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
double EdgeInternal<dimension, BasisType, DataType, EdgeDataType>::IntegrationLambda(const uint dof,
                                                                                     const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_fact[dof][gp];
    }

    return integral;
}

template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
double EdgeInternal<dimension, BasisType, DataType, EdgeDataType>::IntegrationLambdaLambda(
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

template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
void EdgeInternal<dimension, BasisType, DataType, EdgeDataType>::ApplyMinv(const std::vector<double>& rhs,
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