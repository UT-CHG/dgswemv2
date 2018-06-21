#ifndef CLASS_EDGE_INTERNAL_HPP
#define CLASS_EDGE_INTERNAL_HPP

namespace Geometry {
template <uint dimension, typename BasisType, typename IntegrationType, typename DataType, typename EdgeDataType>
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
    Array2D<double> lambda_gp;
    Array2D<double> int_lambda_fact;

  public:
    EdgeInternal(const RawBoundary<dimension, DataType>& raw_boundary_in,
                 const RawBoundary<dimension, DataType>& raw_boundary_ex);

    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);
    double IntegrationLambda(const uint dof, const std::vector<double>& u_gp);
};

template <uint dimension, typename BasisType, typename IntegrationType, typename DataType, typename EdgeDataType>
EdgeInternal<dimension, BasisType, IntegrationType, DataType, EdgeDataType>::EdgeInternal(
    const RawBoundary<dimension, DataType>& raw_boundary_in,
    const RawBoundary<dimension, DataType>& raw_boundary_ex)
    : data_in(raw_boundary_in.data),
      data_ex(raw_boundary_ex.data),
      bound_id_in(raw_boundary_in.bound_id),
      bound_id_ex(raw_boundary_ex.bound_id) {
    // *** //
    uint p = std::max(raw_boundary_in.p, raw_boundary_ex.p);

    IntegrationType integration;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * p + 1);

    BasisType basis;

    this->lambda_gp = basis.GetPhi(p, integration_rule.second);

    // transfrom gp to master coord in
    std::vector<Point<dimension + 1>> z_master_in =
        raw_boundary_in.master.BoundaryToMasterCoordinates(this->bound_id_in, integration_rule.second);

    std::vector<double> surface_J = raw_boundary_in.shape.GetSurfaceJ(this->bound_id_in, z_master_in);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_lambda_fact = this->lambda_gp;
        for (uint dof = 0; dof < this->int_lambda_fact.size(); dof++) {
            for (uint gp = 0; gp < this->int_lambda_fact[dof].size(); gp++) {
                this->int_lambda_fact[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->surface_normal_in =
            Array2D<double>(integration_rule.first.size(),
                            *raw_boundary_in.shape.GetSurfaceNormal(this->bound_id_in, z_master_in).begin());

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
}

template <uint dimension, typename BasisType, typename IntegrationType, typename DataType, typename EdgeDataType>
void EdgeInternal<dimension, BasisType, IntegrationType, DataType, EdgeDataType>::ComputeUgp(
    const std::vector<double>& u,
    std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->lambda_gp[dof][gp];
        }
    }
}

template <uint dimension, typename BasisType, typename IntegrationType, typename DataType, typename EdgeDataType>
double EdgeInternal<dimension, BasisType, IntegrationType, DataType, EdgeDataType>::IntegrationLambda(
    const uint dof,
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_fact[dof][gp];
    }

    return integral;
}
}

#endif