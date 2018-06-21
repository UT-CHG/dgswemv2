#ifndef CLASS_EDGE_BOUNDARY_HPP
#define CLASS_EDGE_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, typename BasisType, typename IntegrationType, typename DataType, typename EdgeDataType>
class EdgeBoundary {
  public:
    EdgeDataType edge_data;

    DataType& data;

    uint bound_id;

    Array2D<double> surface_normal;

  private:
    Array2D<double> lambda_gp;
    Array2D<double> int_lambda_fact;

  public:
    EdgeBoundary(const RawBoundary<dimension, DataType>& raw_boundary);

    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);
    double IntegrationLambda(const uint dof, const std::vector<double>& u_gp);
};

template <uint dimension, typename BasisType, typename IntegrationType, typename DataType, typename EdgeDataType>
EdgeBoundary<dimension, BasisType, IntegrationType, DataType, EdgeDataType>::EdgeBoundary(
    const RawBoundary<dimension, DataType>& raw_boundary)
    : data(raw_boundary.data), bound_id(raw_boundary.bound_id) {
    // *** //
    IntegrationType integration;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule =
        integration.GetRule(2 * raw_boundary.p + 1);

    BasisType basis;

    this->lambda_gp = basis.GetPhi(raw_boundary.p, integration_rule.second);

    std::vector<Point<dimension + 1>> z_master =
        raw_boundary.master.BoundaryToMasterCoordinates(this->bound_id, integration_rule.second);

    std::vector<double> surface_J = raw_boundary.shape.GetSurfaceJ(this->bound_id, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_lambda_fact = this->lambda_gp;
        for (uint dof = 0; dof < this->int_lambda_fact.size(); dof++) {
            for (uint gp = 0; gp < this->int_lambda_fact[dof].size(); gp++) {
                this->int_lambda_fact[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->surface_normal = Array2D<double>(integration_rule.first.size(),
                                               *raw_boundary.shape.GetSurfaceNormal(this->bound_id, z_master).begin());
    }
}

template <uint dimension, typename BasisType, typename IntegrationType, typename DataType, typename EdgeDataType>
void EdgeBoundary<dimension, BasisType, IntegrationType, DataType, EdgeDataType>::ComputeUgp(
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
double EdgeBoundary<dimension, BasisType, IntegrationType, DataType, EdgeDataType>::IntegrationLambda(
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