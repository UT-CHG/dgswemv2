#ifndef CLASS_EDGE_BOUNDARY_HPP
#define CLASS_EDGE_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
class EdgeBoundary {
  public:
    EdgeDataType edge_data;

    DataType& data;

    uint bound_id;

    Master::Master<dimension + 1>& master;
    Shape::Shape<dimension + 1>& shape;

    Array2D<double> surface_normal;

  private:
    Array2D<double> lambda_gp;
    Array2D<double> int_lambda_fact;

  public:
    template <typename BoundaryType>
    EdgeBoundary(const BoundaryType& bound);

    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);
    double IntegrationLambda(const uint dof, const std::vector<double>& u_gp);
};

template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
template <typename BoundaryType>
EdgeBoundary<dimension, BasisType, DataType, EdgeDataType>::EdgeBoundary(const BoundaryType& bound)
    : data(bound.data), bound_id(bound.bound_id), master(bound.master), shape(bound.shape) {
    // *** //
    typename BoundaryType::BoundaryIntegrationType integration;

    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule =
        integration.GetRule(2 * this->master.p + 1);

    BasisType basis;

    this->lambda_gp = basis.GetPhi(this->master.p, integration_rule.second);

    std::vector<Point<dimension + 1>> z_master =
        this->master.BoundaryToMasterCoordinates(this->bound_id, integration_rule.second);

    std::vector<double> surface_J = this->shape.GetSurfaceJ(this->bound_id, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_lambda_fact = this->lambda_gp;
        for (uint dof = 0; dof < this->int_lambda_fact.size(); dof++) {
            for (uint gp = 0; gp < this->int_lambda_fact[dof].size(); gp++) {
                this->int_lambda_fact[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->surface_normal = Array2D<double>(integration_rule.first.size(),
                                               *this->shape.GetSurfaceNormal(this->bound_id, z_master).begin());
    }

    uint ndof = this->lambda_gp.size();
    uint ngp  = integration_rule.first.size();

    this->edge_data = EdgeDataType(ndof, ngp);
}

template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
void EdgeBoundary<dimension, BasisType, DataType, EdgeDataType>::ComputeUgp(const std::vector<double>& u,
                                                                            std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->lambda_gp[dof][gp];
        }
    }
}

template <uint dimension, typename BasisType, typename DataType, typename EdgeDataType>
double EdgeBoundary<dimension, BasisType, DataType, EdgeDataType>::IntegrationLambda(const uint dof,
                                                                                     const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_lambda_fact[dof][gp];
    }

    return integral;
}
}

#endif