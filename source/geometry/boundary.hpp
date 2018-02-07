#ifndef CLASS_BOUNDARY_HPP
#define CLASS_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, class IntegrationType, class DataType, class BoundaryType>
class Boundary {
  public:
    using BoundaryIntegrationType = IntegrationType;

  public:
    BoundaryType boundary_condition;

    uint bound_id;
    DataType* data = nullptr;

    Array2D<double> surface_normal;

  private:
    Array2D<double> phi_gp;
    std::vector<double> int_fact;
    Array2D<double> int_fact_phi;

  public:
    Boundary(const RawBoundary<dimension, DataType>& raw_boundary,
             const BoundaryType& boundary_condition = BoundaryType());

    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);
    double Integration(const std::vector<double>& u_gp);
    double IntegrationPhi(const uint dof, const std::vector<double>& u_gp);

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned);
#endif
};

template <uint dimension, class IntegrationType, class DataType, class BoundaryType>
Boundary<dimension, IntegrationType, DataType, BoundaryType>::Boundary(
    const RawBoundary<dimension, DataType>& raw_boundary,
    const BoundaryType& boundary_condition)
    : boundary_condition(std::move(boundary_condition)), data(&raw_boundary.data) {
    IntegrationType integration;
    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule =
        integration.GetRule(2 * raw_boundary.p);

    std::vector<Point<dimension + 1>> z_master =
        raw_boundary.master.BoundaryToMasterCoordinates(raw_boundary.bound_id, integration_rule.second);

    this->phi_gp = raw_boundary.basis.GetPhi(raw_boundary.p, z_master);

    std::vector<double> surface_J = raw_boundary.shape.GetSurfaceJ(raw_boundary.bound_id, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_fact = integration_rule.first;
        for (uint gp = 0; gp < this->int_fact.size(); gp++) {
            this->int_fact[gp] *= surface_J[0];
        }

        this->int_fact_phi = this->phi_gp;
        for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {
            for (uint gp = 0; gp < this->int_fact_phi[dof].size(); gp++) {
                this->int_fact_phi[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->surface_normal =
            Array2D<double>(integration_rule.first.size(),
                            *raw_boundary.shape.GetSurfaceNormal(raw_boundary.bound_id, z_master).begin());
    }

    this->bound_id = raw_boundary.bound_id;
    this->data->set_ngp_boundary(raw_boundary.bound_id, integration_rule.first.size());
}

template <uint dimension, class IntegrationType, class DataType, class BoundaryType>
inline void Boundary<dimension, IntegrationType, DataType, BoundaryType>::ComputeUgp(const std::vector<double>& u,
                                                                                     std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->phi_gp[dof][gp];
        }
    }
}

template <uint dimension, class IntegrationType, class DataType, class BoundaryType>
inline double Boundary<dimension, IntegrationType, DataType, BoundaryType>::Integration(
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact[gp];
    }

    return integral;
}

template <uint dimension, class IntegrationType, class DataType, class BoundaryType>
inline double Boundary<dimension, IntegrationType, DataType, BoundaryType>::IntegrationPhi(
    const uint dof,
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_phi[dof][gp];
    }

    return integral;
}

#ifdef HAS_HPX
template <uint dimension, class IntegrationType, class DataType, class BoundaryType>
template<typename Archive>
void Boundary<dimension, IntegrationType, DataType, BoundaryType>::serialize(Archive& ar, unsigned) {
    ar & boundary_condition
       & bound_id
       & surface_normal
       & phi_gp
       & int_fact
       & int_fact_phi;
}
#endif
}
#endif
