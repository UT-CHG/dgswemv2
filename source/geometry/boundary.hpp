#ifndef CLASS_BOUNDARY_HPP
#define CLASS_BOUNDARY_HPP

namespace Geometry {
template <uint dimension, class DataType>
class RawBoundary {
  public:
    uint p;
    uint nbound;

    DataType& data;

    Basis::Basis<dimension + 1>& basis;
    Master::Master<dimension + 1>& master;
    Shape::Shape<dimension + 1>& shape;

    RawBoundary(uint p,
                uint nbound,
                DataType& data,
                Basis::Basis<dimension + 1>& basis,
                Master::Master<dimension + 1>& master,
                Shape::Shape<dimension + 1>& shape)
        : p(p), nbound(nbound), data(data), basis(basis), master(master), shape(shape) {}
};

template <uint dimension, class IntegrationType, class DataType, class BoundaryType>
class Boundary {
  public:
    BoundaryType boundary_condition;

    uint nbound;
    DataType& data;

    Array2D<double> surface_normal;

  private:
    Array2D<double> phi_gp;
    Array2D<double> int_fact_phi;

  public:
    Boundary(const RawBoundary<dimension, DataType>&, BoundaryType boundary_condition = BoundaryType());

    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);
    double IntegrationPhi(uint, const std::vector<double>&);

  public:
    using BoundaryIntegrationType = IntegrationType;
};

template <uint dimension, class IntegrationType, class DataType, class BoundaryType>
Boundary<dimension, IntegrationType, DataType, BoundaryType>::Boundary(
    const RawBoundary<dimension, DataType>& raw_boundary,
    BoundaryType boundary_condition)
    : boundary_condition(std::move(boundary_condition)), data(raw_boundary.data) {
    IntegrationType integration;
    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule =
        integration.GetRule(2 * raw_boundary.p);

    std::vector<Point<dimension + 1>> z_master =
        raw_boundary.master.BoundaryToMasterCoordinates(raw_boundary.nbound, integration_rule.second);

    this->phi_gp = raw_boundary.basis.GetPhi(raw_boundary.p, z_master);

    std::vector<double> surface_J = raw_boundary.shape.GetSurfaceJ(raw_boundary.nbound, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_fact_phi = this->phi_gp;
        for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {
            for (uint gp = 0; gp < this->int_fact_phi[dof].size(); gp++) {
                this->int_fact_phi[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->surface_normal = Array2D<double>(
            integration_rule.first.size(), *raw_boundary.shape.GetSurfaceNormal(raw_boundary.nbound, z_master).begin());
    }

    this->nbound = raw_boundary.nbound;
    this->data.set_ngp_boundary(raw_boundary.nbound, integration_rule.first.size());
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
inline double Boundary<dimension, IntegrationType, DataType, BoundaryType>::IntegrationPhi(
    uint phi_n,
    const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_phi[phi_n][gp];
    }

    return integral;
}

template <uint dimension, class IntegrationType, class DataType>
class Interface {
  public:
    uint nbound_in;
    uint nbound_ex;
    DataType& data_in;
    DataType& data_ex;

    Array2D<double> surface_normal;

  private:
    Array2D<double> phi_gp_in;
    Array2D<double> phi_gp_ex;
    Array2D<double> int_fact_phi_in;
    Array2D<double> int_fact_phi_ex;

  public:
    Interface(const RawBoundary<dimension, DataType>&, const RawBoundary<dimension, DataType>&);

    void ComputeUgpIN(const std::vector<double>&, std::vector<double>&);
    void ComputeUgpEX(const std::vector<double>&, std::vector<double>&);
    double IntegrationPhiIN(uint, const std::vector<double>&);
    double IntegrationPhiEX(uint, const std::vector<double>&);
};

template <uint dimension, class IntegrationType, class DataType>
Interface<dimension, IntegrationType, DataType>::Interface(const RawBoundary<dimension, DataType>& raw_boundary_in,
                                                           const RawBoundary<dimension, DataType>& raw_boundary_ex)
    : data_in(raw_boundary_in.data), data_ex(raw_boundary_ex.data) {
    uint p = std::max(raw_boundary_in.p, raw_boundary_ex.p);

    IntegrationType integration;
    std::pair<std::vector<double>, std::vector<Point<dimension>>> integration_rule = integration.GetRule(2 * p);

    std::vector<Point<dimension + 1>> z_master =
        raw_boundary_ex.master.BoundaryToMasterCoordinates(raw_boundary_ex.nbound, integration_rule.second);
    this->phi_gp_ex = raw_boundary_ex.basis.GetPhi(raw_boundary_ex.p, z_master);

    z_master = raw_boundary_in.master.BoundaryToMasterCoordinates(raw_boundary_in.nbound, integration_rule.second);
    this->phi_gp_in = raw_boundary_in.basis.GetPhi(raw_boundary_in.p, z_master);

    std::vector<double> surface_J = raw_boundary_in.shape.GetSurfaceJ(raw_boundary_in.nbound, z_master);

    if (surface_J.size() == 1) {  // constant Jacobian
        this->int_fact_phi_in = this->phi_gp_in;
        for (uint dof = 0; dof < this->int_fact_phi_in.size(); dof++) {
            for (uint gp = 0; gp < this->int_fact_phi_in[dof].size(); gp++) {
                this->int_fact_phi_in[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->int_fact_phi_ex = this->phi_gp_ex;
        for (uint dof = 0; dof < this->int_fact_phi_ex.size(); dof++) {
            for (uint gp = 0; gp < this->int_fact_phi_ex[dof].size(); gp++) {
                this->int_fact_phi_ex[dof][gp] *= integration_rule.first[gp] * surface_J[0];
            }
        }

        this->surface_normal =
            Array2D<double>(integration_rule.first.size(),
                            *raw_boundary_in.shape.GetSurfaceNormal(raw_boundary_in.nbound, z_master).begin());
    }

    this->nbound_in = raw_boundary_in.nbound;
    this->nbound_ex = raw_boundary_ex.nbound;
    this->data_in.set_ngp_boundary(raw_boundary_in.nbound, integration_rule.first.size());
    this->data_ex.set_ngp_boundary(raw_boundary_ex.nbound, integration_rule.first.size());
}

template <uint dimension, class IntegrationType, class DataType>
inline void Interface<dimension, IntegrationType, DataType>::ComputeUgpIN(const std::vector<double>& u,
                                                                          std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->phi_gp_in[dof][gp];
        }
    }
}

template <uint dimension, class IntegrationType, class DataType>
inline double Interface<dimension, IntegrationType, DataType>::IntegrationPhiIN(uint phi_n,
                                                                                const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_phi_in[phi_n][gp];
    }

    return integral;
}

template <uint dimension, class IntegrationType, class DataType>
inline void Interface<dimension, IntegrationType, DataType>::ComputeUgpEX(const std::vector<double>& u,
                                                                          std::vector<double>& u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->phi_gp_ex[dof][gp];
        }
    }
}

template <uint dimension, class IntegrationType, class DataType>
inline double Interface<dimension, IntegrationType, DataType>::IntegrationPhiEX(uint phi_n,
                                                                                const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_phi_ex[phi_n][gp];
    }

    return integral;
}
}

#endif