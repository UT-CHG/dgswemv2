#ifndef CLASS_ELEMENT_HPP
#define CLASS_ELEMENT_HPP

namespace Geometry {
template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
class Element {
  public:
    DataType data;

  private:
    uint ID;

    MasterType& master;
    ShapeType   shape;

    std::vector<uint>  node_ID;
    std::vector<uint>  neighbor_ID;
    std::vector<uchar> boundary_type;

    std::vector<Point<dimension>> gp_global_coordinates;

    Array3D<double> dpsi_fact;
    Array3D<double> dphi_fact;

    bool                const_J;
    std::vector<double> int_fact;
    Array2D<double>     int_fact_phi;
    Array3D<double>     int_fact_dphi;

    std::pair<bool, Array2D<double>> m_inv;

  public:
    Element(const uint                           ID,
            MasterType&                          master,
            const std::vector<Point<dimension>>& nodal_coordinates,
            const std::vector<uint>&             node_ID,
            const std::vector<uint>&             neighbor_ID,
            const std::vector<uchar>&            boundary_type);

    void CreateRawBoundaries(
        std::map<uint, std::map<uint, RawBoundary<dimension - 1, DataType>>>& pre_interfaces,
        std::map<uchar, std::vector<RawBoundary<dimension - 1, DataType>>>&   pre_boundaries,
        std::map<uint, std::map<uint, RawBoundary<dimension - 1, DataType>>>& pre_distributed_boundaries);

    uint        GetID() { return this->ID; }
    MasterType& GetMaster() { return this->master; }
    ShapeType&  GetShape() { return this->shape; }

    std::vector<uint>&  GetNodeID() { return this->node_ID; }
    std::vector<uchar>& GetBoundaryType() { return this->boundary_type; }

    template <typename F>
    void L2Projection(const F& f, std::vector<double>& projection);
    void L2Projection(const std::vector<double>& nodal_values, std::vector<double>& projection);

    void ProjectBasisToLinear(const std::vector<double>& u, std::vector<double>& u_lin);
    void ProjectLinearToBasis(const std::vector<double>& u_lin, std::vector<double>& u);

    template <typename F>
    void ComputeFgp(const F& f, std::vector<double>& f_gp);
    void ComputeUgp(const std::vector<double>& u, std::vector<double>& u_gp);
    void ComputeDUgp(const uint dir, const std::vector<double>& u, std::vector<double>& du_gp);

    void ComputeLinearUgp(const std::vector<double>& u_lin, std::vector<double>& u_lin_gp);
    void ComputeLinearDUgp(const uint dir, const std::vector<double>& u_lin, std::vector<double>& du_lin_gp);
    void ComputeLinearUbaryctr(const std::vector<double>& u_lin, double& u_lin_baryctr);
    void ComputeLinearUmidpts(const std::vector<double>& u_lin, std::vector<double>& u_lin_midpts);
    void ComputeLinearUvrtx(const std::vector<double>& u_lin, std::vector<double>& u_lin_vrtx);

    double Integration(const std::vector<double>& u_gp);
    double IntegrationPhi(const uint dof, const std::vector<double>& u_gp);
    double IntegrationDPhi(const uint dir, const uint dof, const std::vector<double>& u_gp);

    void ApplyMinv(const std::vector<double>& rhs, std::vector<double>& solution);

    void InitializeVTK(std::vector<Point<3>>& points, Array2D<uint>& cells);
    void WriteCellDataVTK(const std::vector<double>& u, std::vector<double>& cell_data);
    void WritePointDataVTK(const std::vector<double>& u, std::vector<double>& point_data);

    template <typename F>
    double ComputeResidualL2(const F& f, const std::vector<double>& u);

  public:
    using ElementMasterType = MasterType;
};

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
Element<dimension, MasterType, ShapeType, DataType>::Element(const uint                           ID,
                                                             MasterType&                          master,
                                                             const std::vector<Point<dimension>>& nodal_coordinates,
                                                             const std::vector<uint>&             node_ID,
                                                             const std::vector<uint>&             neighbor_ID,
                                                             const std::vector<uchar>&            boundary_type)
    : ID(ID),
      master(master),
      shape(ShapeType(nodal_coordinates)),
      node_ID(std::move(node_ID)),
      neighbor_ID(std::move(neighbor_ID)),
      boundary_type(std::move(boundary_type)) {
    // GLOBAL COORDINATES OF GPS
    this->gp_global_coordinates = this->shape.LocalToGlobalCoordinates(this->master.integration_rule.second);

    // DEFORMATION
    std::vector<double> det_J = this->shape.GetJdet(this->master.integration_rule.second);
    Array3D<double>     J_inv = this->shape.GetJinv(this->master.integration_rule.second);

    this->const_J = (det_J.size() == 1);

    if (const_J) {  // constant Jacobian
        // DIFFERENTIATION FACTORS
        this->dpsi_fact.resize(this->master.dpsi.size());
        for (uint dof = 0; dof < this->master.dpsi.size(); dof++) {
            this->dpsi_fact[dof].resize(dimension);
            for (uint dir = 0; dir < dimension; dir++) {
                this->dpsi_fact[dof][dir].reserve(this->master.dphi_gp[dof][dir].size());
                for (uint gp = 0; gp < this->master.dphi_gp[dof][dir].size(); gp++) {
                    double dpsi = 0;
                    for (uint z = 0; z < dimension; z++) {
                        dpsi += this->master.dpsi[dof][z] * J_inv[z][dir][0];
                    }
                    this->dpsi_fact[dof][dir].push_back(dpsi);
                }
            }
        }

        this->dphi_fact.resize(this->master.dphi_gp.size());
        for (uint dof = 0; dof < this->master.dphi_gp.size(); dof++) {
            this->dphi_fact[dof].resize(dimension);
            for (uint dir = 0; dir < dimension; dir++) {
                this->dphi_fact[dof][dir].reserve(this->master.dphi_gp[dof][dir].size());
                for (uint gp = 0; gp < this->master.dphi_gp[dof][dir].size(); gp++) {
                    double dphi = 0;
                    for (uint z = 0; z < dimension; z++) {
                        dphi += this->master.dphi_gp[dof][z][gp] * J_inv[z][dir][0];
                    }
                    this->dphi_fact[dof][dir].push_back(dphi);
                }
            }
        }

        // INTEGRATION OVER ELEMENT FACTORS
        this->int_fact = this->master.integration_rule.first;
        for (uint gp = 0; gp < this->int_fact.size(); gp++) {
            this->int_fact[gp] *= std::abs(det_J[0]);
        }

        this->int_fact_phi = this->master.int_fact_phi;
        for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {
            for (uint gp = 0; gp < this->int_fact_phi[dof].size(); gp++) {
                this->int_fact_phi[dof][gp] *= std::abs(det_J[0]);
            }
        }

        this->int_fact_dphi.resize(this->master.int_fact_dphi.size());
        for (uint dof = 0; dof < this->master.int_fact_dphi.size(); dof++) {
            this->int_fact_dphi[dof].resize(dimension);
            for (uint dir = 0; dir < dimension; dir++) {
                this->int_fact_dphi[dof][dir].reserve(this->master.int_fact_dphi[dof][dir].size());
                for (uint gp = 0; gp < this->master.int_fact_dphi[dof][dir].size(); gp++) {
                    double int_dphi = 0;
                    for (uint z = 0; z < dimension; z++) {
                        int_dphi += this->master.int_fact_dphi[dof][z][gp] * J_inv[z][dir][0];
                    }
                    int_dphi *= std::abs(det_J[0]);
                    this->int_fact_dphi[dof][dir].push_back(int_dphi);
                }
            }
        }

        // MASS MATRIX
        this->m_inv = this->master.m_inv;
        for (uint i = 0; i < this->m_inv.second.size(); i++) {
            for (uint j = 0; j < this->m_inv.second[i].size(); j++) {
                this->m_inv.second[i][j] /= std::abs(det_J[0]);
            }
        }
    } else {
        // Placeholder for nonconstant Jacobian
    }

    this->data.set_nvrtx(this->master.nvrtx);
    this->data.set_nbound(this->master.nbound);
    this->data.set_ndof(this->master.phi_gp.size());
    this->data.set_ngp_internal((*this->master.phi_gp.begin()).size());
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
void Element<dimension, MasterType, ShapeType, DataType>::CreateRawBoundaries(
    std::map<uint, std::map<uint, RawBoundary<dimension - 1, DataType>>>& pre_interfaces,
    std::map<uchar, std::vector<RawBoundary<dimension - 1, DataType>>>&   pre_boundaries,
    std::map<uint, std::map<uint, RawBoundary<dimension - 1, DataType>>>& pre_distributed_boundaries) {
    //Separate from Args
    Basis::Basis<dimension>*   my_basis  = (Basis::Basis<dimension>*)(&this->master.basis);
    Master::Master<dimension>* my_master = (Master::Master<dimension>*)(&this->master);
    Shape::Shape<dimension>*   my_shape  = (Shape::Shape<dimension>*)(&this->shape);

    for (uint bound_id = 0; bound_id < this->boundary_type.size(); bound_id++) {
        if (this->boundary_type[bound_id] == INTERNAL) {
            pre_interfaces[this->ID].emplace(
                std::make_pair(this->neighbor_ID[bound_id],
                               RawBoundary<dimension - 1, DataType>(
                                   this->master.p, bound_id, this->data, *my_basis, *my_master, *my_shape)));
        } else if (this->boundary_type[bound_id] == DISTRIBUTED) {
            pre_distributed_boundaries[this->ID].emplace(
                std::make_pair(bound_id,
                               RawBoundary<dimension - 1, DataType>(
                                   this->master.p, bound_id, this->data, *my_basis, *my_master, *my_shape)));
        } else {
            pre_boundaries[this->boundary_type[bound_id]].emplace_back(RawBoundary<dimension - 1, DataType>(
                this->master.p, bound_id, this->data, *my_basis, *my_master, *my_shape));
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F>
inline void Element<dimension, MasterType, ShapeType, DataType>::L2Projection(const F&             f,
                                                                              std::vector<double>& projection) {
    std::vector<double> rhs;

    std::vector<double> f_vals(this->gp_global_coordinates.size());

    this->ComputeFgp(f, f_vals);

    for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {
        rhs.push_back(this->IntegrationPhi(dof, f_vals));
    }

    this->ApplyMinv(rhs, projection);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::L2Projection(const std::vector<double>& nodal_values,
                                                                              std::vector<double>&       projection) {
    std::vector<double> rhs;

    std::vector<double> interpolation =
        this->shape.InterpolateNodalValues(nodal_values, this->master.integration_rule.second);

    for (uint dof = 0; dof < this->int_fact_phi.size(); dof++) {
        rhs.push_back(this->IntegrationPhi(dof, interpolation));
    }

    this->ApplyMinv(rhs, projection);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ProjectBasisToLinear(const std::vector<double>& u,
                                                                                      std::vector<double>& u_lin) {
    if (const_J) {
        this->master.basis.ProjectBasisToLinear(u, u_lin);
    } else {
        // Placeholder for nonconstant Jacobian
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ProjectLinearToBasis(const std::vector<double>& u_lin,
                                                                                      std::vector<double>&       u) {
    if (const_J) {
        this->master.basis.ProjectLinearToBasis(u_lin, u);
    } else {
        // Placeholder for nonconstant Jacobian
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeFgp(const F& f, std::vector<double>& f_gp) {
    for (uint gp = 0; gp < f_gp.size(); gp++) {
        f_gp[gp] = f(this->gp_global_coordinates[gp]);
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeUgp(const std::vector<double>& u,
                                                                            std::vector<double>&       u_gp) {
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->master.phi_gp[dof][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeDUgp(const uint                 dir,
                                                                             const std::vector<double>& u,
                                                                             std::vector<double>&       du_gp) {
    std::fill(du_gp.begin(), du_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < du_gp.size(); gp++) {
            du_gp[gp] += u[dof] * this->dphi_fact[dof][dir][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUgp(const std::vector<double>& u_lin,
                                                                                  std::vector<double>&       u_lin_gp) {
    std::fill(u_lin_gp.begin(), u_lin_gp.end(), 0.0);

    for (uint dof = 0; dof < u_lin.size(); dof++) {
        for (uint gp = 0; gp < u_lin_gp.size(); gp++) {
            u_lin_gp[gp] += u_lin[dof] * this->master.psi_gp[dof][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearDUgp(const uint                 dir,
                                                                                   const std::vector<double>& u_lin,
                                                                                   std::vector<double>& du_lin_gp) {
    std::fill(du_lin_gp.begin(), du_lin_gp.end(), 0.0);

    for (uint dof = 0; dof < u_lin.size(); dof++) {
        for (uint gp = 0; gp < du_lin_gp.size(); gp++) {
            du_lin_gp[gp] += u_lin[dof] * this->dpsi_fact[dof][dir][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUbaryctr(const std::vector<double>& u_lin,
                                                                                       double& u_lin_baryctr) {
    this->master.ComputeLinearUbaryctr(u_lin, u_lin_baryctr);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUmidpts(
    const std::vector<double>& u_lin,
    std::vector<double>&       u_lin_midpts) {
    this->master.ComputeLinearUmidpts(u_lin, u_lin_midpts);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUvrtx(const std::vector<double>& u_lin,
                                                                                    std::vector<double>& u_lin_vrtx) {
    this->master.ComputeLinearUvrtx(u_lin, u_lin_vrtx);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline double Element<dimension, MasterType, ShapeType, DataType>::Integration(const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact[gp];
    }

    return integral;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline double Element<dimension, MasterType, ShapeType, DataType>::IntegrationPhi(const uint                 dof,
                                                                                  const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_phi[dof][gp];
    }

    return integral;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline double Element<dimension, MasterType, ShapeType, DataType>::IntegrationDPhi(const uint                 dir,
                                                                                   const uint                 dof,
                                                                                   const std::vector<double>& u_gp) {
    double integral = 0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact_dphi[dof][dir][gp];
    }

    return integral;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::ApplyMinv(const std::vector<double>& rhs,
                                                                           std::vector<double>&       solution) {
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

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
void Element<dimension, MasterType, ShapeType, DataType>::InitializeVTK(std::vector<Point<3>>& points,
                                                                        Array2D<uint>&         cells) {
    this->shape.GetVTK(points, cells);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::WriteCellDataVTK(const std::vector<double>& u,
                                                                                  std::vector<double>& cell_data) {
    Array2D<double> temp = this->master.phi_postprocessor_cell;

    for (uint dof = 0; dof < temp.size(); dof++) {
        for (uint cell = 0; cell < temp[dof].size(); cell++) {
            temp[dof][cell] *= u[dof];
        }
    }

    for (uint dof = 1; dof < temp.size(); dof++) {
        for (uint cell = 0; cell < temp[dof].size(); cell++) {
            temp[0][cell] += temp[dof][cell];
        }
    }

    cell_data.insert(cell_data.end(), temp[0].begin(), temp[0].end());
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
inline void Element<dimension, MasterType, ShapeType, DataType>::WritePointDataVTK(const std::vector<double>& u,
                                                                                   std::vector<double>& point_data) {
    Array2D<double> temp = this->master.phi_postprocessor_point;

    for (uint dof = 0; dof < temp.size(); dof++) {
        for (uint pt = 0; pt < temp[dof].size(); pt++) {
            temp[dof][pt] *= u[dof];
        }
    }

    for (uint dof = 1; dof < temp.size(); dof++) {
        for (uint pt = 0; pt < temp[dof].size(); pt++) {
            temp[0][pt] += temp[dof][pt];
        }
    }

    point_data.insert(point_data.end(), temp[0].begin(), temp[0].end());
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F>
double Element<dimension, MasterType, ShapeType, DataType>::ComputeResidualL2(const F&                   f,
                                                                              const std::vector<double>& u) {
    std::pair<std::vector<double>, std::vector<Point<2>>> rule = this->master.integration.GetRule(20);
    // At this point we use maximum possible p for Dunavant integration

    Array2D<double> Phi = this->master.basis.GetPhi(this->master.p, rule.second);

    std::vector<double> u_gp(rule.first.size());
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < this->data.get_ndof(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += Phi[dof][gp] * u[dof];
        }
    }

    std::vector<Point<2>> gp_global = this->shape.LocalToGlobalCoordinates(rule.second);

    std::vector<double> f_gp(rule.first.size());

    for (uint gp = 0; gp < f_gp.size(); gp++) {
        f_gp[gp] = f(gp_global[gp]);
    }

    std::vector<double> sq_diff(rule.first.size());

    for (uint gp = 0; gp < sq_diff.size(); gp++) {
        sq_diff[gp] = pow((f_gp[gp] - u_gp[gp]), 2);
    }

    double L2 = 0;

    if (const_J) {
        for (uint gp = 0; gp < sq_diff.size(); gp++) {
            L2 += sq_diff[gp] * rule.first[gp];
        }

        L2 *= std::abs(this->shape.GetJdet(rule.second)[0]);
    } else {
        // Placeholder for nonconstant Jacobian
    }

    return L2;
}
}

#endif
