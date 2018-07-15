#ifndef CLASS_ELEMENT_HPP
#define CLASS_ELEMENT_HPP

namespace Geometry {
template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
class Element {
  public:
    DataType data;

  private:
    uint ID;

    MasterType* master = nullptr;
    ShapeType shape;

    std::vector<uint> node_ID;
    std::vector<uint> neighbor_ID;
    std::vector<uchar> boundary_type;

    std::vector<Point<dimension>> gp_global_coordinates;

    /* psi_gp stored in shape */   // nodal basis, i.e. shape functions
    /* chi_gp stored in master */  // linear basis
    /* phi_gp stroed in master */  // modal basis

    /* dpsi_gp stored in shape */  // nodal basis, i.e. shape functions
    Array3D<double> dchi_gp;       // linear basis
    Array3D<double> dphi_gp;       // modal basis

    bool const_J;
    std::vector<double> int_fact;
    Array2D<double> int_phi_fact;
    Array3D<double> int_phi_phi_fact;
    Array3D<double> int_dphi_fact;
    Array4D<double> int_phi_dphi_fact;

    std::pair<bool, Array2D<double>> m_inv;

  public:
    Element() = default;
    Element(const uint ID,
            MasterType& master,
            std::vector<Point<3>>&& nodal_coordinates,
            std::vector<uint>&& node_ID,
            std::vector<uint>&& neighbor_ID,
            std::vector<uchar>&& boundary_type);

    uint GetID() { return this->ID; }
    MasterType& GetMaster() { return *this->master; }
    ShapeType& GetShape() { return this->shape; }
    std::vector<uint>& GetNodeID() { return this->node_ID; }
    std::vector<uchar>& GetBoundaryType() { return this->boundary_type; }

    void SetMaster(MasterType& master) { this->master = &master; };

    void Initialize();
    void CreateRawBoundaries(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundary<dimension - 1, DataType>>>&
                                 pre_specialized_interfaces);

    template <typename F, typename T>
    void L2Projection(const F& f, std::vector<T>& projection);
    template <typename T>
    void L2Projection(const std::vector<T>& nodal_values, std::vector<T>& projection);

    template <typename T>
    void ProjectBasisToLinear(const std::vector<T>& u, std::vector<T>& u_lin);
    template <typename T>
    void ProjectLinearToBasis(const std::vector<T>& u_lin, std::vector<T>& u);

    template <typename F, typename T>
    void ComputeFgp(const F& f, std::vector<T>& f_gp);
    template <typename T>
    void ComputeUgp(const std::vector<T>& u, std::vector<T>& u_gp);
    template <typename T>
    void ComputeDUgp(const uint dir, const std::vector<T>& u, std::vector<T>& du_gp);
    template <typename T>
    void ComputeDUgp(const std::vector<T>& u, std::vector<StatVector<T, dimension>>& du_gp);

    template <typename T>
    void ComputeLinearUgp(const std::vector<T>& u_lin, std::vector<T>& u_lin_gp);
    template <typename T>
    void ComputeLinearDUgp(const uint dir, const std::vector<T>& u_lin, std::vector<T>& du_lin_gp);
    template <typename T>
    void ComputeLinearUbaryctr(const std::vector<T>& u_lin, T& u_lin_baryctr);
    template <typename T>
    void ComputeLinearUmidpts(const std::vector<T>& u_lin, std::vector<T>& u_lin_midpts);
    template <typename T>
    void ComputeLinearUvrtx(const std::vector<T>& u_lin, std::vector<T>& u_lin_vrtx);

    template <typename T>
    void ComputeNodalUgp(const std::vector<T>& u_nodal, std::vector<T>& u_nodal_gp);
    template <typename T>
    void ComputeNodalDUgp(const uint dir, const std::vector<T>& u_nodal, std::vector<T>& du_nodal_gp);
    template <typename T>
    void ComputeNodalDUgp(const std::vector<T>& u_nodal, std::vector<StatVector<T, dimension>>& du_nodal_gp);

    template <typename T>
    T Integration(const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhi(const uint dof, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhiPhi(const uint dof_i, const uint dof_j, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationDPhi(const uint dir, const uint dof, const std::vector<T>& u_gp);
    template <typename T>
    T IntegrationPhiDPhi(const uint dof_i, const uint dir_j, const uint dof_j, const std::vector<T>& u_gp);

    template <typename T>
    void ApplyMinv(const std::vector<T>& rhs, std::vector<T>& solution);

    void InitializeVTK(std::vector<Point<3>>& points, Array2D<uint>& cells);
    template <typename T>
    void WriteCellDataVTK(const std::vector<T>& u, std::vector<T>& cell_data);
    template <typename T>
    void WritePointDataVTK(const std::vector<T>& u, std::vector<T>& point_data);

    template <typename F, typename T>
    T ComputeResidualL2(const F& f, const std::vector<T>& u);

  public:
    using ElementMasterType = MasterType;

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & data  
            & ID
            & shape
            & node_ID
            & neighbor_ID
            & boundary_type;
        // clang-format on
    }
#endif
};

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
Element<dimension, MasterType, ShapeType, DataType>::Element(const uint ID,
                                                             MasterType& master,
                                                             std::vector<Point<3>>&& nodal_coordinates,
                                                             std::vector<uint>&& node_ID,
                                                             std::vector<uint>&& neighbor_ID,
                                                             std::vector<uchar>&& boundary_type)
    : ID(ID),
      master(&master),
      shape(ShapeType(std::move(nodal_coordinates))),
      node_ID(std::move(node_ID)),
      neighbor_ID(std::move(neighbor_ID)),
      boundary_type(std::move(boundary_type)) {
    this->Initialize();
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
void Element<dimension, MasterType, ShapeType, DataType>::Initialize() {
    // GLOBAL COORDINATES OF GPS
    this->gp_global_coordinates = this->shape.LocalToGlobalCoordinates(this->master->integration_rule.second);

    // DEFORMATION
    std::vector<double> det_J = this->shape.GetJdet(this->master->integration_rule.second);
    Array3D<double> J_inv     = this->shape.GetJinv(this->master->integration_rule.second);

    this->const_J = (det_J.size() == 1);

    // Compute factors to expand nodal values and derivatives of nodal values
    this->shape.psi_gp  = this->shape.GetPsi(this->master->integration_rule.second);
    this->shape.dpsi_gp = this->shape.GetDPsi(this->master->integration_rule.second);

    if (const_J) {  // constant Jacobian
        // DIFFERENTIATION FACTORS
        this->dchi_gp.resize(this->master->dchi_gp.size());
        for (uint dof = 0; dof < this->master->dchi_gp.size(); dof++) {
            this->dchi_gp[dof].resize(dimension);
            for (uint dir = 0; dir < dimension; dir++) {
                this->dchi_gp[dof][dir].reserve(this->master->dphi_gp[dof][dir].size());
                for (uint gp = 0; gp < this->master->dphi_gp[dof][dir].size(); gp++) {
                    double dchi = 0;
                    for (uint z = 0; z < dimension; z++) {
                        dchi += this->master->dchi_gp[dof][z] * J_inv[z][dir][0];
                    }
                    this->dchi_gp[dof][dir].push_back(dchi);
                }
            }
        }

        this->dphi_gp.resize(this->master->dphi_gp.size());
        for (uint dof = 0; dof < this->master->dphi_gp.size(); dof++) {
            this->dphi_gp[dof].resize(dimension);
            for (uint dir = 0; dir < dimension; dir++) {
                this->dphi_gp[dof][dir].reserve(this->master->dphi_gp[dof][dir].size());
                for (uint gp = 0; gp < this->master->dphi_gp[dof][dir].size(); gp++) {
                    double dphi = 0;
                    for (uint z = 0; z < dimension; z++) {
                        dphi += this->master->dphi_gp[dof][z][gp] * J_inv[z][dir][0];
                    }
                    this->dphi_gp[dof][dir].push_back(dphi);
                }
            }
        }

        // INTEGRATION OVER ELEMENT FACTORS
        this->int_fact = this->master->integration_rule.first;
        for (uint gp = 0; gp < this->int_fact.size(); gp++) {
            this->int_fact[gp] *= std::abs(det_J[0]);
        }

        this->int_phi_fact = this->master->int_phi_fact;
        for (uint dof = 0; dof < this->int_phi_fact.size(); dof++) {
            for (uint gp = 0; gp < this->int_phi_fact[dof].size(); gp++) {
                this->int_phi_fact[dof][gp] *= std::abs(det_J[0]);
            }
        }

        this->int_phi_phi_fact.resize(this->master->phi_gp.size());
        for (uint dof_i = 0; dof_i < this->master->phi_gp.size(); dof_i++) {
            this->int_phi_phi_fact[dof_i] = this->int_phi_fact;
            for (uint dof_j = 0; dof_j < this->master->phi_gp.size(); dof_j++) {
                for (uint gp = 0; gp < this->int_phi_phi_fact[dof_i][dof_j].size(); gp++) {
                    this->int_phi_phi_fact[dof_i][dof_j][gp] *= this->master->phi_gp[dof_i][gp];
                }
            }
        }

        this->int_dphi_fact.resize(this->master->int_dphi_fact.size());
        for (uint dof = 0; dof < this->master->int_dphi_fact.size(); dof++) {
            this->int_dphi_fact[dof].resize(dimension);
            for (uint dir = 0; dir < dimension; dir++) {
                this->int_dphi_fact[dof][dir].reserve(this->master->int_dphi_fact[dof][dir].size());
                for (uint gp = 0; gp < this->master->int_dphi_fact[dof][dir].size(); gp++) {
                    double int_dphi = 0;
                    for (uint z = 0; z < dimension; z++) {
                        int_dphi += this->master->int_dphi_fact[dof][z][gp] * J_inv[z][dir][0];
                    }
                    int_dphi *= std::abs(det_J[0]);
                    this->int_dphi_fact[dof][dir].push_back(int_dphi);
                }
            }
        }

        this->int_phi_dphi_fact.resize(this->master->phi_gp.size());
        for (uint dof_i = 0; dof_i < this->master->phi_gp.size(); dof_i++) {
            this->int_phi_dphi_fact[dof_i] = this->int_dphi_fact;
            for (uint dof_j = 0; dof_j < this->master->phi_gp.size(); dof_j++) {
                for (uint dir = 0; dir < dimension; dir++) {
                    for (uint gp = 0; gp < this->int_phi_dphi_fact[dof_i][dof_j][dir].size(); gp++) {
                        this->int_phi_dphi_fact[dof_i][dof_j][dir][gp] *= this->master->phi_gp[dof_i][gp];
                    }
                }
            }
        }

        // MASS MATRIX
        this->m_inv = this->master->m_inv;
        for (uint i = 0; i < this->m_inv.second.size(); i++) {
            for (uint j = 0; j < this->m_inv.second[i].size(); j++) {
                this->m_inv.second[i][j] /= std::abs(det_J[0]);
            }
        }
    } else {
        // Placeholder for nonconstant Jacobian
    }

    this->data.set_nnode(this->shape.nodal_coordinates.size());
    this->data.set_nvrtx(this->master->nvrtx);
    this->data.set_nbound(this->master->nbound);
    this->data.set_ndof(this->master->phi_gp.size());
    this->data.set_ngp_internal((*this->master->phi_gp.begin()).size());
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
void Element<dimension, MasterType, ShapeType, DataType>::CreateRawBoundaries(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundary<dimension - 1, DataType>>>& raw_boundaries) {
    // *** //
    Basis::Basis<dimension>* my_basis    = (Basis::Basis<dimension>*)(&this->master->basis);
    Master::Master<dimension>* my_master = (Master::Master<dimension>*)(this->master);
    Shape::Shape<dimension>* my_shape    = (Shape::Shape<dimension>*)(&this->shape);

    for (uint bound_id = 0; bound_id < this->boundary_type.size(); bound_id++) {
        std::vector<uint> bound_node_ID = this->shape.GetBoundaryNodeID(bound_id, this->node_ID);

        if (is_internal(this->boundary_type[bound_id])) {
            raw_boundaries[this->boundary_type[bound_id]].emplace(
                std::pair<uint, uint>{this->ID, this->neighbor_ID[bound_id]},
                RawBoundary<dimension - 1, DataType>(
                    this->master->p, bound_id, bound_node_ID, this->data, *my_basis, *my_master, *my_shape));
        } else {
            raw_boundaries[this->boundary_type[bound_id]].emplace(
                std::pair<uint, uint>{this->ID, bound_id},
                RawBoundary<dimension - 1, DataType>(
                    this->master->p, bound_id, bound_node_ID, this->data, *my_basis, *my_master, *my_shape));
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F, typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::L2Projection(const F& f, std::vector<T>& projection) {
    std::vector<T> rhs;

    std::vector<T> f_vals(this->gp_global_coordinates.size());

    this->ComputeFgp(f, f_vals);

    for (uint dof = 0; dof < this->int_phi_fact.size(); dof++) {
        rhs.push_back(this->IntegrationPhi(dof, f_vals));
    }

    this->ApplyMinv(rhs, projection);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::L2Projection(const std::vector<T>& nodal_values,
                                                                              std::vector<T>& projection) {
    std::vector<T> rhs;

    std::vector<T> interpolation(this->data.get_ngp_internal());

    this->ComputeNodalUgp(nodal_values, interpolation);

    for (uint dof = 0; dof < this->int_phi_fact.size(); dof++) {
        rhs.push_back(this->IntegrationPhi(dof, interpolation));
    }

    this->ApplyMinv(rhs, projection);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ProjectBasisToLinear(const std::vector<T>& u,
                                                                                      std::vector<T>& u_lin) {
    if (const_J) {
        this->master->basis.ProjectBasisToLinear(u, u_lin);
    } else {
        // Placeholder for nonconstant Jacobian
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ProjectLinearToBasis(const std::vector<T>& u_lin,
                                                                                      std::vector<T>& u) {
    if (const_J) {
        this->master->basis.ProjectLinearToBasis(u_lin, u);
    } else {
        // Placeholder for nonconstant Jacobian
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F, typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeFgp(const F& f, std::vector<T>& f_gp) {
    for (uint gp = 0; gp < f_gp.size(); gp++) {
        f_gp[gp] = f(this->gp_global_coordinates[gp]);
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeUgp(const std::vector<T>& u,
                                                                            std::vector<T>& u_gp) {
    assert(this->master);
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += u[dof] * this->master->phi_gp[dof][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeDUgp(const uint dir,
                                                                             const std::vector<T>& u,
                                                                             std::vector<T>& du_gp) {
    std::fill(du_gp.begin(), du_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < du_gp.size(); gp++) {
            du_gp[gp] += u[dof] * this->dphi_gp[dof][dir][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeDUgp(const std::vector<T>& u,
                                                                             std::vector<StatVector<T, dimension>>& du_gp) {
    std::fill(du_gp.begin(), du_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < du_gp.size(); gp++) {
            for (uint dir = 0; dir < dimension; dir++) {
                du_gp[gp][dir] += u[dof] * this->dphi_gp[dof][dir][gp];
            }
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUgp(const std::vector<T>& u_lin,
                                                                                  std::vector<T>& u_lin_gp) {
    std::fill(u_lin_gp.begin(), u_lin_gp.end(), 0.0);

    for (uint dof = 0; dof < u_lin.size(); dof++) {
        for (uint gp = 0; gp < u_lin_gp.size(); gp++) {
            u_lin_gp[gp] += u_lin[dof] * this->master->chi_gp[dof][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearDUgp(const uint dir,
                                                                                   const std::vector<T>& u_lin,
                                                                                   std::vector<T>& du_lin_gp) {
    std::fill(du_lin_gp.begin(), du_lin_gp.end(), 0.0);

    for (uint dof = 0; dof < u_lin.size(); dof++) {
        for (uint gp = 0; gp < du_lin_gp.size(); gp++) {
            du_lin_gp[gp] += u_lin[dof] * this->dchi_gp[dof][dir][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUbaryctr(const std::vector<T>& u_lin,
                                                                                       T& u_lin_baryctr) {
    this->master->ComputeLinearUbaryctr(u_lin, u_lin_baryctr);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUmidpts(const std::vector<T>& u_lin,
                                                                                      std::vector<T>& u_lin_midpts) {
    this->master->ComputeLinearUmidpts(u_lin, u_lin_midpts);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUvrtx(const std::vector<T>& u_lin,
                                                                                    std::vector<T>& u_lin_vrtx) {
    this->master->ComputeLinearUvrtx(u_lin, u_lin_vrtx);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeNodalUgp(const std::vector<T>& u_nodal,
                                                                                 std::vector<T>& u_nodal_gp) {
    std::fill(u_nodal_gp.begin(), u_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_nodal.size(); dof++) {
        for (uint gp = 0; gp < u_nodal_gp.size(); gp++) {
            u_nodal_gp[gp] += u_nodal[dof] * this->shape.psi_gp[dof][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeNodalDUgp(const uint dir,
                                                                                  const std::vector<T>& u_nodal,
                                                                                  std::vector<T>& du_nodal_gp) {
    std::fill(du_nodal_gp.begin(), du_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_nodal.size(); dof++) {
        for (uint gp = 0; gp < du_nodal_gp.size(); gp++) {
            du_nodal_gp[gp] += u_nodal[dof] * this->shape.dpsi_gp[dof][dir][gp];
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeNodalDUgp(
    const std::vector<T>& u_nodal,
    std::vector<StatVector<T, dimension>>& du_nodal_gp) {
    // *** //
    std::fill(du_nodal_gp.begin(), du_nodal_gp.end(), 0.0);

    for (uint dof = 0; dof < u_nodal.size(); dof++) {
        for (uint gp = 0; gp < du_nodal_gp.size(); gp++) {
            for (uint dir = 0; dir < dimension; dir++) {
                du_nodal_gp[gp][dir] += u_nodal[dof] * this->shape.dpsi_gp[dof][dir][gp];
            }
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline T Element<dimension, MasterType, ShapeType, DataType>::Integration(const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_fact[gp];
    }

    return integral;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline T Element<dimension, MasterType, ShapeType, DataType>::IntegrationPhi(const uint dof,
                                                                             const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_fact[dof][gp];
    }

    return integral;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline T Element<dimension, MasterType, ShapeType, DataType>::IntegrationPhiPhi(const uint dof_i,
                                                                                const uint dof_j,
                                                                                const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_phi_fact[dof_i][dof_j][gp];
    }

    return integral;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline T Element<dimension, MasterType, ShapeType, DataType>::IntegrationDPhi(const uint dir,
                                                                              const uint dof,
                                                                              const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_dphi_fact[dof][dir][gp];
    }

    return integral;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline T Element<dimension, MasterType, ShapeType, DataType>::IntegrationPhiDPhi(const uint dof_i,
                                                                                 const uint dir_j,
                                                                                 const uint dof_j,
                                                                                 const std::vector<T>& u_gp) {
    T integral;

    integral = 0.0;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_dphi_fact[dof_i][dof_j][dir_j][gp];
    }

    return integral;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ApplyMinv(const std::vector<T>& rhs,
                                                                           std::vector<T>& solution) {
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
                                                                        Array2D<uint>& cells) {
    this->shape.GetVTK(points, cells);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::WriteCellDataVTK(const std::vector<T>& u,
                                                                                  std::vector<T>& cell_data) {
    T temp;

    for (uint cell = 0; cell < this->master->phi_postprocessor_cell[0].size(); cell++) {
        temp = 0.0;

        for (uint dof = 0; dof < u.size(); dof++) {
            temp += u[dof] * this->master->phi_postprocessor_cell[dof][cell];
        }

        cell_data.push_back(temp);
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::WritePointDataVTK(const std::vector<T>& u,
                                                                                   std::vector<T>& point_data) {
    T temp;

    for (uint pt = 0; pt < this->master->phi_postprocessor_point[0].size(); pt++) {
        temp = 0.0;

        for (uint dof = 0; dof < u.size(); dof++) {
            temp += u[dof] * this->master->phi_postprocessor_point[dof][pt];
        }

        point_data.push_back(temp);
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F, typename T>
T Element<dimension, MasterType, ShapeType, DataType>::ComputeResidualL2(const F& f, const std::vector<T>& u) {
    std::pair<std::vector<double>, std::vector<Point<2>>> rule = this->master->integration.GetRule(20);
    // At this point we use maximum possible p for Dunavant integration

    Array2D<double> Phi = this->master->basis.GetPhi(this->master->p, rule.second);

    std::vector<T> u_gp(rule.first.size());
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < this->data.get_ndof(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += Phi[dof][gp] * u[dof];
        }
    }

    std::vector<Point<2>> gp_global = this->shape.LocalToGlobalCoordinates(rule.second);

    std::vector<T> f_gp(rule.first.size());

    for (uint gp = 0; gp < f_gp.size(); gp++) {
        f_gp[gp] = f(gp_global[gp]);
    }

    std::vector<T> sq_diff(rule.first.size());

    for (uint gp = 0; gp < sq_diff.size(); gp++) {
        sq_diff[gp] = (f_gp[gp] - u_gp[gp]) * (f_gp[gp] - u_gp[gp]);
    }

    T L2;

    L2 = 0;

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
