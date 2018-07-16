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

    DynVector<uint> node_ID;
    DynVector<uint> neighbor_ID;
    DynVector<uchar> boundary_type;

    DynVector<Point<dimension>> gp_global_coordinates;

    bool const_J;

    /* psi_gp stored in shape */   // nodal basis, i.e. shape functions
    /* chi_gp stored in master */  // linear basis
    /* phi_gp stroed in master */  // modal basis

    /* dpsi_gp stored in shape */                      // nodal basis, i.e. shape functions
    StatVector<DynMatrix<double>, dimension> dchi_gp;  // linear basis
    StatVector<DynMatrix<double>, dimension> dphi_gp;  // modal basis

    DynVector<double> int_fact;
    DynMatrix<double> int_phi_fact;
    DynMatrix<double> int_phi_phi_fact;
    StatVector<DynMatrix<double>, dimension> int_dphi_fact;
    StatVector<DynMatrix<double>, dimension> int_phi_dphi_fact;

    std::pair<bool, DynMatrix<double>> m_inv;

  public:
    Element() = default;
    Element(const uint ID,
            MasterType& master,
            DynVector<Point<3>>&& nodal_coordinates,
            DynVector<uint>&& node_ID,
            DynVector<uint>&& neighbor_ID,
            DynVector<uchar>&& boundary_type);

    uint GetID() { return this->ID; }
    MasterType& GetMaster() { return *this->master; }
    ShapeType& GetShape() { return this->shape; }
    DynVector<uint>& GetNodeID() { return this->node_ID; }
    DynVector<uchar>& GetBoundaryType() { return this->boundary_type; }

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
                                                             DynVector<Point<3>>&& nodal_coordinates,
                                                             DynVector<uint>&& node_ID,
                                                             DynVector<uint>&& neighbor_ID,
                                                             DynVector<uchar>&& boundary_type)
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
    DynVector<double> det_J = this->shape.GetJdet(this->master->integration_rule.second);
    DynVector<StatMatrix<double, dimension, dimension>> J_inv =
        this->shape.GetJinv(this->master->integration_rule.second);

    this->const_J = (det_J.size() == 1);

    // Compute factors to expand nodal values and derivatives of nodal values
    this->shape.psi_gp  = this->shape.GetPsi(this->master->integration_rule.second);
    this->shape.dpsi_gp = this->shape.GetDPsi(this->master->integration_rule.second);

    if (const_J) {  // constant Jacobian
        // DIFFERENTIATION FACTORS
        this->dchi_gp = this->master->dchi_gp;
        for (uint dir = 0; dir < dimension; dir++) {
            for (uint gp = 0; gp < this->master->ngp; gp++) {
                for (uint dof = 0; dof < this->master->nvrtx; dof++) {
                    double dchi = 0;
                    for (uint z = 0; z < dimension; z++) {
                        dchi += this->master->dchi_gp[z](gp, dof) * J_inv[0](z, dir);
                    }
                    this->dchi_gp[dir](gp, dof) = dchi;
                }
            }
        }

        this->dphi_gp = this->master->dphi_gp;
        for (uint dir = 0; dir < dimension; dir++) {
            for (uint gp = 0; gp < this->master->ngp; gp++) {
                for (uint dof = 0; dof < this->master->ndof; dof++) {
                    double dphi = 0.0;
                    for (uint z = 0; z < dimension; z++) {
                        dphi += this->master->dphi_gp[z](gp, dof) * J_inv[0](z, dir);
                    }
                    this->dphi_gp[dir](gp, dof) = dphi;
                }
            }
        }

        // INTEGRATION OVER ELEMENT FACTORS
        this->int_fact = this->master->integration_rule.first * std::abs(det_J[0]);

        this->int_phi_fact = this->master->int_phi_fact * std::abs(det_J[0]);

        this->int_phi_phi_fact.resize(std::pow(this->master->ndof, 2), this->master->ngp);
        for (uint dof_i = 0; dof_i < this->master->ndof; dof_i++) {
            for (uint dof_j = 0; dof_j < this->master->ndof; dof_j++) {
                uint lookup = this->master->ndof * dof_i + dof_j;
                for (uint gp = 0; gp < this->master->ngp; gp++) {
                    this->int_phi_phi_fact(lookup, gp) =
                        this->master->phi_gp(gp, dof_i) * this->int_phi_fact(dof_j, gp);
                }
            }
        }

        this->int_dphi_fact = this->master->int_dphi_fact;
        for (uint dir = 0; dir < dimension; dir++) {
            for (uint dof = 0; dof < this->master->ndof; dof++) {
                for (uint gp = 0; gp < this->master->ngp; gp++) {
                    double int_dphi = 0;
                    for (uint z = 0; z < dimension; z++) {
                        int_dphi += this->master->int_dphi_fact[z](dof, gp) * J_inv[0](z, dir);
                    }
                    int_dphi *= std::abs(det_J[0]);
                    this->int_dphi_fact[dir](dof, gp) = int_dphi;
                }
            }
        }

        for (uint dir = 0; dir < dimension; dir++) {
            this->int_phi_dphi_fact[dir].resize(std::pow(this->master->ndof, 2), this->master->ngp);
            for (uint dof_i = 0; dof_i < this->master->ndof; dof_i++) {
                for (uint dof_j = 0; dof_j < this->master->ndof; dof_j++) {
                    uint lookup = this->master->ndof * dof_i + dof_j;
                    for (uint gp = 0; gp < this->master->ngp; gp++) {
                        this->int_phi_dphi_fact[dir](lookup, gp) =
                            this->master->phi_gp(gp, dof_i) * this->int_dphi_fact[dir](dof_j, gp);
                    }
                }
            }
        }

        // MASS MATRIX
        this->m_inv = this->master->m_inv;
        this->m_inv.second /= std::abs(det_J[0]);
    } else {
        // Placeholder for nonconstant Jacobian
    }

    this->data.set_nnode(this->shape.nodal_coordinates.size());
    this->data.set_nvrtx(this->master->nvrtx);
    this->data.set_nbound(this->master->nbound);
    this->data.set_ndof(this->master->ndof);
    this->data.set_ngp_internal(this->master->ngp);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
void Element<dimension, MasterType, ShapeType, DataType>::CreateRawBoundaries(
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundary<dimension - 1, DataType>>>& raw_boundaries) {
    // *** //
    Basis::Basis<dimension>* my_basis    = (Basis::Basis<dimension>*)(&this->master->basis);
    Master::Master<dimension>* my_master = (Master::Master<dimension>*)(this->master);
    Shape::Shape<dimension>* my_shape    = (Shape::Shape<dimension>*)(&this->shape);

    for (uint bound_id = 0; bound_id < this->boundary_type.size(); bound_id++) {
        DynVector<uint> bound_node_ID = this->shape.GetBoundaryNodeID(bound_id, this->node_ID);

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

    for (uint dof = 0; dof < this->master->ndof; dof++) {
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

    for (uint dof = 0; dof < this->master->ndof; dof++) {
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
            u_gp[gp] += u[dof] * this->master->phi_gp(gp, dof);
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
            du_gp[gp] += u[dof] * this->dphi_gp[dir](gp, dof);
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ComputeDUgp(
    const std::vector<T>& u,
    std::vector<StatVector<T, dimension>>& du_gp) {
    std::fill(du_gp.begin(), du_gp.end(), 0.0);

    for (uint dof = 0; dof < u.size(); dof++) {
        for (uint gp = 0; gp < du_gp.size(); gp++) {
            for (uint dir = 0; dir < dimension; dir++) {
                du_gp[gp][dir] += u[dof] * this->dphi_gp[dir](gp, dof);
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
            u_lin_gp[gp] += u_lin[dof] * this->master->chi_gp(gp, dof);
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
            du_lin_gp[gp] += u_lin[dof] * this->dchi_gp[dir](gp, dof);
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
            u_nodal_gp[gp] += u_nodal[dof] * this->shape.psi_gp(gp, dof);
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
            du_nodal_gp[gp] += u_nodal[dof] * this->shape.dpsi_gp[dir](gp, dof);
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
                du_nodal_gp[gp][dir] += u_nodal[dof] * this->shape.dpsi_gp[dir](gp, dof);
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
        integral += u_gp[gp] * this->int_phi_fact(dof, gp);
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

    uint lookup = this->master->ndof * dof_i + dof_j;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_phi_fact(lookup, gp);
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
        integral += u_gp[gp] * this->int_dphi_fact[dir](dof, gp);
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

    uint lookup = this->master->ndof * dof_i + dof_j;

    for (uint gp = 0; gp < u_gp.size(); gp++) {
        integral += u_gp[gp] * this->int_phi_dphi_fact[dir_j](lookup, gp);
    }

    return integral;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::ApplyMinv(const std::vector<T>& rhs,
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

    for (uint cell = 0; cell < N_DIV * N_DIV; cell++) {
        temp = 0.0;

        for (uint dof = 0; dof < u.size(); dof++) {
            temp += u[dof] * this->master->phi_postprocessor_cell(cell, dof);
        }

        cell_data.push_back(temp);
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename T>
inline void Element<dimension, MasterType, ShapeType, DataType>::WritePointDataVTK(const std::vector<T>& u,
                                                                                   std::vector<T>& point_data) {
    T temp;

    for (uint pt = 0; pt < (N_DIV + 1) * (N_DIV + 2) / 2; pt++) {
        temp = 0.0;

        for (uint dof = 0; dof < u.size(); dof++) {
            temp += u[dof] * this->master->phi_postprocessor_point(pt, dof);
        }

        point_data.push_back(temp);
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F, typename T>
T Element<dimension, MasterType, ShapeType, DataType>::ComputeResidualL2(const F& f, const std::vector<T>& u) {
    std::pair<DynVector<double>, DynVector<Point<2>>> rule = this->master->integration.GetRule(20);
    // At this point we use maximum possible p for Dunavant integration

    DynMatrix<double> Phi = this->master->basis.GetPhi(this->master->p, rule.second);

    std::vector<T> u_gp(rule.first.size());
    std::fill(u_gp.begin(), u_gp.end(), 0.0);

    for (uint dof = 0; dof < this->data.get_ndof(); dof++) {
        for (uint gp = 0; gp < u_gp.size(); gp++) {
            u_gp[gp] += Phi(gp, dof) * u[dof];
        }
    }

    DynVector<Point<2>> gp_global = this->shape.LocalToGlobalCoordinates(rule.second);

    DynVector<T> f_gp(rule.first.size());

    for (uint gp = 0; gp < f_gp.size(); gp++) {
        f_gp[gp] = f(gp_global[gp]);
    }

    DynVector<T> sq_diff(rule.first.size());

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
