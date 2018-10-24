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

    bool const_J;

    /* psi_gp stored in shape */   // nodal basis, i.e. shape functions
    /* chi_gp stored in master */  // linear basis
    /* phi_gp stroed in master */  // modal basis

    /* dpsi_gp stored in shape */                      // nodal basis, i.e. shape functions
    std::array<DynMatrix<double>, dimension> dchi_gp;  // linear basis
    std::array<DynMatrix<double>, dimension> dphi_gp;  // modal basis

    DynVector<double> int_fact;
    DynMatrix<double> int_phi_fact;
    DynMatrix<double> int_phi_phi_fact;
    std::array<DynMatrix<double>, dimension> int_dphi_fact;
    std::array<DynMatrix<double>, dimension> int_phi_dphi_fact;

    DynMatrix<double> m_inv;

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

    template <typename F>
    DynMatrix<double> L2ProjectionF(const F& f);
    template <typename InputArrayType>
    decltype(auto) L2ProjectionNode(const InputArrayType& nodal_values);

    template <typename InputArrayType>
    DynMatrix<double> ProjectBasisToLinear(const InputArrayType& u);
    template <typename InputArrayType>
    DynMatrix<double> ProjectLinearToBasis(const uint ndof, const InputArrayType& u_lin);

    template <typename F>
    DynMatrix<double> ComputeFgp(const F& f);
    template <typename InputArrayType>
    decltype(auto) ComputeUgp(const InputArrayType& u);
    template <typename InputArrayType>
    decltype(auto) ComputeDUgp(const uint dir, const InputArrayType& u);

    template <typename InputArrayType>
    decltype(auto) ComputeLinearUgp(const InputArrayType& u_lin);
    template <typename InputArrayType>
    decltype(auto) ComputeLinearDUgp(const uint dir, const InputArrayType& u_lin);
    template <typename InputArrayType>
    decltype(auto) ComputeLinearUbaryctr(const InputArrayType& u_lin);
    template <typename InputArrayType>
    decltype(auto) ComputeLinearUmidpts(const InputArrayType& u_lin);
    template <typename InputArrayType>
    decltype(auto) ComputeLinearUvrtx(const InputArrayType& u_lin);

    template <typename InputArrayType>
    decltype(auto) ComputeNodalUgp(const InputArrayType& u_nodal);
    template <typename InputArrayType>
    decltype(auto) ComputeNodalDUgp(const uint dir, const InputArrayType& u_nodal);

    template <typename InputArrayType>
    decltype(auto) Integration(const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhi(const uint dof, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhi(const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiPhi(const uint dof_i, const uint dof_j, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationDPhi(const uint dir, const uint dof, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationDPhi(const uint dir, const InputArrayType& u_gp);
    template <typename InputArrayType>
    decltype(auto) IntegrationPhiDPhi(const uint dof_i, const uint dir_j, const uint dof_j, const InputArrayType& u_gp);

    template <typename InputArrayType>
    decltype(auto) ApplyMinv(const InputArrayType& rhs);

    void InitializeVTK(std::vector<Point<3>>& points, Array2D<uint>& cells);
    template <typename InputArrayType, typename OutputArrayType>
    void WriteCellDataVTK(const InputArrayType& u, AlignedVector<OutputArrayType>& cell_data);
    template <typename InputArrayType, typename OutputArrayType>
    void WritePointDataVTK(const InputArrayType& u, AlignedVector<OutputArrayType>& point_data);

    template <typename F, typename InputArrayType>
    double ComputeResidualL2(const F& f, const InputArrayType& u);

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
    DynVector<double> det_J = this->shape.GetJdet(this->master->integration_rule.second);
    AlignedVector<StatMatrix<double, dimension, dimension>> J_inv =
        this->shape.GetJinv(this->master->integration_rule.second);

    this->const_J = (det_J.size() == 1);

    // Compute factors to expand nodal values and derivatives of nodal values
    this->shape.psi_gp  = this->shape.GetPsi(this->master->integration_rule.second);
    this->shape.dpsi_gp = this->shape.GetDPsi(this->master->integration_rule.second);

    if (const_J) {  // constant Jacobian
        // DIFFERENTIATION FACTORS
        this->dchi_gp = this->master->dchi_gp;
        for (uint dir = 0; dir < dimension; ++dir) {
            for (uint gp = 0; gp < this->master->ngp; ++gp) {
                for (uint dof = 0; dof < this->master->nvrtx; ++dof) {
                    double dchi = 0;
                    for (uint z = 0; z < dimension; ++z) {
                        dchi += this->master->dchi_gp[z](dof, gp) * J_inv[0](z, dir);
                    }
                    this->dchi_gp[dir](dof, gp) = dchi;
                }
            }
        }

        this->dphi_gp = this->master->dphi_gp;
        for (uint dir = 0; dir < dimension; ++dir) {
            for (uint gp = 0; gp < this->master->ngp; ++gp) {
                for (uint dof = 0; dof < this->master->ndof; ++dof) {
                    double dphi = 0.0;
                    for (uint z = 0; z < dimension; ++z) {
                        dphi += this->master->dphi_gp[z](dof, gp) * J_inv[0](z, dir);
                    }
                    this->dphi_gp[dir](dof, gp) = dphi;
                }
            }
        }

        // INTEGRATION OVER ELEMENT FACTORS
        this->int_fact = this->master->integration_rule.first * std::abs(det_J[0]);

        this->int_phi_fact = this->master->int_phi_fact * std::abs(det_J[0]);

        this->int_phi_phi_fact.resize(this->master->ngp, std::pow(this->master->ndof, 2));
        for (uint dof_i = 0; dof_i < this->master->ndof; ++dof_i) {
            for (uint dof_j = 0; dof_j < this->master->ndof; ++dof_j) {
                uint lookup = this->master->ndof * dof_i + dof_j;
                for (uint gp = 0; gp < this->master->ngp; ++gp) {
                    this->int_phi_phi_fact(gp, lookup) =
                        this->master->phi_gp(dof_i, gp) * this->int_phi_fact(gp, dof_j);
                }
            }
        }

        this->int_dphi_fact = this->master->int_dphi_fact;
        for (uint dir = 0; dir < dimension; ++dir) {
            for (uint dof = 0; dof < this->master->ndof; ++dof) {
                for (uint gp = 0; gp < this->master->ngp; ++gp) {
                    double int_dphi = 0;
                    for (uint z = 0; z < dimension; ++z) {
                        int_dphi += this->master->int_dphi_fact[z](gp, dof) * J_inv[0](z, dir);
                    }
                    int_dphi *= std::abs(det_J[0]);
                    this->int_dphi_fact[dir](gp, dof) = int_dphi;
                }
            }
        }

        for (uint dir = 0; dir < dimension; ++dir) {
            this->int_phi_dphi_fact[dir].resize(this->master->ngp, std::pow(this->master->ndof, 2));
            for (uint dof_i = 0; dof_i < this->master->ndof; ++dof_i) {
                for (uint dof_j = 0; dof_j < this->master->ndof; ++dof_j) {
                    uint lookup = this->master->ndof * dof_i + dof_j;
                    for (uint gp = 0; gp < this->master->ngp; ++gp) {
                        this->int_phi_dphi_fact[dir](gp, lookup) =
                            this->master->phi_gp(dof_i, gp) * this->int_dphi_fact[dir](gp, dof_j);
                    }
                }
            }
        }

        // MASS MATRIX
        this->m_inv = this->master->m_inv / std::abs(det_J[0]);
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

    for (uint bound_id = 0; bound_id < this->boundary_type.size(); ++bound_id) {
        std::vector<uint> bound_node_ID = this->shape.GetBoundaryNodeID(bound_id, this->node_ID);

        if (is_internal(this->boundary_type[bound_id])) {
            raw_boundaries[this->boundary_type[bound_id]].emplace(
                std::pair<uint, uint>{this->ID, this->neighbor_ID[bound_id]},
                RawBoundary<dimension - 1, DataType>(
                    this->master->p, bound_id, std::move(bound_node_ID), this->data, *my_basis, *my_master, *my_shape));
        } else {
            raw_boundaries[this->boundary_type[bound_id]].emplace(
                std::pair<uint, uint>{this->ID, bound_id},
                RawBoundary<dimension - 1, DataType>(
                    this->master->p, bound_id, std::move(bound_node_ID), this->data, *my_basis, *my_master, *my_shape));
        }
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F>
inline DynMatrix<double> Element<dimension, MasterType, ShapeType, DataType>::L2ProjectionF(const F& f) {
    // projection(q, dof) = f_values(q, gp) * int_phi_fact(gp, dof) * m_inv(dof, dof)
    DynMatrix<double> projection = this->ComputeFgp(f) * this->int_phi_fact * this->m_inv;

    return projection;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::L2ProjectionNode(
    const InputArrayType& nodal_values) {
    // projection(q, dof) = nodal_values(q, node) * psi_gp(node, gp) * int_phi_fact(gp, dof) * m_inv(dof, dof)
    return nodal_values * this->shape.psi_gp * this->int_phi_fact * this->m_inv;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline DynMatrix<double> Element<dimension, MasterType, ShapeType, DataType>::ProjectBasisToLinear(
    const InputArrayType& u) {
    if (const_J) {
        return this->master->basis.ProjectBasisToLinear(u);
    } else {
        return DynMatrix<double>();
        // Placeholder for nonconstant Jacobian
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline DynMatrix<double> Element<dimension, MasterType, ShapeType, DataType>::ProjectLinearToBasis(
    const uint ndof,
    const InputArrayType& u_lin) {
    if (const_J) {
        return this->master->basis.ProjectLinearToBasis(ndof, u_lin);
    } else {
        return DynMatrix<double>();
        // Placeholder for nonconstant Jacobian
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F>
inline DynMatrix<double> Element<dimension, MasterType, ShapeType, DataType>::ComputeFgp(const F& f) {
    uint nvar = f(this->gp_global_coordinates[0]).size();
    uint ngp  = this->gp_global_coordinates.size();

    DynMatrix<double> f_vals(nvar, ngp);

    for (uint gp = 0; gp < this->gp_global_coordinates.size(); ++gp) {
        column(f_vals, gp) = f(this->gp_global_coordinates[gp]);
    }

    return f_vals;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ComputeUgp(const InputArrayType& u) {
    // u_gp(q, gp) = u(q, dof) * phi_gp(dof, gp)
    return u * this->master->phi_gp;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ComputeDUgp(const uint dir,
                                                                                       const InputArrayType& u) {
    // du_gp(q, gp) = u(q, dof) * dphi_gp[dir](dof, gp)
    return u * this->dphi_gp[dir];
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUgp(
    const InputArrayType& u_lin) {
    // u_lin_gp(q, gp) = u_lin(q, dof) * chi_gp(dof, gp)
    return u_lin * this->master->chi_gp;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearDUgp(
    const uint dir,
    const InputArrayType& u_lin) {
    // du_lin_gp(q, gp) = du(q, dof) * chi_gp[dir](dof, gp)
    return u_lin * this->dchi_gp[dir];
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUbaryctr(
    const InputArrayType& u_lin) {
    return this->master->ComputeLinearUbaryctr(u_lin);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUmidpts(
    const InputArrayType& u_lin) {
    return this->master->ComputeLinearUmidpts(u_lin);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ComputeLinearUvrtx(
    const InputArrayType& u_lin) {
    return this->master->ComputeLinearUvrtx(u_lin);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ComputeNodalUgp(
    const InputArrayType& u_nodal) {
    // u_nodal_gp(q, gp) = u_nodal(q, dof) * psi_gp(dof, gp)
    return u_nodal * this->shape.psi_gp;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ComputeNodalDUgp(
    const uint dir,
    const InputArrayType& u_nodal) {
    // du_nodal_gp(q, gp) = u_nodal(q, dof) * dpsi_gp[dir](dof, gp)
    return u_nodal * this->shape.dpsi_gp[dir];
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::Integration(const InputArrayType& u_gp) {
    // integral[q] = u_gp(q, gp) * this->int_fact[gp]
    return u_gp * this->int_fact;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::IntegrationPhi(const uint dof,
                                                                                          const InputArrayType& u_gp) {
    // integral[q] = u_gp(q, gp) * this->int_phi_fact(gp, dof)
    return u_gp * column(this->int_phi_fact, dof);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::IntegrationPhi(const InputArrayType& u_gp) {
    // integral(q, dof) = u_gp(q, gp) * this->int_phi_fact(gp, dof)
    return u_gp * this->int_phi_fact;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::IntegrationPhiPhi(
    const uint dof_i,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] = u_gp(q, gp) * this->int_phi_phi_fact(gp, lookup)
    uint lookup = this->master->ndof * dof_i + dof_j;

    return u_gp * column(this->int_phi_phi_fact, lookup);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::IntegrationDPhi(const uint dir,
                                                                                           const uint dof,
                                                                                           const InputArrayType& u_gp) {
    // integral[q] =  u_gp(q, gp) * this->int_dphi_fact[dir](gp. dof)
    return u_gp * column(this->int_dphi_fact[dir], dof);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::IntegrationDPhi(const uint dir,
                                                                                           const InputArrayType& u_gp) {
    // integral(q, dof) =  u_gp(q, gp) * this->int_dphi_fact[dir](gp. dof)
    return u_gp * this->int_dphi_fact[dir];
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::IntegrationPhiDPhi(
    const uint dof_i,
    const uint dir_j,
    const uint dof_j,
    const InputArrayType& u_gp) {
    // integral[q] = u_gp(q, gp) * this->int_phi_dphi_fact[dir_j](lookup, gp)
    uint lookup = this->master->ndof * dof_i + dof_j;

    return u_gp * column(this->int_phi_dphi_fact[dir_j], lookup);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType>
inline decltype(auto) Element<dimension, MasterType, ShapeType, DataType>::ApplyMinv(const InputArrayType& rhs) {
    // solution(q, dof) = rhs(q, dof) * this->m_inv(dof, dof)
    return rhs * this->m_inv;
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
void Element<dimension, MasterType, ShapeType, DataType>::InitializeVTK(std::vector<Point<3>>& points,
                                                                        Array2D<uint>& cells) {
    this->shape.GetVTK(points, cells);
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType, typename OutputArrayType>
inline void Element<dimension, MasterType, ShapeType, DataType>::WriteCellDataVTK(
    const InputArrayType& u,
    AlignedVector<OutputArrayType>& cell_data) {
    // cell_data[q] = u(q, dof) * phi_postprocessor_cell(dof, cell)
    OutputArrayType temp;

    for (uint cell = 0; cell < columns(this->master->phi_postprocessor_cell); ++cell) {
        temp = u * column(this->master->phi_postprocessor_cell, cell);

        cell_data.emplace_back(std::move(temp));
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename InputArrayType, typename OutputArrayType>
inline void Element<dimension, MasterType, ShapeType, DataType>::WritePointDataVTK(
    const InputArrayType& u,
    AlignedVector<OutputArrayType>& point_data) {
    // point_data[q] = u(q, dof) * phi_postprocessor_point(pt, dof)
    OutputArrayType temp;

    for (uint pt = 0; pt < columns(this->master->phi_postprocessor_point); ++pt) {
        temp = u * column(this->master->phi_postprocessor_point, pt);

        point_data.emplace_back(std::move(temp));
    }
}

template <uint dimension, typename MasterType, typename ShapeType, typename DataType>
template <typename F, typename InputArrayType>
double Element<dimension, MasterType, ShapeType, DataType>::ComputeResidualL2(const F& f, const InputArrayType& u) {
    // At this point we use maximum possible p for Dunavant integration
    std::pair<DynVector<double>, std::vector<Point<2>>> rule = this->master->integration.GetRule(20);

    // get u_gp
    DynMatrix<double> phi_gp = this->master->basis.GetPhi(this->master->p, rule.second);
    DynMatrix<double> u_gp   = u * phi_gp;

    // get true_gp
    std::vector<Point<dimension>> gp_global = this->shape.LocalToGlobalCoordinates(rule.second);

    uint nvar = f(*(gp_global.begin())).size();
    uint ngp  = gp_global.size();

    DynMatrix<double> true_gp(nvar, ngp);

    for (uint gp = 0; gp < ngp; ++gp) {
        column(true_gp, gp) = f(gp_global[gp]);
    }

    // find square difference betwee u_gp and true_gp
    DynMatrix<double> diff = true_gp - u_gp;
    DynMatrix<double> sq_diff(nvar, ngp);

    for (uint var = 0; var < nvar; ++var) {
        row(sq_diff, var) = vec_cw_mult(row(diff, var), row(diff, var));
    }

    DynVector<double> L2;

    // integrate over element
    if (const_J) {
        L2 = sq_diff * rule.first * std::abs(this->shape.GetJdet(rule.second)[0]);
    } else {
        // Placeholder for nonconstant Jacobian
    }

    return L2[0];
}
}

#endif
