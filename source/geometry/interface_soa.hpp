#ifndef GEOMETRY_INTERFACE_SOA_HPP
#define GEOMETRY_INTERFACE_SOA_HPP

namespace Geometry {

namespace detail {

template <typename InterfaceType>
struct InterfaceSoAHelper {
    constexpr static uint dimension = InterfaceType::dimension;

    //fixme: add phi_gp_in/ex to flip elements when needed
    // n_dof x ngp
    StatVector<DynMatrix<double>,3> phi_gp_bdry;

    StatVector<DynMatrix<double>,3> int_phi_fact;
    StatVector<DiagonalMatrix<double>,3> abs_J;

    StatVector<SparseMatrix<double>, 3> scatter_in;
    StatVector<SparseMatrix<double>, 3> scatter_ex;

    StatVector<SparseMatrix<double>, 3> gather_in;
    StatVector<SparseMatrix<double>, 3> gather_ex;

    template <typename ArrayType>
    decltype(auto) ComputeUgpBdry(const ArrayType& u, uint bdry_id) {
        return u * phi_gp_bdry[bdry_id];
    }

    template <typename ArrayType>
    decltype(auto) IntegratePhiBdry(const ArrayType& u, uint bdry_id) {
        return (abs_J[bdry_id] * u) * int_phi_fact[bdry_id];
    }

    InterfaceSoAHelper() = default;

    template <typename ElementContainer>
    InterfaceSoAHelper(uint p, ElementContainer& elt_container) {

        //Get **interface** integration rule
        typename InterfaceType::InterfaceIntegrationType integration;

        std::pair<DynVector<double>, std::vector<Point<dimension>>> integration_rule
            = integration.GetRule(2 * p + 1);

        uint ngp = integration_rule.second.size();
        auto& master = elt_container.GetMaster();
        for ( uint side = 0; side < 3; ++side ) {
            std::vector<Point<dimension + 1>> z_master =
                master.BoundaryToMasterCoordinates(side, integration_rule.second);

            this->phi_gp_bdry[side] = master.basis.GetPhi(p, z_master);

            this->int_phi_fact[side] = transpose(this->phi_gp_bdry[side]);
            for (uint dof = 0; dof < master.ndof; ++dof) {
                for (uint gp = 0; gp < ngp; ++gp) {
                    this->int_phi_fact[side](gp, dof) *= integration_rule.first[gp];
                }
            }
        }
    }

    void reserve(uint ninterfaces, uint nelements) {
        for ( uint side = 0; side < 3; ++side ) {
            scatter_in[side].resize(ninterfaces, nelements);
            scatter_in[side].reserve(nelements);

            scatter_ex[side].resize(ninterfaces, nelements);
            scatter_ex[side].reserve(nelements);

            gather_in[side].resize(nelements, ninterfaces);
            gather_in[side].reserve(nelements);

            gather_ex[side].resize(nelements, ninterfaces);
            gather_ex[side].reserve(nelements);

            abs_J[side] = DiagonalMatrix<double>(nelements,nelements);
            for ( uint i = 0; i < nelements; ++i ) {
                abs_J[side](i,i) = 0.;
            }
        }
    }

    void finalize_initialization() {
        /* fixme: shrinkToFit() causes memory leaks
        for ( uint side = 0; side < 3; ++side ) {
            scatter_in[side].shrinkToFit();
            scatter_ex[side].shrinkToFit();

            gather_in[side].shrinkToFit();
            gather_ex[side].shrinkToFit();
            }*/
    }
};

}

template <typename InterfaceType, typename InterfaceDataSoAType, typename ElementContainers>
class InterfaceSoA;

template <typename InterfaceType, typename InterfaceDataSoAType, typename... ElementContainer>
class InterfaceSoA<InterfaceType,
                   InterfaceDataSoAType,
                   std::tuple<ElementContainer...>> {
public:
    using AccessorType = InterfaceType;
    using specialization_t = typename InterfaceType::specialization_t;
    using ElementContainers = std::tuple<ElementContainer...>;

    constexpr static uint dimension = AccessorType::dimension;
    InterfaceDataSoAType data;
    StatVector<DynMatrix<double, SO::ColumnMajor>,dimension+1> surface_normal;

    InterfaceSoA() = default;
    InterfaceSoA(const uint p, ElementContainers& element_data_) : element_data(&element_data_) {
        uint helper_index = 0;
        Utilities::for_each_in_tuple( *element_data, [this, p, &helper_index](auto& elt_container) {
                this->helpers[helper_index++] = detail::InterfaceSoAHelper<InterfaceType>(p, elt_container);
            });
    }

    void reserve(uint ninterfaces) {
        assert(element_data);

        uint helper_index = 0u;
        Utilities::for_each_in_tuple( *element_data, [this, &helper_index, ninterfaces](auto& elt_container) {
                this->helpers[helper_index++].reserve(ninterfaces, elt_container.size());
            });

        data = InterfaceDataSoAType(ninterfaces,this->Getngp());

        surface_normal[GlobalCoord::x].resize(ninterfaces,this->Getngp());
        surface_normal[GlobalCoord::y].resize(ninterfaces,this->Getngp());

    }

    void finalize_initialization() {
        assert(element_data);

        uint helper_index = 0;
        Utilities::for_each_in_tuple( *element_data, [this, &helper_index] (auto& elt_container) {
                this->helpers[helper_index++].finalize_initialization();
            });
    }

    template <typename... Args>
    AccessorType at(uint intfc_index,
                    Args&&... args) {
        uint helper_index = 0u;

        AccessorType intfc(std::forward<Args>(args)...);

        assert(element_data);
        //GetLocalIndex returns a negative number if pointer is not found
        Utilities::for_each_in_tuple( *element_data, [this, intfc_index, &helper_index, &intfc] (auto& elt_container) {
                                          auto& helper = this->helpers[helper_index++];

                                          uint elt_indx_in = elt_container.GetLocalIndex(intfc.data_in);
                                          if ( elt_indx_in >= 0 ) {
                                              helper.scatter_in[intfc.bound_id_in].insert(intfc_index,
                                                                                          elt_indx_in, 1.);

                                              helper.gather_in[intfc.bound_id_in].insert(elt_indx_in,
                                                                                         intfc_index, 1.);

                                              helper.abs_J[intfc.bound_id_in](elt_indx_in, elt_indx_in) =
                                                  intfc.GetAbsJ();
                                          }

                                          uint elt_indx_ex = elt_container.GetLocalIndex(intfc.data_ex);
                                          if ( elt_indx_ex >= 0 ) {
                                              helper.scatter_ex[intfc.bound_id_ex].insert(intfc_index,
                                                                                          elt_indx_ex, 1.);

                                              helper.gather_ex[intfc.bound_id_ex].insert(elt_indx_ex,
                                                                                         intfc_index, 1.);

                                              helper.abs_J[intfc.bound_id_ex](elt_indx_ex, elt_indx_ex) =
                                                  intfc.GetAbsJ();
                                          }
            });

        return intfc;
    }

    void SetNormal(uint index, const std::array<DynRowVector<double>,dimension+1>& surface_normal_) {
        assert(surface_normal_[0].size() == this->Getngp());
        for ( uint dim = 0; dim < dimension+1; ++dim ) {
            for ( uint gp = 0; gp < this->Getngp(); ++gp ) {
                this->surface_normal[dim](index,gp) = surface_normal_[dim][gp];
            }
        }
    }

    void SetElementData(ElementContainers& element_data_) {
        this->element_data = &element_data_;
    }

    ElementContainers& GetElementData() {
        assert(element_data);
        return *element_data;
    }

    template <typename ArrayType>
    decltype(auto) ComputeUgpBdry(const ArrayType& u, uint bdry_id, const uint element_type_index) {
        return helpers[element_type_index].ComputeUgpBdry(u, bdry_id);
    }

    template <typename ArrayType>
    decltype(auto) IntegratePhiBdry(const ArrayType& u, uint bdry_id, const uint element_type_index) {
        return helpers[element_type_index].IntegratePhiBdry(u, bdry_id);
    }

    template <typename ArrayType>
    decltype(auto) GatherIn(uint element_type_index, uint side, const ArrayType& array) {
        return this->helpers[element_type_index].gather_in[side] * array;
    }

    template <typename ArrayType>
    decltype(auto) GatherEx(uint element_type_index, uint side, const ArrayType& array) {
        return reverse_columns(this->helpers[element_type_index].gather_ex[side] * array);
    }

    template <typename ArrayType>
    decltype(auto) ScatterIn(uint element_type_index, uint side, const ArrayType& array) {
        return this->helpers[element_type_index].scatter_in[side] * array;
    }

    template <typename ArrayType>
    decltype(auto) ScatterEx(uint element_type_index, uint side, const ArrayType& array) {
        return reverse_columns(this->helpers[element_type_index].scatter_ex[side] * array);
    }

    uint Getninterfaces() const {
        return rows(surface_normal[0]);
    }

    uint Getngp() const {
        return columns(helpers[0].phi_gp_bdry[0]);
    }
private:
    constexpr static size_t n_element_types = sizeof...(ElementContainer);

    //n_interfaces x n_dimensions
    std::tuple<ElementContainer...> * element_data = nullptr;
    std::array<detail::InterfaceSoAHelper<AccessorType>, n_element_types> helpers;
};
}

#endif