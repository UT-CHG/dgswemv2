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

    //n_interfaces x n_dimensions
    DynMatrix<double> surface_normal;

    StatVector<SparseMatrix<double>, 3> scatter_in;
    StatVector<SparseMatrix<double>, 3> scatter_ex;

    StatVector<SparseMatrix<double>, 3> gather_in;
    StatVector<SparseMatrix<double>, 3> gather_ex;

    template <typename ArrayType>
    decltype(auto) ComputeUgpBdry(const ArrayType& u, uint bdry_id) {
        return u * phi_gp_bdry[bdry_id];
    }

    InterfaceSoAHelper() = default;

    template <typename ElementContainer>
    InterfaceSoAHelper(uint p, const ElementContainer& elt_container) {

        { //assemble phi_gp

            //Get **interface** integration rule
            typename InterfaceType::InterfaceIntegrationType integration;

            std::pair<DynVector<double>, std::vector<Point<dimension>>> integration_rule
                = integration.GetRule(2 * p + 1);

            const auto& master = elt_container.GetMaster();
            for ( uint side = 0; side < 3; ++side ) {
                std::vector<Point<dimension + 1>> z_master =
                    master.BoundaryToMasterCoordinates(side, integration_rule.second);

                this->phi_gp_bdry[side] = master.basis.GetPhi(p, z_master);
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

        }
    }

    void finalize_initialization() {
        for ( uint side = 0; side < 3; ++side ) {
            scatter_in[side].shrinkToFit();
            scatter_ex[side].shrinkToFit();

            gather_in[side].shrinkToFit();
            gather_ex[side].shrinkToFit();
        }
    }
};

}

template <typename InterfaceType, typename ElementContainers>
class InterfaceSoA;

template <typename InterfaceType, typename... ElementContainer>
class InterfaceSoA<InterfaceType,
                   std::tuple<ElementContainer...>> {
public:
    using AccessorType = InterfaceType;
    using ElementContainers = std::tuple<ElementContainer...>;
    constexpr static uint dimension = AccessorType::dimension;

    InterfaceSoA() = default;
    InterfaceSoA(ElementContainers& element_data_) : element_data(&element_data_) {}

    void reserve(uint ninterfaces) {
        assert(element_data);

        uint helper_index = 0u;
        Utilities::for_each_in_tuple( *element_data, [this, &helper_index, ninterfaces](auto& elt_container) {
                this->helpers[helper_index++].reserve(ninterfaces, elt_container.size());
            });
    }

    void finalize_initialization() {
        assert(element_data);

        uint helper_index = 0;
        Utilities::for_each_in_tuple( *element_data, [this, &helper_index] (auto& elt_container) {
                this->helpers[helper_index++].finalize_initialization();
            });
    }

    template <typename RawBoundary, typename... Args>
    AccessorType at(uint intfc_index,
                    RawBoundary&& raw_boundary_in,
                    RawBoundary&& raw_boundary_ex,
                    Args&&... args) {
        uint helper_index = 0u;

        assert(element_data);
        //GetLocalIndex returns a negative number if pointer is not found
        Utilities::for_each_in_tuple( *element_data, [this, intfc_index, &helper_index,
                                                      &raw_boundary_in, &raw_boundary_ex] (auto& elt_container) {
                                          auto& helper = this->helpers[helper_index++];

                                          ptrdiff_t elt_indx_in = elt_container.GetLocalIndex(raw_boundary_in.data);
                                          if ( elt_indx_in >= 0 ) {
                                              helper.scatter_in[raw_boundary_in.bound_id].insert(intfc_index,
                                                                                                 elt_indx_in, 1.);

                                              helper.gather_in[raw_boundary_in.bound_id].insert(elt_indx_in,
                                                                                                intfc_index, 1.);
                                          }

                                          ptrdiff_t elt_indx_ex = elt_container.GetLocalIndex(raw_boundary_ex.data);
                                          if ( elt_indx_ex >= 0 ) {
                                              helper.scatter_ex[raw_boundary_ex.bound_id].insert(intfc_index,
                                                                                                 elt_indx_ex, 1.);

                                              helper.gather_ex[raw_boundary_ex.bound_id].insert(elt_indx_ex,
                                                                                                intfc_index, 1.);
                                          }
            });

        return AccessorType(std::move(raw_boundary_in), std::move(raw_boundary_ex), std::forward<Args>(args)...);
    }

    ElementContainers& GetElementData(uint var) {
        assert(element_data);
        return *element_data;
    }

    template <typename ElementContainerType, typename ArrayType>
    decltype(auto) ComputeUgpBdry(const ArrayType& u, uint bdry_id) {
        constexpr int indx = Utilities::index<ElementContainerType, ElementContainers>::value;

        return helpers[indx].ComputeUgpBdry(u, bdry_id);
    }

private:
    constexpr static size_t n_element_types = sizeof...(ElementContainer);

    std::tuple<ElementContainer...> * element_data = nullptr;
    std::array<detail::InterfaceSoAHelper<AccessorType>, n_element_types> helpers;
};
}

#endif