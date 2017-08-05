#ifndef MESH_HPP
#define MESH_HPP

#include "../utilities/heterogeneous_containers.hpp"

#include "mesh_utilities.hpp"

namespace Geometry {
// Since elements types already come in a tuple. We can use specialization
// to get easy access to the parameter packs for the element and edge types.
template <typename ElementTypeTuple, typename InteriorEdgeTypeTuple, typename BoundaryEdgeTypeTuple>
class Mesh;

template <typename... Elements, typename... Interfaces, typename... Boundaries>
class Mesh<std::tuple<Elements...>, std::tuple<Interfaces...>, std::tuple<Boundaries...>> {
  private:
    using MasterElementTypes = typename make_master_type<std::tuple<Elements...>>::type;
    using ElementContainer = Utilities::HeterogeneousMap<Elements...>;
    using InterfaceContainer = Utilities::HeterogeneousVector<Interfaces...>;
    using BoundaryContainer = Utilities::HeterogeneousVector<Boundaries...>;

    MasterElementTypes masters;
    ElementContainer elements;
    InterfaceContainer interfaces;
    BoundaryContainer boundaries;

  public:
    Mesh(uint p) : masters(master_maker<MasterElementTypes>::construct_masters(p)) {}

    uint GetNumberElements() { return this->elements.size(); }
    uint GetNumberInterfaces() { return this->interfaces.size(); }
    uint GetNumberBoundaries() { return this->boundaries.size(); }

    template <typename T, typename... Args>
    void CreateElement(uint n, Args&&... args) {
        using MasterType = typename T::master_element_type;

        MasterType& master_elt = std::get<Utilities::index<MasterType, MasterElementTypes>::value>(masters);

        this->elements.template emplace<T>(n, create<T>(n, master_elt, std::forward<Args>(args)...));
    }

    template <typename T, typename... Args>
    void CreateInterface(Args&&... args) {
        this->interfaces.template emplace_back<T>(create<T>(std::forward<Args>(args)...));
    }

    template <typename T, typename... Args>
    void CreateBoundary(Args&&... args) {
        this->boundaries.template emplace_back<T>(create<T>(std::forward<Args>(args)...));
    }

    template <typename F>
    void CallForEachElement(const F& f) {
        Utilities::for_each_in_tuple(elements.data, [&f](auto& element_map) {
            std::for_each(element_map.begin(), element_map.end(), [&f](auto& pair) { f(pair.second); });
        });
    }

    template <typename F>
    void CallForEachInterface(F f) {
        Utilities::for_each_in_tuple(interfaces.data, [f](auto& interface_vector) {
            std::for_each(interface_vector.begin(), interface_vector.end(), f);
        });
    }

    template <typename F>
    void CallForEachBoundary(F f) {
        Utilities::for_each_in_tuple(boundaries.data, [f](auto& boundary_vector) {
            std::for_each(boundary_vector.begin(), boundary_vector.end(), f);
        });
    }
};
}

#endif