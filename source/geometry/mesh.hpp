#ifndef MESH_HPP
#define MESH_HPP

#include "general_definitions.hpp"

#include "utilities//heterogeneous_containers.hpp"
#include "mesh_utilities.hpp"

namespace Geometry {
// Since elements types already come in a tuple. We can use specialization
// to get easy access to the parameter packs for the element and edge types.
template <typename ElementTypeTuple,
          typename InterfaceTypeTuple,
          typename BoundaryTypeTuple,
          typename DistributedBoundaryTypeTuple>
class Mesh;

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
class Mesh<std::tuple<Elements...>,
           std::tuple<Interfaces...>,
           std::tuple<Boundaries...>,
           std::tuple<DistributedBoundaries...>> {
  public:
    using MasterElementTypes           = typename make_master_type<std::tuple<Elements...>>::type;
    using ElementContainer             = Utilities::HeterogeneousMap<Elements...>;
    using InterfaceContainer           = Utilities::HeterogeneousVector<Interfaces...>;
    using BoundaryContainer            = Utilities::HeterogeneousVector<Boundaries...>;
    using DistributedBoundaryContainer = Utilities::HeterogeneousVector<DistributedBoundaries...>;

  private:
    uint p;

    MasterElementTypes masters;
    ElementContainer elements;
    InterfaceContainer interfaces;
    BoundaryContainer boundaries;
    DistributedBoundaryContainer distributed_boundaries;

    std::string mesh_name;

  public:
    Mesh() = default;
    Mesh(const uint p) : p(p), masters(master_maker<MasterElementTypes>::construct_masters(p)) {}

    Mesh& operator=(const Mesh&) = delete;
    Mesh(Mesh&)                  = delete;
    Mesh& operator=(Mesh&&) = default;
    Mesh(Mesh&&)            = default;

    std::string GetMeshName() { return this->mesh_name; }
    void SetMeshName(const std::string& mesh_name) { this->mesh_name = mesh_name; }

    MasterElementTypes& GetMasters() { return this->masters; }
    void SetMasters() { this->masters = master_maker<MasterElementTypes>::construct_masters(p); }

    uint GetNumberElements() { return this->elements.size(); }
    uint GetNumberInterfaces() { return this->interfaces.size(); }
    uint GetNumberBoundaries() { return this->boundaries.size(); }
    uint GetNumberDistributedBoundaries() { return this->distributed_boundaries.size(); }

    template <typename ElementType, typename... Args>
    void CreateElement(const uint ID, Args&&... args);
    template <typename InterfaceType, typename... Args>
    void CreateInterface(Args&&... args);
    template <typename BoundaryType, typename... Args>
    void CreateBoundary(Args&&... args);
    template <typename DistributedBoundaryType, typename... Args>
    void CreateDistributedBoundary(Args&&... args);

    template <typename F>
    void CallForEachElement(const F& f);
    template <typename F>
    void CallForEachInterface(const F& f);
    template <typename F>
    void CallForEachBoundary(const F& f);
    template <typename F>
    void CallForEachDistributedBoundary(const F& f);

    template <typename ElementType, typename F>
    void CallForEachElementOfType(const F& f);
    template <typename InterfaceType, typename F>
    void CallForEachInterfaceOfType(const F& f);
    template <typename BoundaryType, typename F>
    void CallForEachBoundaryOfType(const F& f);
    template <typename DistributedBoundaryType, typename F>
    void CallForEachDistributedBoundaryOfType(const F& f);

#ifdef HAS_HPX
    template <typename Archive>
    void serialize(Archive& ar, unsigned) {
        // clang-format off
        ar  & mesh_name
            & p;
        // clang-format on
        Utilities::for_each_in_tuple(elements.data, [&ar](auto& element_map) { ar& element_map; });
    }
#endif
};

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename ElementType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CreateElement(const uint ID, Args&&... args) {
    using MasterType = typename ElementType::ElementMasterType;

    MasterType& master_elt = std::get<Utilities::index<MasterType, MasterElementTypes>::value>(this->masters);

    this->elements.template emplace<ElementType>(ID, ElementType(ID, master_elt, std::forward<Args>(args)...));
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename InterfaceType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CreateInterface(Args&&... args) {
    this->interfaces.template emplace_back<InterfaceType>(std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename BoundaryType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CreateBoundary(Args&&... args) {
    this->boundaries.template emplace_back<BoundaryType>(std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename DistributedBoundaryType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CreateDistributedBoundary(Args&&... args) {
    this->distributed_boundaries.template emplace_back<DistributedBoundaryType>(std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CallForEachElement(const F& f) {
    Utilities::for_each_in_tuple(this->elements.data, [&f](auto& element_map) {
        std::for_each(element_map.begin(), element_map.end(), [&f](auto& pair) { f(pair.second); });
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CallForEachInterface(const F& f) {
    Utilities::for_each_in_tuple(this->interfaces.data, [&f](auto& interface_vector) {
        std::for_each(interface_vector.begin(), interface_vector.end(), f);
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CallForEachBoundary(const F& f) {
    Utilities::for_each_in_tuple(this->boundaries.data, [&f](auto& boundary_vector) {
        std::for_each(boundary_vector.begin(), boundary_vector.end(), f);
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CallForEachDistributedBoundary(const F& f) {
    Utilities::for_each_in_tuple(this->distributed_boundaries.data, [&f](auto& distributed_boundaries_vector) {
        std::for_each(distributed_boundaries_vector.begin(), distributed_boundaries_vector.end(), f);
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename ElementType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CallForEachElementOfType(const F& f) {
    auto& el_container =
        std::get<Utilities::index<ElementType, typename ElementContainer::TupleType>::value>(this->elements.data);

    std::for_each(el_container.begin(), el_container.end(), f);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename InterfaceType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CallForEachInterfaceOfType(const F& f) {
    auto& intface_container =
        std::get<Utilities::index<InterfaceType, typename InterfaceContainer::TupleType>::value>(this->interfaces.data);

    std::for_each(intface_container.begin(), intface_container.end(), f);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename BoundaryType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CallForEachBoundaryOfType(const F& f) {
    auto& bound_container =
        std::get<Utilities::index<BoundaryType, typename BoundaryContainer::TupleType>::value>(this->boundaries.data);

    std::for_each(bound_container.begin(), bound_container.end(), f);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries>
template <typename DistributedBoundaryType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>>::CallForEachDistributedBoundaryOfType(const F& f) {
    auto& dbound_container =
        std::get<Utilities::index<DistributedBoundaryType, typename DistributedBoundaryContainer::TupleType>::value>(
            this->distributed_boundaries.data);

    std::for_each(dbound_container.begin(), dbound_container.end(), f);
}
}

#endif