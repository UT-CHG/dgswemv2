#ifndef MESH_HPP
#define MESH_HPP

#include "general_definitions.hpp"
#include "utilities/heterogeneous_containers.hpp"
#include "utilities/is_vectorized.hpp"
#include "mesh_utilities.hpp"

#include "geometry/element_container.hpp"

namespace Geometry {
// Since elements types already come in a tuple. We can use specialization
// to get easy access to the parameter packs for the element and edge types.
template <typename ElementTypeTuple,
          typename InterfaceTypeTuple,
          typename BoundaryTypeTuple,
          typename DistributedBoundaryTypeTuple,
          typename DataSoAType
 >
class Mesh;

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
class Mesh<std::tuple<Elements...>,
           std::tuple<Interfaces...>,
           std::tuple<Boundaries...>,
           std::tuple<DistributedBoundaries...>,
           DataSoAType> {
  public:
    using ElementContainers            = std::tuple<ElementContainer<Elements,DataSoAType>...>;
    using InterfaceContainer           = Utilities::HeterogeneousVector<Interfaces...>;
    using BoundaryContainer            = Utilities::HeterogeneousVector<Boundaries...>;
    using DistributedBoundaryContainer = Utilities::HeterogeneousVector<DistributedBoundaries...>;

  private:
    uint p;

    ElementContainers elements;
    InterfaceContainer interfaces;
    BoundaryContainer boundaries;
    DistributedBoundaryContainer distributed_boundaries;

    std::string mesh_name;

  public:
    Mesh() = default;
    Mesh(const uint p) : p(p), elements(container_maker<ElementContainers>::construct_containers(p)) {}

    Mesh& operator=(const Mesh&) = delete;
    Mesh(Mesh&)                  = delete;
    Mesh& operator=(Mesh&&) = default;
    Mesh(Mesh&&)            = default;

    std::string GetMeshName() { return this->mesh_name; }
    void SetMeshName(const std::string& mesh_name) { this->mesh_name = mesh_name; }

    uint GetNumberElements() {
        size_t sz = 0U;
        Utilities::for_each_in_tuple(elements, [&sz](auto& elt_container) {
                sz += elt_container.size();
            });
        return sz;
    }
    uint GetNumberInterfaces() { return this->interfaces.size(); }
    uint GetNumberBoundaries() { return this->boundaries.size(); }
    uint GetNumberDistributedBoundaries() { return this->distributed_boundaries.size(); }

    void reserve(uint nstages, const std::array<uint,sizeof...(Elements)>& n_elements) {
        Utilities::for_each_in_tuple(elements, [nstages, &n_elements](auto& element_container) {
                using ContainerType = typename std::remove_reference<decltype(element_container)>::type;
                uint index_of_type = Utilities::index<ContainerType,ElementContainers>::value;
                element_container.reserve(nstages, n_elements[index_of_type]);
            });
    }

    void finalize_initialization() {
        Utilities::for_each_in_tuple(elements, [](auto& element_container) {
                element_container.finalize_initialization();
            });
    }

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

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename ElementType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CreateElement(const uint ID, Args&&... args) {
    using ElementContainerType = ElementContainer<ElementType, DataSoAType>;

    ElementContainerType& elt_container = std::get<Utilities::index<ElementContainerType, ElementContainers>::value>(this->elements);

    elt_container.CreateElement(ID, std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename InterfaceType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CreateInterface(Args&&... args) {
    this->interfaces.template emplace_back<InterfaceType>(std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename BoundaryType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CreateBoundary(Args&&... args) {
    this->boundaries.template emplace_back<BoundaryType>(std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename DistributedBoundaryType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CreateDistributedBoundary(Args&&... args) {
    this->distributed_boundaries.template emplace_back<DistributedBoundaryType>(std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CallForEachElement(const F& f) {
    Utilities::for_each_in_tuple(this->elements, [&f](auto& element_container) {
            element_container.CallForEachElement(f, Utilities::is_vectorized<F>{});
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CallForEachInterface(const F& f) {
    Utilities::for_each_in_tuple(this->interfaces.data, [&f](auto& interface_vector) {
        std::for_each(interface_vector.begin(), interface_vector.end(), f);
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CallForEachBoundary(const F& f) {
    Utilities::for_each_in_tuple(this->boundaries.data, [&f](auto& boundary_vector) {
        std::for_each(boundary_vector.begin(), boundary_vector.end(), f);
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CallForEachDistributedBoundary(const F& f) {
    Utilities::for_each_in_tuple(this->distributed_boundaries.data, [&f](auto& distributed_boundaries_vector) {
        std::for_each(distributed_boundaries_vector.begin(), distributed_boundaries_vector.end(), f);
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename ElementType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CallForEachElementOfType(const F& f) {
    auto& elt_container =
        std::get<Utilities::index<ElementType, ElementContainers>::value>(this->elements.data);

    elt_container.CallForEachElement(f, Utilities::is_vectorized<F>{});
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename InterfaceType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CallForEachInterfaceOfType(const F& f) {
    auto& intface_container =
        std::get<Utilities::index<InterfaceType, typename InterfaceContainer::TupleType>::value>(this->interfaces.data);

    std::for_each(intface_container.begin(), intface_container.end(), f);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename BoundaryType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CallForEachBoundaryOfType(const F& f) {
    auto& bound_container =
        std::get<Utilities::index<BoundaryType, typename BoundaryContainer::TupleType>::value>(this->boundaries.data);

    std::for_each(bound_container.begin(), bound_container.end(), f);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType>
template <typename DistributedBoundaryType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType>::CallForEachDistributedBoundaryOfType(const F& f) {
    auto& dbound_container =
        std::get<Utilities::index<DistributedBoundaryType, typename DistributedBoundaryContainer::TupleType>::value>(
            this->distributed_boundaries.data);

    std::for_each(dbound_container.begin(), dbound_container.end(), f);
}
}

#endif