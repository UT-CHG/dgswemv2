#ifndef MESH_HPP
#define MESH_HPP

#include "general_definitions.hpp"
#include "utilities/heterogeneous_containers.hpp"
#include "utilities/is_vectorized.hpp"
#include "mesh_utilities.hpp"

#include "geometry/element_container.hpp"
#include "geometry/interface_container.hpp"
#include "geometry/boundary_container.hpp"

namespace Geometry {
// Since elements types already come in a tuple. We can use specialization
// to get easy access to the parameter packs for the element and edge types.
template <typename ElementTypeTuple,
          typename InterfaceTypeTuple,
          typename BoundaryTypeTuple,
          typename DistributedBoundaryTypeTuple,
          typename DataSoAType,
          typename InterfaceDataSoAType
 >
class Mesh;

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
class Mesh<std::tuple<Elements...>,
           std::tuple<Interfaces...>,
           std::tuple<Boundaries...>,
           std::tuple<DistributedBoundaries...>,
           DataSoAType,
           InterfaceDataSoAType> {
  public:
    using ElementContainers             = std::tuple<ElementContainer<Elements,DataSoAType>...>;
    using InterfaceContainers           = std::tuple<InterfaceContainer<Interfaces,InterfaceDataSoAType,ElementContainers>...>;
    using BoundaryContainers            = std::tuple<BoundaryContainer<Boundaries>...>;
    using DistributedBoundaryContainers = std::tuple<BoundaryContainer<DistributedBoundaries>...>;

  private:
    uint p;

    ElementContainers elements;
    InterfaceContainers interfaces;
    BoundaryContainers boundaries;
    DistributedBoundaryContainers distributed_boundaries;

    std::string mesh_name;

  public:
    Mesh() = default;
    Mesh(const uint p) : p(p), elements(container_maker<ElementContainers>::construct_containers(p)) {
        this->interfaces = container_maker<InterfaceContainers>::construct_containers(p, this->elements);
    }
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

    uint GetNumberInterfaces() {
        size_t sz = 0U;
        Utilities::for_each_in_tuple(interfaces, [&sz](auto& intfc_container) {
                sz += intfc_container.size();
            });
        return sz;
    }

    uint GetNumberBoundaries() {
        size_t sz = 0U;
        Utilities::for_each_in_tuple(boundaries, [&sz](auto& bdry_container) {
                sz += bdry_container.size();
            });
        return sz;
    }

    uint GetNumberDistributedBoundaries() {
      size_t sz = 0U;
      Utilities::for_each_in_tuple(distributed_boundaries, [&sz](auto& distr_bdry_container) {
	  sz += distr_bdry_container.size();
	});
      return sz;
    }

    void reserve(uint nstages, const std::array<uint,sizeof...(Elements)>& n_elements) {
        uint max_edge_gp = 0u;
        Utilities::for_each_in_tuple(interfaces, [this, &max_edge_gp](auto& interface_container) {
                using intfc_container_t = typename std::decay<decltype(interface_container)>::type;
                using integration_t = typename intfc_container_t::InterfaceType::InterfaceIntegrationType;

                integration_t rule;
                max_edge_gp = std::max(max_edge_gp, rule.GetNumGP(2 * p + 1));
            });

        Utilities::for_each_in_tuple(elements, [nstages, &n_elements, max_edge_gp](auto& element_container) {
                using ContainerType = typename std::remove_reference<decltype(element_container)>::type;
                uint index_of_type = Utilities::index<ContainerType,ElementContainers>::value;
                element_container.reserve(nstages, n_elements[index_of_type], max_edge_gp);
            });
    }

    template <typename InterfaceType>
    void reserve_interfaces(size_t n_interfaces) {
        using InterfaceContainerType = InterfaceContainer<InterfaceType, InterfaceDataSoAType, ElementContainers>;

        auto& intface_container =
            std::get<Utilities::index<InterfaceContainerType, InterfaceContainers>::value>(this->interfaces);

        intface_container.SetElementData(this->elements);
        intface_container.reserve(n_interfaces);
    }

    template <typename BoundaryType>
    void reserve_boundaries(size_t n_boundaries) {
        using BoundaryContainerType = BoundaryContainer<BoundaryType>;

        auto& bdry_container =
            std::get<Utilities::index<BoundaryContainerType, BoundaryContainers>::value>(this->boundaries);

        bdry_container.reserve(n_boundaries);
    }

  template <typename DBType>
  void reserve_distributed_boundaries(size_t n_distributed_boundaries) {
    using DistributedBoundaryContainerType = BoundaryContainer<DBType>;
    auto& dbdry_container =
      std::get<Utilities::index<DistributedBoundaryContainerType, DistributedBoundaryContainers>::value>(this->distributed_boundaries);

    dbdry_container.reserve(n_distributed_boundaries);
  }

    void finalize_initialization() {
        Utilities::for_each_in_tuple(elements, [](auto& element_container) {
                element_container.finalize_initialization();
            });

        Utilities::for_each_in_tuple(interfaces, [](auto& interface_container) {
                interface_container.finalize_initialization();
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

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename ElementType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CreateElement(const uint ID, Args&&... args) {
    using ElementContainerType = ElementContainer<ElementType, DataSoAType>;

    ElementContainerType& elt_container = std::get<Utilities::index<ElementContainerType, ElementContainers>::value>(this->elements);

    elt_container.CreateElement(ID, std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename InterfaceType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CreateInterface(Args&&... args) {
    using InterfaceContainerType = InterfaceContainer<InterfaceType, InterfaceDataSoAType, ElementContainers>;

    InterfaceContainerType& intfc_container = std::get<Utilities::index<InterfaceContainerType, InterfaceContainers>::value>(this->interfaces);

    intfc_container.CreateInterface(std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename BoundaryType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CreateBoundary(Args&&... args) {
    using BoundaryContainerType = BoundaryContainer<BoundaryType>;

    BoundaryContainerType& bdry_container = std::get<Utilities::index<BoundaryContainerType, BoundaryContainers>::value>(this->boundaries);

    bdry_container.CreateBoundary(std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename DistributedBoundaryType, typename... Args>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CreateDistributedBoundary(Args&&... args) {
    using DistributedBoundaryContainerType = BoundaryContainer<DistributedBoundaryType>;

    DistributedBoundaryContainerType& distr_bdry_container =
      std::get<Utilities::index<DistributedBoundaryContainerType, DistributedBoundaryContainers>::value>(distributed_boundaries);
    distr_bdry_container.CreateBoundary(std::forward<Args>(args)...);
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CallForEachElement(const F& f) {
    Utilities::for_each_in_tuple(this->elements, [&f](auto& element_container) {
            using element_t = typename std::decay<decltype(element_container)>::type::ElementType;

            element_container.CallForEachElement(f, Utilities::is_vectorized<F, element_t>{});
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CallForEachInterface(const F& f) {
    Utilities::for_each_in_tuple(this->interfaces, [&f](auto& interface_container) {
            using interface_t = typename std::decay<decltype(interface_container)>::type::InterfaceType;

            interface_container.CallForEachInterface(f, Utilities::is_vectorized<F,interface_t>{});
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CallForEachBoundary(const F& f) {
    Utilities::for_each_in_tuple(this->boundaries, [&f](auto& boundary_container) {
            boundary_container.CallForEachBoundary(f, Utilities::is_vectorized<F>{});
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CallForEachDistributedBoundary(const F& f) {
    Utilities::for_each_in_tuple(this->distributed_boundaries, [&f](auto& distributed_boundary_container) {
	distributed_boundary_container.CallForEachBoundary(f, Utilities::is_vectorized<F>{});
    });
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename ElementType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CallForEachElementOfType(const F& f) {
    auto& elt_container =
        std::get<Utilities::index<ElementType, ElementContainers>::value>(this->elements);

    using element_t = typename std::decay<decltype(elt_container)>::type::ElementType;

    elt_container.CallForEachElement(f, Utilities::is_vectorized<F, element_t>{});
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename InterfaceType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CallForEachInterfaceOfType(const F& f) {
    auto& intface_container =
        std::get<Utilities::index<InterfaceType, InterfaceContainers>::value>(this->interfaces);

    intface_container.CallForEachInterface(f, Utilities::is_vectorized<F>{});
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename BoundaryType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CallForEachBoundaryOfType(const F& f) {
    auto& bdry_container =
        std::get<Utilities::index<BoundaryType, BoundaryContainers>::value>(this->boundaries);

    bdry_container.CallForEachBoundary(f, Utilities::is_vectorized<F>{});
}

template <typename... Elements, typename... Interfaces, typename... Boundaries, typename... DistributedBoundaries, typename DataSoAType, typename InterfaceDataSoAType>
template <typename DistributedBoundaryType, typename F>
void Mesh<std::tuple<Elements...>,
          std::tuple<Interfaces...>,
          std::tuple<Boundaries...>,
          std::tuple<DistributedBoundaries...>,
          DataSoAType,
          InterfaceDataSoAType>::CallForEachDistributedBoundaryOfType(const F& f) {
    auto& dbound_container =
        std::get<Utilities::index<DistributedBoundaryType, DistributedBoundaryContainers>::value>(
            this->distributed_boundaries);

    dbound_container.CallForEachBoundry(f, Utilities::is_vectorized<F>{});
}
}

#endif