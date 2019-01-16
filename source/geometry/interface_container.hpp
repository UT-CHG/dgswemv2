#ifndef INTERFACE_CONTAINER_HPP
#define INTERFACE_CONTAINER_HPP

#include "geometry/interface.hpp"
#include "geometry/interface_soa.hpp"

namespace Geometry {

template <typename Interface, typename InterfaceDataSoAType, typename ElementContainers>
class InterfaceContainer;

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType,
          typename InterfaceDataSoAType, typename... ElementContainer>
class InterfaceContainer<Interface<dimension, IntegrationType, DataType, SpecializationType>,
                         InterfaceDataSoAType,
                         std::tuple<ElementContainer...> > {
public:
    using ElementContainers = std::tuple<ElementContainer...>;
    using InterfaceType = Interface<dimension, IntegrationType, DataType, SpecializationType>;

    InterfaceContainer() = default;
    InterfaceContainer(uint p, ElementContainers& elt_data_) : capacity(0UL), size_(0UL),  interface_data(p, elt_data_) {}

    void reserve(uint ninterfaces) {
        this->capacity = ninterfaces;

        interface_accessors.reserve(ninterfaces);

        interface_data.reserve(ninterfaces);
    }

    void finalize_initialization() {
        interface_data.finalize_initialization();

        for ( uint index = 0; index < interface_accessors.size(); ++index ) {
            interface_data.SetNormal(index, interface_accessors[index].surface_normal_in);
        }
    }

    void SetElementData(ElementContainers& element_data_) {
        this->interface_data.SetElementData(element_data_);
    }

    template <typename... Args>
    void CreateInterface(Args&&... args) {
        if ( size_ >= capacity ) {
            throw std::runtime_error{"Not enough interfaces reserved. Emplacing new interfaces will cause accessors to be invalidated\n"};
        }

        interface_accessors.emplace_back(std::move(interface_data.at(size_++,
                                                                     std::forward<Args>(args)...)));
    }

    template <typename F>
    void CallForEachInterface(const F& f, std::false_type /*is not vectorized*/) {
        std::for_each(interface_accessors.begin(), interface_accessors.end(), [&f](InterfaceType& intfc) {
                f(intfc);
            });
    }

    template <typename F>
    void CallForEachInterface(const F& f, std::true_type /*is vectorized*/) {
        f(interface_data);
    }

    size_t size() const { return size_; }
private:
    constexpr static size_t n_element_types = sizeof...(ElementContainer);
    size_t capacity;
    size_t size_;

    AlignedVector<InterfaceType> interface_accessors;
    InterfaceSoA<InterfaceType, InterfaceDataSoAType, ElementContainers> interface_data;
};


}

#endif