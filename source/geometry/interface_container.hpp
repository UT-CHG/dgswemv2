#ifndef INTERFACE_CONTAINER_HPP
#define INTERFACE_CONTAINER_HPP

#include "geometry/interface.hpp"
#include "geometry/interface_soa.hpp"

namespace Geometry {

template <typename Interface, typename ElementContainers>
class InterfaceContainer;

template <uint dimension, typename IntegrationType, typename DataType, typename SpecializationType,
          typename... ElementContainer>
class InterfaceContainer<Interface<dimension, IntegrationType, DataType, SpecializationType>,
                         std::tuple<ElementContainer...> > {
public:
    using ElementContainers = std::tuple<ElementContainer...>;
    using InterfaceType = Interface<dimension, IntegrationType, DataType, SpecializationType>;

    InterfaceContainer() = default;
    InterfaceContainer(ElementContainers& elt_data_) : capacity(0UL), size_(0UL),  interface_data(elt_data_) {}

    void reserve(uint ninterfaces) {
        this->capacity = ninterfaces;

        interface_accessors.reserve(ninterfaces);

//        interface_data.reserve(ninterfaces, ngp);
    }

    template <typename... Args>
    void CreateInterface(Args&&... args) {
        if ( size_++ >= capacity ) {
            throw std::runtime_error{"Not enough interfaces reserved. Emplacing new interfaces will cause accessors to be invalidated\n"};
        }

        interface_accessors.emplace_back(std::forward<Args>(args)...);
    }

    template <typename F>
    void CallForEachInterface(const F& f, std::false_type /*is not vectorized*/) {
        std::for_each(interface_accessors.begin(), interface_accessors.end(), [&f](InterfaceType& intfc) {
                f(intfc);
            });
    }

    size_t size() const { return size_; }
private:


    AlignedVector<InterfaceType> interface_accessors;
    InterfaceSoA<InterfaceType, ElementContainers> interface_data;

    size_t capacity;
    size_t size_;
};


}

#endif